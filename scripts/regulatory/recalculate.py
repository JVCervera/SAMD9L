"""Recompute molecular analyses from exported model predictions and the AlphaFold monomer."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import chromatin
import geometry
from read_events import roster

HERE=Path(__file__).resolve().parent


def bh(values):
    p=np.asarray(values,float)
    return stats.false_discovery_control(p,method='bh')


def phase_groups(events):
    valid=events[events.sequence]
    confirmed={'confirmed_cis','confirmed_cis_both'}
    reported=confirmed|{'reported_cis'}
    statuses=[confirmed,reported,reported|{'inferred_cis'},set(valid.phase)]
    labels=['confirmed','confirmed_reported','plus_inferred','all_modelled']
    return {label+str(len(ids)):ids for label,ph in zip(labels,statuses) for ids in [set(valid.loc[valid.phase.isin(ph),'event'])]}


def gene_tests(pairs, groups, out, prefix=''):
    tests,backgrounds=[],[]
    for label, ids in groups.items():
        eligible=pairs[pairs.pair_id.isin(ids)]
        for scope,d in eligible.groupby('scope'):
            for metric, null in [('median_rescue_to_WT',0),('fraction_closer_to_WT',0.5)]:
                family=[]
                for modality,sub in d.groupby('modality'):
                    by_bg=sub.groupby('background')[metric].median().dropna()
                    scores=by_bg-null
                    n_nonzero=int(scores.ne(0).sum())
                    p=float(stats.wilcoxon(scores).pvalue) if n_nonzero else 1.
                    family.append([label,scope,modality,metric,sub.pair_id.nunique(),len(by_bg),n_nonzero,by_bg.median(),null,p])
                    backgrounds.extend([[label,scope,modality,metric,b,v] for b,v in by_bg.items()])
                tests.extend(r+[v] for r,v in zip(family,bh([r[-1] for r in family])))
    cols=['subset','scope','modality','endpoint','n_events_with_scores','n_backgrounds','n_nonzero','median','null','P','q_BH4']
    pd.DataFrame(tests,columns=cols).to_csv(out/(prefix+'gene_tests.csv'),index=False)
    pd.DataFrame(backgrounds,columns=['subset','scope','modality','endpoint','background','value']).to_csv(out/(prefix+'gene_backgrounds.csv'),index=False)


def matched_rna(analysis_dir,groups,out,prefix=''):
    src=analysis_dir/'portable_inputs'
    frozen=pd.read_csv(src/'pair_control_metrics.csv.gz')
    assigned=pd.read_csv(src/'control_assignments.csv')
    all_ids=set.union(*groups.values())
    assigned=assigned[assigned.observed_pair_id.isin(all_ids)]
    assert not assigned.duplicated(['observed_pair_id','control_somatic_variant_id']).any()
    assert assigned.groupby('observed_pair_id').size().eq(5).all() and not assigned.observed_somatic_variant_id.eq(assigned.control_somatic_variant_id).any()
    allowed=set(zip(assigned.observed_pair_id,assigned.control_pair_id))
    primary=frozen[(frozen.modality=='RNA expression')&(frozen.selection=='all')&frozen.pair_id.isin(all_ids)]
    assert set(zip(primary.pair_id,primary.control_id))<=allowed
    tests,events,backgrounds=[],[],[]
    for label,ids in groups.items():
        for scope in ['native_all','native_hematopoietic']:
            sub=primary[primary.pair_id.isin(ids)&primary.scope.eq(scope)]
            ev=sub.groupby(['pair_id','background'],as_index=False).agg(advantage=('median_advantage','median'),fraction=('fraction_advantage','median'),observed_reduction=('observed_reduction','median'),control_reduction=('control_reduction','median'),n_controls=('control_id','nunique'))
            by_bg=ev.groupby('background').advantage.median()
            n_pos=int(by_bg.gt(0).sum());n_neg=int(by_bg.lt(0).sum())
            p=stats.binomtest(n_pos,n_pos+n_neg).pvalue if n_pos+n_neg else 1.
            tests.append([label,scope,len(ev),len(by_bg),by_bg.median(),n_pos,n_neg,p,int(ev.n_controls.min())])
            events.append(ev.assign(subset=label,scope=scope)); backgrounds.extend([[label,scope,b,v] for b,v in by_bg.items()])
    cols=['subset','scope','n_events','n_backgrounds','median_advantage','positive_backgrounds','negative_backgrounds','P_sign','min_controls']
    table=pd.DataFrame(tests,columns=cols)
    table['q_BH8']=bh(table.P_sign)
    table.to_csv(out/(prefix+'rna_tests.csv'),index=False)
    pd.concat(events).to_csv(out/(prefix+'rna_events.csv'),index=False)
    pd.DataFrame(backgrounds,columns=['subset','scope','background','advantage']).to_csv(out/(prefix+'rna_backgrounds.csv'),index=False)
    return {'n_assignments':len(assigned),'n_alternatives':assigned.control_pair_id.nunique()}


def isolated(analysis_dir,out):
    scores=pd.read_csv(analysis_dir/'source_data/S12_AG_isolated_scores.csv')
    normalization = pd.read_csv(HERE.parents[1] / 'data/protein_normalizations.csv').set_index('variant_id')
    corrected_positions = scores.variant_id_project.map(normalization.residue)
    corrected = corrected_positions.notna()
    scores.loc[corrected, 'AA_position'] = corrected_positions[corrected].astype(int)
    univ=pd.read_csv(analysis_dir/'source_data/S9_AG_variant_universe.csv')
    exact=univ[univ.Variant_ID.isin(scores.variant_id_project)].dropna(subset=['GRCh38_allele'])
    assert len(exact)==138 and exact.GRCh38_allele.nunique()==136
    assert not scores.duplicated(['dataset','variant_id_project','modality']).any()
    counts=scores.groupby(['dataset','modality']).variant_id_project.nunique().unstack(0)
    assert counts.shape==(9,2) and counts.GERMLINE.eq(88).all() and counts.SOMATIC.eq(50).all()
    shared_ids=set(exact.loc[exact.GRCh38_allele.duplicated(False),'Variant_ID'])
    assert len(shared_ids)==4
    models=[]
    for spec,shared_removed,quadratic in [('all_linear',False,False),('exclude_shared_linear',True,False),('all_quadratic',False,True),('exclude_shared_quadratic',True,True)]:
        subset = scores
        if shared_removed:
            subset = subset.loc[~subset.variant_id_project.isin(shared_ids)]
        family=[]
        for modality,sub in subset.groupby('modality'):
            pos=(sub.AA_position.to_numpy()-1000)/1000
            cov=[np.ones(len(sub)),sub.dataset.ne('GERMLINE'),pos,sub.truncating]
            if quadratic:cov.append(pos**2)
            x=np.column_stack(cov).astype(float);y=sub.heme_median_signed_q.to_numpy()
            inverse=np.linalg.inv(x.T@x);beta=inverse@x.T@y
            residual=y-x@beta; leverage=np.sum(x*(x@inverse),axis=1)
            adjusted=x*(residual/(1-leverage))[:,None]
            cov_hc3=inverse@adjusted.T@adjusted@inverse;se=np.sqrt(cov_hc3[1,1])
            p=2*stats.norm.sf(abs(beta[1]/se));radius=stats.norm.ppf(.975)*se
            family.append([spec,modality,len(sub),beta[1],beta[1]-radius,beta[1]+radius,p])
        models.extend(row+[q] for row,q in zip(family,bh([r[-1] for r in family])))
    reg=pd.DataFrame(models,columns=['specification','modality','N','Adjusted_beta','CI_low','CI_high','P','q_BH9'])
    reg.loc[reg.specification.ne('all_linear')].to_csv(out/'S4J_sensitivity.csv', index=False)
    reg = reg.loc[reg.specification.eq('all_linear')]
    reg.to_csv(out/'origin_regressions.csv', index=False)
    sub=scores[scores.dataset.eq('GERMLINE')&scores.modality.eq('SPLICE_JUNCTIONS')]
    pos=[]
    for c in ['BM_p90_abs_q','PBMC_p90_abs_q','CD34_p90_abs_q','MYELOID_p90_abs_q','OTHER_HEME_p90_abs_q','delta_p90_q']:
        rho,p=stats.spearmanr(sub.AA_position,sub[c]);pos.append([c,len(sub),rho,p])
    pos=pd.DataFrame(pos,columns=['Context','N_primary_alleles','rho','P']);pos['q_BH6']=bh(pos.P);pos.to_csv(out/'junction_position.csv',index=False)
    scores.groupby(['dataset','modality'],as_index=False).agg(n_alleles=('variant_id_project','nunique'),n_heme_higher=('heme_p90_gt_other','sum')).to_csv(out/'tissue_summary.csv',index=False)
    return scores.variant_id_project.nunique()


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--archive',type=Path,required=True,help='Directory containing the analysis inputs')
    cli.add_argument('--tables',type=Path,required=True,help='Build event roster from Supplementary Tables S1/S2')
    cli.add_argument('--output',type=Path,default=HERE/'results')
    args=cli.parse_args(); out=args.output;out.mkdir(parents=True,exist_ok=True)
    analysis=args.archive/'analysis';events=roster(args.tables,args.archive)
    events.to_csv(out/'events.csv',index=False)
    assert events.event.is_unique and (~events.sequence|events.structural).all()
    groups=phase_groups(events); labels=list(groups); eligible=groups[labels[-1]]
    inputs=analysis/'portable_inputs'
    manifest=pd.read_csv(inputs/'phase_manifest.csv'); modelled=manifest[manifest.pair_id.isin(eligible)]
    assert set(modelled.pair_id)==eligible
    pairs=pd.read_csv(analysis/'source_data/Recovered_gene_centred_pair_summaries.csv')
    gene_tests(pairs,{k:groups[k] for k in [labels[0],labels[1],labels[-1]]},out)
    mapping=matched_rna(analysis,groups,out)
    mixed=set(events.loc[events.phase.eq('confirmed_cis_both'),'event'])
    if mixed:
        exclusive={k:ids-mixed for k,ids in groups.items()}
        gene_tests(pairs,{k:exclusive[k] for k in [labels[0],labels[1],labels[-1]]},out,'omit_dual_configuration_')
        matched_rna(analysis,exclusive,out,'omit_dual_configuration_')
    scores=isolated(analysis,out)
    summary={'n_observations':len(events),'patients_all':events.patient.nunique(),'structural_events':int(events.structural.sum()),
             'structural_backgrounds':events.loc[events.structural,'germline'].nunique(),'sequence_events':len(eligible),
             'sequence_patients':modelled.patient_id.nunique(),'sequence_backgrounds':modelled.germline.nunique(),
             'sequence_distinct_acquired_alleles':modelled.somatic_AG_input.nunique(),'phase_subsets':{k:len(v) for k,v in groups.items()},
             'phase_counts_in_sequence':events[events.sequence].phase.value_counts().to_dict(),'isolated_allele_entries':scores,'RNA_mapping':mapping}
    print('Calculating structural reassociation.',flush=True)
    geometry.run(events,analysis,out)
    provisional=events.loc[events.event.eq('SGR056')]
    if len(provisional):
        without=out/'without_provisional_4B'
        without.mkdir(exist_ok=True)
        geometry.run(events[~events.event.isin(provisional.event)],analysis,without)
    print('Calculating chromatin reassociation.',flush=True)
    summary['chromatin_mapping']=chromatin.run(events,analysis,out,bh)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
