"""Recompute clinical and allele-level Letter results from the supplementary tables."""
import argparse
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
S1, S2, S3 = 'S1 Cohort', 'S2 Variants', 'S3 Clinical associations'
ENDPOINTS = [('clinical',8),('first',9),('hematologic',11),('mixed',21)]
PHENOTYPES = ['Hematologic predominant','Mixed: Hematologic + Neurologic','Neurologic predominant',
    'Mixed: Hematologic + Neurologic + Immune/inflammatory','Mixed: Hematologic + Immune/inflammatory',
    'Immune/inflammatory predominant','Mixed: Neurologic + Immune/inflammatory',
    'Unaffected/limited phenotype','Uncertain/insufficient']
DIAGNOSES = ['MDS/RCC','Cytopenia-only / transient','BMF/AA','AML','MDS/AML progression',
             'MPN/PMF','MDS/MPN','No hematologic disease reported','Other/uncertain','Unknown']


def age_group(age):
    return 'Unknown' if pd.isna(age) else ('Pediatric' if age<=18 else 'Adult')


def effect(age,position):
    rho,p=stats.spearmanr(age,position)
    return {'n':len(age),'rho':float(rho),'p':float(p)}


def save_json(path,data):
    def encode(item):
        if isinstance(item,np.generic):return item.item()
        raise TypeError(type(item))
    path.write_text(json.dumps(data,default=encode,allow_nan=False,indent=2,ensure_ascii=False)+'\n')


def exact_sahoo(aa,loss):
    ranks=stats.rankdata(aa);normal=~loss;n=sum(normal);center=n*(len(ranks)+1)/2
    observed=abs(ranks[normal].sum()-center);assignments=list(itertools.combinations(ranks,n))
    extreme=sum(abs(sum(x)-center)>=observed for x in assignments)
    return {'loss':int(sum(loss)),'normal':int(n),'loss_median':float(np.median(aa[loss])),
            'normal_median':float(np.median(aa[normal])),'extreme':int(extreme),
            'assignments':len(assignments),'p':float(extreme/len(assignments))}


def family_results(d,out,patch):
    support=pd.read_csv(ROOT/'data/source_support.csv').rename(columns={'Patient_Master_ID':'patient'})
    d=d.merge(support,on='patient',validate='one_to_one');known=d.dropna(subset=['family'])
    units=known.groupby('family').agg(n=('patient','size'),age=('age','median'),position=('position','median'))
    proteins=d.groupby('protein').agg(age=('age','median'),position=('position','median'))
    sets=[('All',d),('Exclude source VUS-only',d[d.source_classification_support!='Source VUS/uncertain without pathogenic label']),
          ('Exclude competing therapy',d[~d.competing_etiology_source_flag]),('Source pathogenic label',d[d.source_classification_support.isin(['Source pathogenic/likely pathogenic label','Source mixed/changed VUS and pathogenic labels'])]),
          ('Known family',known),('Unresolved family',d[d.family.isna()]),('Family medians',units),('Protein medians',proteins)]
    result=[]
    for i,(label,x) in enumerate(sets,95):
        r=effect(x.age,x.position);result.append(dict(label=label,**r))
        for col,k in zip('CDE',['n','rho','p']):patch[f'{col}{i}']=r[k]
    rng=np.random.default_rng(20261004);idx=rng.integers(len(units),size=(10000,len(units)))
    aa=stats.rankdata(units.position.to_numpy()[idx],axis=1);ages=stats.rankdata(units.age.to_numpy()[idx],axis=1)
    aa-=aa.mean(axis=1,keepdims=True);ages-=ages.mean(axis=1,keepdims=True)
    denominator=np.linalg.norm(aa,axis=1)*np.linalg.norm(ages,axis=1)
    rr=np.sum(aa*ages,axis=1)/np.where(denominator==0,np.nan,denominator)
    lo,hi=np.nanquantile(rr,[.025,.975]);patch.update({'F101':lo,'G101':hi})
    result[6].update(ci_low=lo,ci_high=hi,valid_replicates=sum(np.isfinite(rr)))
    patch['A93']=(f'Individual analyses retain related carriers. Among {len(d)} eligible individuals, {len(known)} have documented family or de novo units and {len(d)-len(known)} lack a resolved unit. The family-median analysis gives each of the {len(units)} documented units one observation, using median age and position and 10,000 unit-level bootstrap samples (seed 20261004). Protein-label grouping is a separate sensitivity. Source support labels follow published classifications. P values are nominal.')
    pd.DataFrame(result).to_csv(out/'family.csv',index=False);units.to_csv(out/'family_units.csv')
    return result


def run(workbook, out, policy):
    out.mkdir(parents=True,exist_ok=True)
    wb=load_workbook(workbook,read_only=True,data_only=True);sheets={s.title:list(s.values) for s in wb}
    s1,s2,s3=(sheets[k] for k in [S1,S2,S3])
    rows1=s1[5:183]
    patients={r[0]:r for r in rows1};assert len(patients)==178 and None not in patients
    rows2=s2
    links=[r for r in rows2 if r[0] in patients and isinstance(r[1],str) and r[1].startswith('INT') and 'modifier' not in str(r[5]).lower()]
    inventory={r[1]:r for r in rows2 if r[0] in ['Germline','Somatic'] and r[1]}
    events=[r for r in rows2 if isinstance(r[0],str) and r[0].startswith('SGR') and r[1] in patients]
    assert len({r[0] for r in events})==len(events)
    bypatient=defaultdict(list)
    for r in links:bypatient[r[0]].append(r)
    assert set(bypatient)==set(patients)
    patch={S1:{},S3:{}};p3=patch[S3];p1=patch[S1]
    summary={'cohort':len(patients),'affected':Counter(r[5] for r in patients.values()),'link_count':len(links),
             'multiple':sum(len(v)>1 for v in bypatient.values()),'phenotypes':Counter(r[7] for r in patients.values()),
             'diagnoses':Counter(r[6] for r in patients.values())}
    for names, field in [(PHENOTYPES, 'phenotypes'), (DIAGNOSES, 'diagnoses')]:
        diff = set(summary[field]) - set(names)
        if diff:
            raise ValueError(f'New categories need explicit mapping: {diff}')
    pri_counts=Counter(r[1] for r in links)
    primary=[r for r in inventory.values() if r[0]=='Germline' and r[2]=='Yes']
    assert {r[1] for r in primary}==set(pri_counts)
    summary['germline']={'unique':len(primary),'types':Counter(r[6] for r in primary),'carriers_per_allele':Counter(pri_counts.values())}
    gnomad=pd.DataFrame([dict(allele=r[1],protein=r[4],position=r[5],af=r[12],status=r[8]) for r in primary])
    positive=gnomad[gnomad.status.str.startswith('Observed')];negative=gnomad[gnomad.status.str.startswith('Not observed')]
    assert len(positive)+len(negative)==len(gnomad) and positive.af.notna().all() and negative.af.isna().all()
    gnomad.to_csv(out/'gnomad.csv',index=False)
    g={'observed':len(positive),'unobserved':len(negative),'median_observed':positive.position.median(),
       'median_unobserved':negative.position.median(),'p':stats.mannwhitneyu(positive.position,negative.position).pvalue}
    summary['gnomad']=g;p3.update({'C29':len(gnomad),'E29':g['p'],'G29':f"{g['observed']}/{g['unobserved']}; medians {g['median_observed']:g}/{g['median_unobserved']:g}"})
    event_counts=Counter(r[1] for r in events)
    exacts=Counter('c.2699A>G' if r[0]=='SGR049' else r[3] for r in events if r[3] and str(r[3]).startswith('c.'))
    somatic=[r for r in inventory.values() if r[0]=='Somatic' and str(r[3]).startswith('c.')]
    if set(exacts)!=set(r[3] for r in somatic):raise ValueError('Exact cDNA inventory/event mismatch: '+str(set(exacts)^set(r[3] for r in somatic)))
    summary['acquired']={'events':len(events),'individuals':len(event_counts),'exact_events':sum(exacts.values()),'exact_alleles':len(exacts),
                         'observations_per_exact_allele':Counter(exacts.values()),'exact_types':Counter(r[6] for r in somatic),
                         'phase_status':Counter(r[6] for r in events),'Arg589Ter':[(r[0],r[1],r[2],r[6]) for r in events if r[3]=='c.1765C>T']}
    clonal=pd.DataFrame([dict(patient=k,chr7=r[16]=='Reported',upd=r[17]=='Reported',intra=r[18]>0,
                              clinical_group=age_group(r[8]),event_count=r[18]) for k,r in patients.items()])
    mismatch=[r.patient for r in clonal.itertuples() if r.event_count!=event_counts[r.patient]]
    if mismatch:raise ValueError('S1 intragenic vs S2 observations: '+str(mismatch))
    clonal['routes']=clonal[['chr7','upd','intra']].sum(axis=1);clonal['any_recorded']=clonal.routes>0
    summary['clonal']={k:int(clonal[k].sum()) for k in ['chr7','upd','intra','any_recorded']}
    summary['clonal'].update(multiple=int((clonal.routes>=2).sum()),all_three=int((clonal.routes==3).sum()))
    for i,group in enumerate(['Adult','Pediatric','Unknown'],35):
        sample=clonal.loc[clonal.clinical_group==group]
        for col,key in zip('BCDE',['chr7','upd','intra','any_recorded']):p3[f'{col}{i}']=int(sample[key].sum())
    clonal.to_csv(out/'clonal.csv',index=False)
    ages=[];reasons=[]
    for pt,r in patients.items():
        for endpoint,col in ENDPOINTS:
            reason='';ll=bypatient[pt]
            if r[5]!='Yes':reason='Not coded affected'
            elif len(ll)!=1:reason='Multiple primary alleles'
            elif pd.isna(r[col]):reason='Age unavailable'
            elif endpoint in policy.get('exclude_ages',{}).get(pt,[]):reason='Curated exclusion for this endpoint'
            elif endpoint=='mixed' and r[22]!='Yes':reason='Mixed age excluded in adjudicated table'
            elif endpoint=='clinical' and r[12]!='Yes':reason='Clinical age excluded in adjudicated table'
            if reason:reasons.append([pt,endpoint,reason]);continue
            v=inventory[ll[0][1]];age=policy.get('override_ages',{}).get(pt,{}).get(endpoint,r[col])
            ages.append(dict(patient=pt,endpoint=endpoint,age=age,protein=ll[0][2],position=v[5],allele=ll[0][1],
                     variant_type=ll[0][4],source='St Jude' if 'St Jude' in r[1] else 'Other',family=r[24],group=age_group(age)))
    ages=pd.DataFrame(ages);ages.to_csv(out/'age_rows.csv',index=False)
    pd.DataFrame(reasons,columns=['patient','endpoint','reason']).to_csv(out/'age_exclusions.csv',index=False)
    summary['age']={}
    for j,(name,_) in enumerate(ENDPOINTS):
        d=ages[ages.endpoint==name];r=effect(d.age,d.position);ped=d[d.group=='Pediatric'];adult=d[d.group=='Adult']
        nped,nad=len(ped),len(adult);tp=(ped.variant_type!='Missense').sum();ta=(adult.variant_type!='Missense').sum()
        r.update(n_ped=nped,n_adult=nad,med_ped=ped.position.median(),med_adult=adult.position.median(),mw=stats.mannwhitneyu(ped.position,adult.position).pvalue,
                 trunc_ped=int(tp),trunc_adult=int(ta),fisher=stats.fisher_exact([[tp,nped-tp],[ta,nad-ta]]).pvalue,
                 unique_ped=ped.allele.nunique(),unique_adult=adult.allele.nunique());summary['age'][name]=r
        for c,k in zip('BCD',['n','rho','p']):p3[f'{c}{6+j}']=r[k]
        row=15+2*j;p3.update({f'C{row}':len(d),f'E{row}':r['mw'],f'G{row}':f"n={nped}/{nad}; median={r['med_ped']:g}/{r['med_adult']:g}",
                   f'C{row+1}':len(d),f'E{row+1}':r['fisher'],f'G{row+1}':f'truncating={tp}/{nped} vs {ta}/{nad}'})
    clinical=ages[ages.endpoint=='clinical'];summary['source']={}
    for i,source in enumerate(['Other','St Jude'],25):
        d=clinical[clinical.source==source];t=effect(d.age,d.position);summary['source'][source]=t
        for c,k in zip('CDE',['n','rho','p']):p3[f'{c}{i}']=t[k]
    summary['family']=family_results(clinical,out,p3)
    inf=pd.DataFrame([dict(patient=k,positive=r[23]=='Reported',truncating=any(a[4]!='Missense' for a in bypatient[k]))
                       for k,r in patients.items() if r[23] in ('Reported','Not source-coded')])
    li=pd.DataFrame([dict(patient=r[0],allele=r[1],truncating=r[4]!='Missense') for r in links]).merge(inf[['patient','positive']],on='patient')
    al=li.groupby('allele')[['positive','truncating']].max().reset_index()
    for row,(name,d) in enumerate([('individuals',inf),('alleles',al)],23):
        a,b=d[d.truncating],d[~d.truncating];t=int(a.positive.sum());m=int(b.positive.sum());odds,p=stats.fisher_exact([[t,len(a)-t],[m,len(b)-m]])
        summary['inflammation_'+name]={'n':len(d),'positive_trunc':t,'all_trunc':len(a),'positive_missense':m,'all_missense':len(b),'p':p}
        p3.update({f'C{row}':len(d),f'D{row}':float(odds) if np.isfinite(odds) else 'Infinite odds ratio (zero cell)',f'E{row}':p,f'G{row}':f'Truncating {t}/{len(a)} vs missense {m}/{len(b)} source-coded positive; Fisher exact two-sided.'})
        d.to_csv(out/f'inflammation_{name}.csv',index=False)
    broad = inf.copy()
    broad["positive"] = broad.patient.map(
        lambda identifier: "immune/inflammatory" in patients[identifier][7].lower()
    )
    truncating = broad[broad.truncating]
    missense = broad[~broad.truncating]
    positive_trunc = int(truncating.positive.sum())
    positive_missense = int(missense.positive.sum())
    odds, p_value = stats.fisher_exact([
        [positive_trunc, len(truncating) - positive_trunc],
        [positive_missense, len(missense) - positive_missense],
    ])
    summary["inflammation_broad"] = {
        "n": len(broad), "positive_trunc": positive_trunc, "all_trunc": len(truncating),
        "positive_missense": positive_missense, "all_missense": len(missense),
        "odds_ratio": odds, "p": p_value,
    }
    broad.to_csv(out / "inflammation_broad.csv", index=False)
    sahoo_rows = []
    for section, cohort in [("D.", "A"), ("E.", "B")]:
        start, = [i for i, row in enumerate(s3) if row and str(row[0]).startswith(section)]
        header = start + 2
        assert s3[header][:3] == ("Patient ID (S1)", "Sahoo UPN", "Source cohort")
        for row in s3[header + 1:]:
            if len(row) <= 2 or row[2] != cohort:
                break
            if row[8] == "Yes":
                sahoo_rows.append(dict(patient=row[0], upn=row[1], cohort=cohort,
                                       baseline_karyotype=row[3], position=row[7],
                                       loss=row[3] != "NK"))
    sahoo = pd.DataFrame(sahoo_rows)
    assert set(sahoo.baseline_karyotype)<=set(['NK','-7','-7/complex'])
    sahoo.to_csv(out/'sahoo_baseline.csv',index=False);summary['sahoo']=[]
    for i,(cohort,exclude) in enumerate(itertools.product(['A','B'],[False,True]),107):
        d=sahoo[sahoo.cohort==cohort]
        if exclude:d=d[d.position!=986]
        t=exact_sahoo(d.position.to_numpy(),d.loss.to_numpy());summary['sahoo'].append(dict(cohort=cohort,exclude_986=exclude,**t))
        for c,k in zip('CDEFGHI',['loss','normal','loss_median','normal_median','extreme','assignments','p']):p3[f'{c}{i}']=t[k]
        if cohort=='A':
            row=i-80;p3.update({f'C{row}':len(d),f'E{row}':t['p'],f'G{row}':f"n loss={t['loss']}; assignments={t['assignments']}; median loss/normal={t['loss_median']:g}/{t['normal_median']:g}"})
    save_json(out/'clinical_summary.json',summary);save_json(out/'table_statistics.json',{sheet:{cell:value for cell,value in cells.items() if isinstance(value,(int,float,np.integer,np.floating)) and not isinstance(value,(bool,np.bool_))} for sheet,cells in patch.items()})
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--workbook',type=Path,required=True)
    parser.add_argument('--output',type=Path,default=ROOT/'results/clinical')
    parser.add_argument('--policy',type=Path,default=ROOT/'data/clinical_policy.json')
    a=parser.parse_args()
    run(a.workbook,a.output,json.loads(a.policy.read_text()))
