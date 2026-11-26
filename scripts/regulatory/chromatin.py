"""Recompute both 48-contrast chromatin families using exported pair-level track summaries."""
import hashlib
import numpy as np
import pandas as pd

METRICS=['median_rescue_to_WT','fraction_closer_to_WT','median_interaction_log2sum','median_abs_interaction_log2sum']


def run(events,analysis_dir,out,bh):
    inputs = analysis_dir/'portable_inputs'
    observed = pd.read_csv(inputs/'observed_pair_level_chromatin_metrics.csv')
    controls = pd.read_csv(inputs/'control_pair_level_chromatin_metrics.csv')
    assigned = pd.read_csv(inputs/'chromatin_control_assignments_270.csv')
    retained = set(events.loc[events.sequence,'event'])
    assigned = assigned[assigned.observed_pair_id.isin(retained)]
    assert set(assigned.observed_pair_id)==retained
    assert assigned.groupby('observed_pair_id').size().eq(5).all()
    audit, results, distributions = [], [], {}
    for family, minimum in [('PRIMARY',1),('EFFECTIVE_MIN3',3)]:
        for control_set in ['FULL','STRICT3']:
            mp = assigned
            if control_set == 'STRICT3':
                mp = mp[mp.match_tier<=1]
                mp = mp.groupby('observed_pair_id').filter(lambda x:len(x)>=3)
            for scope in ['native_hematopoietic','all_hematopoietic']:
                for mod in ['DNASE','CHIP_HISTONE','CHIP_TF']:
                    def select(d):
                        return d[(d.analysis_scope==scope)&(d.output_type==mod)&(d.track_set=='strong_p90')&
                                 (d.mask_id=='PRIMARY_PAIR_UNION')&(d.variant_subset=='all_pairs')]
                    obs = select(observed); alt=select(controls).set_index('pair_id'); evs=[]
                    for _,row in obs[obs.pair_id.isin(retained)].sort_values('pair_id').iterrows():
                        ids=mp.loc[mp.observed_pair_id==row.pair_id,'control_pair_id'].tolist()
                        ids=[i for i in ids if i in alt.index and np.isfinite(alt.loc[i,METRICS].to_numpy(float)).all()]
                        vec=row[METRICS].to_numpy(float)
                        if len(ids)>=minimum and np.isfinite(vec).all():
                            assert len(ids)==len(set(ids))
                            evs.append((row.pair_id,row.background,vec,alt.loc[ids,METRICS].to_numpy(float),ids))
                    keys=sorted(set(e[1] for e in evs))
                    stat=np.median(np.stack([np.median(np.stack([e[2] for e in evs if e[1]==k]),axis=0) for k in keys]),axis=0)
                    key='|'.join([family,control_set,scope,mod])
                    rng=np.random.default_rng(np.random.SeedSequence([20260910,int.from_bytes(hashlib.sha256(key.encode()).digest()[:4],'little')]))
                    draw=np.stack([e[3][rng.integers(0,len(e[3]),size=10000)] for e in evs],axis=1)
                    dist=np.median(np.stack([np.median(draw[:,[i for i,e in enumerate(evs) if e[1]==k],:],axis=1) for k in keys],axis=1),axis=1)
                    for e in evs:audit.append([family,control_set,scope,mod,e[0],e[1],len(e[4]),';'.join(e[4])])
                    for k,metric in enumerate(METRICS):
                        x=dist[:,k];center=np.median(x);pid=key+'|'+metric
                        p=(1+np.count_nonzero(np.abs(x-center)>=abs(stat[k]-center)))/10001
                        distributions[pid]=x
                        results.append([pid,family,control_set,scope,mod,metric,'PRIMARY_PAIR_UNION','strong_p90','all_pairs',len(evs),len(keys),minimum,stat[k],center,np.quantile(x,.025),np.quantile(x,.975),stat[k]-center,p,48])
    columns=['test_id','analysis_family','control_set','analysis_scope','output_type','metric','mask_id','track_set','variant_subset','n_events','n_backgrounds','min_effective_controls','observed_stat','null_median','null_q025','null_q975','effect_vs_null_median','p_empirical_two_sided','family_n_tests']
    d=pd.DataFrame(results,columns=columns); d['q_family_global']=d.groupby('analysis_family').p_empirical_two_sided.transform(bh)
    d.to_csv(out/'chromatin_tests.csv',index=False)
    pd.DataFrame(audit,columns=['family','control_set','scope','modality','event','background','used_controls','control_ids']).to_csv(out/'chromatin_control_audit.csv',index=False)
    np.savez_compressed(out/'chromatin_null.npz',**distributions)
    return {'assignments':len(assigned),'unique_controls':assigned.control_pair_id.nunique()}
