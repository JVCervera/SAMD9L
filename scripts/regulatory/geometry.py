"""Measure reference-monomer distances and rerun the event-slot reassociation."""
from pathlib import Path
import json
import numpy as np
import pandas as pd


def read_atoms(path):
    atoms, ca, plddt = {}, {}, {}
    for line in path.read_text().splitlines():
        if not line.startswith('ATOM  ') or line[21] != 'A' or line[16] not in (' ', 'A'):
            continue
        name, position = line[12:16].strip(), int(line[22:26])
        if line[76:78].strip() in ('H','D') or name.startswith('H'):
            continue
        xyz = np.array([float(line[k:k+8]) for k in (30,38,46)])
        atoms.setdefault(position, []).append(xyz)
        if name == 'CA':
            ca[position], plddt[position] = xyz, float(line[60:66])
    return {k:np.array(v) for k,v in atoms.items()}, ca, plddt


def run(events, analysis_dir, out):
    models = analysis_dir / 'model_inputs'
    atoms, ca, plddt = read_atoms(models/'AF-Q8IVG5-F1-model_v6.pdb')
    pae = np.array(json.loads((models/'AF-Q8IVG5-F1-predicted_aligned_error_v6.json').read_text())[0]['predicted_aligned_error'])
    d = events[events.structural].copy()
    g,s = d.germline_pos.to_numpy(int), d.somatic_pos.to_numpy(int)
    n = len(d)
    matrix = {k:np.empty((n,n)) for k in ('ca','heavy','pae')}
    for i,a in enumerate(g):
        for j,b in enumerate(s):
            matrix['ca'][i,j] = np.linalg.norm(ca[a]-ca[b])
            matrix['heavy'][i,j] = np.sqrt(((atoms[a][:,None,:]-atoms[b][None,:,:])**2).sum(axis=2)).min()
            matrix['pae'][i,j] = max(pae[a-1,b-1],pae[b-1,a-1])
    v = {k:np.diag(m) for k,m in matrix.items()}
    conf = np.minimum([plddt[a] for a in g],[plddt[b] for b in s])
    pair = pd.DataFrame({'SGR': d.event,'Patient':d.patient,'Source patient':d.source_patient,
                          'Germline':d.germline_label,'Germline aa':g,'Second hit':d.somatic_label,'Second-hit aa':s,
                          'Sequence distance aa':abs(g-s),'Cα distance Å':v['ca'],'Min heavy-atom distance Å':v['heavy'],
                          'pLDDT germline':[plddt[a] for a in g],'pLDDT second hit':[plddt[b] for b in s],'pLDDT min':conf,
                          'PAE g→s Å':pae[g-1,s-1],'PAE s→g Å':pae[s-1,g-1],'PAE mean Å':(pae[g-1,s-1]+pae[s-1,g-1])/2,
                          'PAE max Å':v['pae'],'Primary confidence':np.where((conf>=70)&(v['pae']<=10),'Yes','No'),
                          'Strict confidence':np.where((conf>=80)&(v['pae']<=5),'Yes','No')})
    pair.to_csv(out/'structural_pairs.csv',index=False)
    reference = pd.read_csv(analysis_dir/'source_data/S24_Structural_pairs.csv').set_index('SGR')
    current = pair.set_index('SGR')
    keep = current.index.intersection(reference.index)
    # The archived validation data use the historical amino-acid anchors.
    # The current supplementary workbook is checked for every revised pair.
    positions = ['Germline aa', 'Second-hit aa']
    matches_archived_position = current.loc[keep, positions].eq(reference.loc[keep, positions]).all(axis=1)
    keep = keep[matches_archived_position]
    for c in reference.select_dtypes('number'):
        np.testing.assert_allclose(current.loc[keep,c], reference.loc[keep,c], rtol=1e-11)
    group_ids = [np.flatnonzero(d.germline.to_numpy()==b) for b in d.germline.drop_duplicates()]
    retained = np.flatnonzero(~d.germline.isin(['R986C','R986H','R1281K']).to_numpy())
    specs = [
        ('Pair-level 3D clustering','Median Cα distance (Å)','ca','median','lower'),
        ('Pair-level 3D clustering','Mean Cα distance (Å)','ca','mean','lower'),
        ('Pair-level 3D clustering','Median heavy-atom distance (Å)','heavy','median','lower'),
        ('Contact thresholds','Pairs Cα <15 Å','ca','count15','upper'),
        ('Contact thresholds','Pairs heavy-atom <15 Å','heavy','count15','upper'),
        ('AlphaFold confidence','Pairs with PAEmax≤5 Å','pae','count5','upper'),
        ('AlphaFold confidence','Pairs with PAEmax≤10 Å','pae','count10','upper'),
        ('Background-collapsed','Median Cα distance (Å)','ca','background','lower'),
        ('Background-collapsed','Median heavy-atom distance (Å)','heavy','background','lower'),
        ('Background-collapsed','Median PAEmax (Å)','pae','background','lower'),
        ('Sensitivity excluding R986C/R986H/R1281K','Median Cα distance (Å)','ca','sensitivity','lower')]
    def reduce(x, kind):
        if kind == 'median': return np.median(x,axis=1)
        if kind == 'mean': return np.mean(x,axis=1)
        if kind == 'count15':return (x<15).sum(axis=1)
        if kind == 'count5': return (x<=5).sum(axis=1)
        if kind == 'count10':return (x<=10).sum(axis=1)
        if kind == 'background': return np.median(np.column_stack([np.median(x[:,ix],axis=1) for ix in group_ids]),axis=1)
        if kind == 'sensitivity': return np.median(x[:,retained],axis=1)
        raise ValueError(kind)
    draws = 100000
    rng = np.random.default_rng(20261004)
    null = np.empty((draws,len(specs)))
    for start in range(0,draws,5000):
        size = min(5000,draws-start)
        ix = np.argsort(rng.random((size,n)),axis=1)
        perm = {k:m[np.arange(n)[None,:],ix] for k,m in matrix.items()}
        for j,(_,_,key,kind,_) in enumerate(specs):
            if kind == 'sensitivity':
                ix2 = np.argsort(rng.random((size,len(retained))),axis=1)
                null[start:start+size,j] = np.median(matrix[key][retained[None,:],retained[ix2]],axis=1)
            else:
                null[start:start+size,j] = reduce(perm[key],kind)
    rows = []
    for j,(analysis,metric,key,kind,tail) in enumerate(specs):
        observed = float(reduce(v[key][None,:],kind)[0])
        sampled = null[:,j]
        extreme = (sampled<=observed).sum() if tail=='lower' else (sampled>=observed).sum()
        p = (1+extreme)/(draws+1)
        eligible = d.iloc[retained] if kind=='sensitivity' else d
        rows.append([analysis,metric,observed,float(sampled.mean() if kind.startswith('count') else np.median(sampled)),p,
                     'Nominal exploratory enrichment.' if p<.05 else 'No statistically significant enrichment of proximity or confidence.',
                     len(eligible),eligible.germline.nunique(),tail])
    pd.DataFrame(rows,columns=['Analysis','Metric','Observed','Null reference','Permutation p','Interpretation','N_configurations','N_backgrounds','Tail']).to_csv(out/'structural_summary.csv',index=False)
    np.savez_compressed(out/'structural_null.npz',draws=null,metrics=[s[1] for s in specs],analyses=[s[0] for s in specs])
