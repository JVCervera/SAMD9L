#!/usr/bin/env python3
"""Run the study's offline statistical workflow and compare with reference results."""
from pathlib import Path
import argparse,subprocess,sys,json,shutil,importlib.metadata
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parent
p=argparse.ArgumentParser()
p.add_argument('--output',type=Path,default=R/'results')
p.add_argument('--figures',action='store_true',help='Also regenerate publication figures and HTML')
p.add_argument('--verify-only',action='store_true',help='Verify saved null distributions without resampling')
a=p.parse_args();O=a.output.resolve();O.mkdir(parents=True,exist_ok=True)
checks=[]
def compare(new,ref,keys,columns):
 n=pd.read_csv(new);r=pd.read_csv(ref)
 if keys:
  n=n.set_index(keys).sort_index();r=r.set_index(keys).sort_index()
  r=r.loc[n.index]
 if len(n)!=len(r):raise ValueError('Row mismatch')
 errors={}
 for c in columns:
  x=n[c].to_numpy(float);y=r[c].to_numpy(float)
  np.testing.assert_allclose(x,y,rtol=1e-8,atol=1e-10,equal_nan=True,err_msg=str(new)+' '+c)
  finite=np.isfinite(x)&np.isfinite(y);errors[c]=float(np.max(abs(x[finite]-y[finite]))) if finite.any() else 0
 checks.append(dict(result=str(new.relative_to(O)),status='PASS',maximum_absolute_errors=errors))
def run(script,label):
 command=[sys.executable,str(R/script),'--output',str(O/label)]
 with (O/(label+'.log')).open('w') as log:
  subprocess.run(command,cwd=O,stdout=log,stderr=subprocess.STDOUT,check=True)
 print('Completed',label,flush=True)
# Validate included chromatin null distributions directly, independently of sampling.
ref=pd.read_csv(R/'reference/chromatin_matched_null_53.csv')
z=np.load(R/'reference/null_distributions.npz',allow_pickle=False)
from analysis.run_match_summary import bh
rows=[]
for row in ref.itertuples():
 v=z[row.test_id];assert v.shape==(10000,)
 center=np.median(v)
 rows.append(dict(test_id=row.test_id,analysis_family=row.analysis_family,null_median=center,null_q025=np.quantile(v,.025),null_q975=np.quantile(v,.975),p_empirical_two_sided=(1+np.count_nonzero(abs(v-center)>=abs(row.observed_stat-center)))/10001))
d=pd.DataFrame(rows);d['q_family_global']=d.groupby('analysis_family').p_empirical_two_sided.transform(bh)
d.to_csv(O/'chromatin_saved_draws.csv',index=False)
compare(O/'chromatin_saved_draws.csv',R/'reference/chromatin_matched_null_53.csv',['test_id'],['null_median','null_q025','null_q975','p_empirical_two_sided','q_family_global'])
if not a.verify_only:
 run('analysis/clinical_statistics/clinical_and_regulatory.py','clinical')
 compare(O/'clinical/clinical_contrasts_full_precision.csv',R/'reference/clinical_contrasts_full_precision.csv',['Analysis','Unit'],['N','Estimate','P'])
 compare(O/'clinical/gene_centered_phase_statistics.csv',R/'reference/gene_centered_phase_statistics.csv',['subset','scope','modality','endpoint'],['n_events_with_scores','n_backgrounds','median','P','q_BH4'])
 compare(O/'clinical/origin_regression_sensitivity.csv',R/'reference/origin_regression_sensitivity.csv',['specification','modality'],['N','Adjusted_beta','CI_low','CI_high','P','q_BH9'])
 run('analysis/run_family.py','family')
 compare(O/'family/family_and_support_sensitivity.csv',R/'reference/family_and_support_sensitivity.csv',['sensitivity'],['N','rho','P_nominal','CI_low','CI_high'])
 run('analysis/run_position_splicing.py','splicing')
 compare(O/'splicing/position_splicing.csv',R/'analysis/source_data/S19_AG_position_splicing.csv',['Context'],['N_primary_alleles','rho','P','q_BH6'])
 run('analysis/run_match_summary.py','rna')
 run('analysis/run_chromatin.py','chromatin')
 compare(O/'chromatin/chromatin_matched_null_53.csv',R/'reference/chromatin_matched_null_53.csv',['test_id'],['n_events','n_backgrounds','observed_stat','null_median','null_q025','null_q975','p_empirical_two_sided','q_family_global'])
 run('analysis/run_structure.py','structure')
 compare(O/'structure/S21_Structural_summary_recomputed.csv',R/'reference/structure_recomputed/S21_Structural_summary_recomputed.csv',['Analysis','Metric'],['Observed','Null reference','Permutation p','N_configurations','N_backgrounds'])
if a.figures:
 shutil.copytree(R/'figure_code',O/'figures',dirs_exist_ok=True)
 for s in ['build_figures.py','build_html.py']:
  with (O/(s+'.log')).open('w') as log:subprocess.run([sys.executable,str(O/'figures'/s)],cwd=O/'figures',stdout=log,stderr=subprocess.STDOUT,check=True)
 for key in json.loads((R/'figure_code/final_legends.json').read_text()):
  for ext in ['pdf','png','svg','html']:assert (O/'figures'/f'{key}.{ext}').stat().st_size>0
 checks.append(dict(result='figures',status='PASS',figures=6))
report=dict(status='PASS',mode='saved-draw verification' if a.verify_only else 'offline reproduction',python=sys.version,versions={n:importlib.metadata.version(n) for n in ['numpy','pandas','scipy','openpyxl','matplotlib','PyMuPDF','Pillow']},external_model_requests=0,checks=checks)
(O/'verification.json').write_text(json.dumps(report,indent=2))
print('PASS. Results:',O)
