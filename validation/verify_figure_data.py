import numpy as np,pandas as pd,re
from pathlib import Path
from lxml import etree as ET
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
import argparse
parser=argparse.ArgumentParser(description="Verify plotted coordinates against the offline reproduction outputs.")
parser.add_argument("output",type=Path,help="Output directory from reproduce.py")
run=parser.parse_args().output
root=Path(__file__).resolve().parents[1]
dir=root/"reference/figures"
reg=run/'regulatory';cl=run/'clinical'
ages=pd.read_csv(cl/'age_rows.csv');roc=pd.read_csv(cl/'roc/roc_input.csv');tests=[]
def svg(n): return ET.parse(str(dir/f'Figure_{n}.svg'))
def axes(n,id):return svg(n).xpath(f'//*[@id="{id}"]')[0]
def xy(us):return np.array([[float(u.get('x')),float(u.get('y'))] for u in us])
def scale(a,k):
 pts=[]
 for g in a.xpath(f'.//*[contains(@id,"{k}tick_")]'):
  ts=g.xpath('.//*[local-name()="text"]');u=g.xpath('.//*[local-name()="use" and @x]')
  if not ts or not u:continue
  try: val=float(''.join(ts[0].itertext()).replace('−','-'))
  except ValueError:continue
  pts.append([float(u[0].get(k)),val])
 if len(pts)<2:raise ValueError([a.get('id'),k,pts])
 return np.polyfit(*np.array(pts).T,1)
def points(a,g):
 return xy(g.xpath('.//*[local-name()="use" and @x]'))
def check(label,actual,want,atol=.001):
 x=np.array(actual,dtype=float);y=np.array(want,dtype=float)
 if x.ndim==1:x=x[:,None];y=y[:,None]
 if x.shape!=y.shape:raise AssertionError((label,x.shape,y.shape))
 dist=cdist(x,y,metric='chebyshev');i,j=linear_sum_assignment(dist);err=dist[i,j].max();assert err<atol,(label,err,atol);print('OK',label,len(x),round(err,6));tests.append((label,len(x),err))
def cmp(n,ax,groups,data,cols):
 a=axes(n,ax)
 for g,df in zip(groups,data):
  x=points(a,a.xpath(f'.//*[@id="{g}"]')[0])
  vs=[np.polyval(scale(axes(n,'axes_1') if n=='S4' else a,dim),x[:,int(dim=='y')]) for dim in cols]
  check(n+' '+g,np.column_stack(vs),df[list(cols.values())])
for n,ax,prefix,endpoint in [('1','B_axes_1','B_','clinical'),('S2','axes_1','','mixed'),('S1','axes_2','','clinical')]:
 cmp(n,ax,[prefix+'PathCollection_'+str(i) for i in [1,2]],[ages[(ages.endpoint==endpoint)&(ages.group==k)] for k in ['Pediatric','Adult']],{'y':'position'} if n=='S1' else {'x':'position','y':'age'})
sahoo=pd.read_csv(cl/'sahoo_baseline.csv')
for i,k in enumerate(['A','B']):
 cmp('S4','axes_'+str(i+1),['PathCollection_'+str(2*i+j) for j in [1,2]],[sahoo[(sahoo.cohort==k)&(sahoo.loss==b)] for b in [False,True]],{'y':'position'})
cmp('S5','A_axes_2',['A_observations_1','A_observations_2'],[roc[roc.chromosome_7_loss==b] for b in [0,1]],{'y':'position_aa'})
d=pd.read_csv(cl/'letter_roc_105/thresholds.csv')
a=axes('S5','B_axes_1');pts=np.vstack([points(a,g) for g in a.xpath('./*[contains(@id,"PathCollection")]')]);pts=np.column_stack([np.polyval(scale(a,k),pts[:,i]) for i,k in enumerate('xy')]);check('S5B residues/Youden',pts,d[['Residue','Youden J']])
a=axes('1','C_axes_2');p=a.xpath('.//*[@id="C_line2d_34"]/*')[0];xy0=np.array(re.findall(r'[ML]\s*([\d.]+)\s+([\d.]+)',p.get('d')),float);true=pd.read_csv(cl/'roc/roc_curve.csv');xy1=np.column_stack([np.polyval(scale(a,k),xy0[:,i]) for i,k in enumerate('xy')]);check('Figure 1C ROC',xy1,true[['false_positive_rate','true_positive_rate']])
raw=root/'archive/analysis'
d=pd.concat([pd.read_csv(raw/'source_data/S12_AG_isolated_scores.csv'),pd.read_csv(raw/'isolated_addendum/S12_dos_alelos_adquiridos.csv')]).drop_duplicates(['dataset','variant_id_project','modality'])
d=d[d.primary_analysis=='Yes'];dn=d[d.modality=='DNASE']
cmp('2','axes_2',['PathCollection_6','PathCollection_7'],[dn[dn.dataset==v] for v in ['GERMLINE','SOMATIC']],{'y':'heme_median_signed_q'})
sp=d[(d.modality=='SPLICE_JUNCTIONS')&(d.dataset=='GERMLINE')]
normalization = pd.read_csv(root/'data/protein_normalizations.csv').set_index('variant_id')
sp = sp.copy()
revised = sp.variant_id_project.map(normalization.residue)
sp.loc[revised.notna(), 'AA_position'] = revised.dropna().astype(int)
cmp('S7','axes_2',['PathCollection_1'],[sp],{'x':'AA_position','y':'delta_p90_q'})
orig=pd.read_csv(reg/'origin_regressions.csv').query('specification=="all_linear"')
cmp('S7','axes_3',['PathCollection_2'],[orig],{'x':'Adjusted_beta'})
gd=pd.read_csv(reg/'gene_backgrounds.csv').query('subset=="all_modelled54"')
for j,groups,mod,endpoint,scopes in [(1,[1,2],'RNA expression','median_rescue_to_WT',['native_all','native_hematopoietic']),(2,[3,4],'RNA expression','fraction_closer_to_WT',['native_all','native_hematopoietic']),(3,[5],'Splice-site usage','median_rescue_to_WT',['native_all']),(4,[6],'Splice-site usage','fraction_closer_to_WT',['native_all'])]:
 cmp('S8','axes_'+str(j),['PathCollection_'+str(i) for i in groups],[gd[(gd.modality==mod)&(gd.scope==s)&(gd.endpoint==endpoint)].dropna(subset=['value']) for s in scopes],{'x':'value'})
chrom=pd.read_csv(reg/'chromatin_tests.csv').query('analysis_family=="PRIMARY"')
cmp('S8','axes_7',['PathCollection_'+str(i) for i in (15,16,17)],[chrom[chrom.output_type==m] for m in ('DNASE','CHIP_HISTONE','CHIP_TF')],{'y':'q_family_global'})
st=pd.read_csv(reg/'structural_pairs.csv')
pts=svg('S6').xpath('//*[starts-with(@id,"structural_SGR")]')
check('S6C structural pairs',np.array([[(float(p.get('cx'))-336.756)/.115476692,(443.72-float(p.get('cy')))/1.213] for p in pts]),st[['Sequence distance aa','Cα distance Å']])
for name,x0,x1 in [('first',38.267717,237.25984),('hematologic',301.03938,489.82679)]:
 for group in ['Pediatric','Adult']:
  a=axes('S3',name+'-'+group.lower());paths=a.xpath('.//*[local-name()="path"]')[::2];arr=np.array([p.get('transform')[7:-1].split(',')[-2:] for p in paths],float)
  check('S3 '+name+group,np.c_[(arr[:,0]-x0)*1600/(x1-x0),(161.5039-arr[:,1])*90/109.7008],ages[(ages.endpoint==name)&(ages.group==group)][['position','age']])
rna=pd.read_csv(reg/'rna_backgrounds.csv');rna['value']=1000*rna.advantage
for i,subset in enumerate(rna.subset.unique()):
 cmp('S8','axes_6',['PathCollection_'+str(7+i*2+j) for j in range(2)],[rna[(rna.subset==subset)&(rna.scope==s)] for s in ('native_all','native_hematopoietic')],{'x':'value'})
def bar(n,g,i):
 p=axes(n,'patch_'+str(i))[0];co=np.array(re.findall(r'[ML]\s*([-\d.]+)\s+([-\d.]+)',p.get('d')),float);return np.polyval(scale(axes(n,'axes_'+str(g[0])),g[1]),co[2,int(g[1]=='y')])
t=pd.read_csv(reg/'tissue_summary.csv').set_index(['dataset','modality']);mods=['RNA_EXPRESSION','CAGE','SPLICE_JUNCTIONS','SPLICE_SITE_USAGE','POLYADENYLATION','ATAC','DNASE','CHIP_HISTONE','CHIP_TF']
for j,origin in enumerate(('GERMLINE','SOMATIC')):
 for i,mod in enumerate(mods):
  row=t.loc[origin,mod];assert abs(bar('S7',(1,'x'),3+9*j+i)-100*row.n_heme_higher/row.n_alleles)<.001
for j,bg in enumerate(['I876T','I891T','R986C','R986H','R1281K','V1512M']):
 row=rna.query('subset=="confirmed12" and scope=="native_all" and background==@bg');assert abs(bar('S7',(4,'y'),32+j)-row.value.iloc[0])<.001
from collections import Counter
import fitz
labels=[''.join(e.itertext()).strip() for e in axes('1','A_axes_3').xpath('.//*[local-name()="text"]')]
coords=pd.read_csv(raw/'source_data/Figure_1_display_coordinates.csv');wanted=[x.replace('(n=','(') for x in coords.Label]+['E685del','K1265N']
print('missing Fig1A allele labels',Counter(wanted)-Counter(labels))
cc=pd.read_csv(cl/'clonal.csv').query('any_recorded');assert cc[['chr7','upd','intra']].sum().tolist()==[60,26,38];print('Figure 2A counts',cc.groupby(['chr7','upd','intra']).size().tolist())

assert not (Counter(wanted)-Counter(labels))
assert sorted(cc.groupby(['chr7','upd','intra']).size())==[6,6,6,8,8,18,30]
a=axes('S6','panel_B')
co=a.xpath('.//*[local-name()="circle"]')[:4];rs=a.xpath('.//*[local-name()="rect"]')[:4]
s=pd.read_csv(reg/'structural_summary.csv');q=s[(s.Analysis.isin(['Pair-level 3D clustering','Background-collapsed']))&s.Metric.isin(['Median Cα distance (Å)','Median heavy-atom distance (Å)'])]
np.testing.assert_allclose([(float(c.get('cx'))-117.354)/2.04094285 for c in co],q.Observed,atol=.001)
np.testing.assert_allclose([(float(c.get('x'))+2.7-117.354)/2.04094285 for c in rs],q['Null reference'],atol=.001)
ev=pd.read_csv(reg/'events.csv').set_index('event')
for p in axes('S6','panel_A').xpath('.//*[starts-with(@id,"SGR")]'):
 name=p.get('id');row=ev.loc[name[:6]]
 x,y=[float(p.get(k)) for k in ('cx','cy')] if p.get('cx') else np.array([v.split(',') for v in p.get('points').split()],float).mean(axis=0)
 assert abs((199.33-y)/.09885625-row.somatic_pos)<.001
 if pd.notna(row.germline_pos):assert abs((x-38.268)/.12437-row.germline_pos)<.001
g=pd.read_csv(cl/'gnomad.csv')
def s3markers(prefix):
 pp=svg('S3').xpath(f'//*[local-name()="path" and starts-with(@d,"{prefix}") and @fill]')
 return np.array([e.get('transform')[7:-1].split(',')[-2:] for e in pp],float)
e=s3markers('M0-1.5C')[::2]
for side,subset in [(0,g[g.af.notna()]),(1,g[g.af.isna()])]:
 vals=(e[(e[:,0]>130)==bool(side),1]);base=609.4488-67.03937
 # Tick spacing is 250 amino acids.
 check('S3E '+str(side),(base-vals)*250/15.62689,subset.position,atol=.01)
f=s3markers('M0-2.061')[::2];obs=g.dropna(subset=['af']).sort_values('position')
pos=(f[:,0]-301.03938)*1600/188.78741
check('S3F detected positions',pos,obs.position)


check('S3F frequency and position',np.c_[pos, -6+(609.4488-80.89917-f[:,1])/40.01120],np.c_[obs.position,np.log10(obs.af)])
from scipy.stats import spearmanr
for endpoint in ['clinical','first','hematologic','mixed']:
 table=ages[ages.endpoint==endpoint];rho,p=spearmanr(table.age,table.position)
 assert f'{rho:.3f}'.replace('-','−') in fitz.open(dir/'Figure_S3.pdf')[0].get_text()
t=cc.groupby('clinical_group')[['chr7','upd','intra']].sum()
np.testing.assert_array_equal(t.loc[['Pediatric','Adult','Unknown']].values,[[50,22,32],[5,2,2],[5,2,4]])
a=ev[ev.somatic_cdna=='c.1765C>T']
assert dict(zip(a.germline_label,a.phase))=={'p.Lys770Glu':'reported_cis','p.Thr936Ile':'unknown','p.Arg986Cys':'inferred_cis','p.Arg1281Lys':'confirmed_cis','p.Thr1474Pro':'reported_cis'}
