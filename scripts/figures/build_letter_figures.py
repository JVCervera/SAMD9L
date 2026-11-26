"""Rebuild the two Letter and eight supplementary figures from the adjudicated tables.

The frozen workbook supplies the clinical and allele records. Statistical annotations
are read from the clinical and regulatory reruns. All figures are re-drawn.
"""
import argparse
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
from matplotlib import font_manager, lines, patches
import numpy as np
import pandas as pd
from openpyxl import load_workbook
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score, roc_curve

COLORS = {"Missense": "#285c64", "Truncating": "#a84826", "In-frame indel": "#745485"}
TEAL, RUST, PLUM = COLORS.values()
DARK, GRAY, LIGHT, GRID = "#20292c", "#61686a", "#a7afac", "#dfe3df"
WIDTH = 175 / 25.4
DOMAINS = [("SAM",12,76,"#b3d2ca"),("AlbA",176,380,"#c1d2a5"),("Sir2-like",395,624,"#d7c0d5"),("P-loop NTPase",629,1040,"#e3bcb6"),("TPR",1056,1178,"#bdd2dc"),("Helical",1179,1520,"#ead6b2"),("OB",1521,1584,"#dbd9a9")]
AA = dict(zip("Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val Ter".split(),"ARNDCQEGHILKMFPSTWYV*"))
MODS = {"ATAC":"ATAC", "CAGE":"CAGE", "DNASE":"DNase", "CHIP_HISTONE":"Histone ChIP", "HISTONE_CHIP":"Histone ChIP", "POLYADENYLATION":"Polyadenylation","RNA_EXPRESSION":"RNA expression","SPLICE_JUNCTIONS":"Splice junctions","SPLICE_SITE_USAGE":"Splice-site usage","SPLICE_SITES":"Splice sites","CHIP_TF":"TF ChIP","TF_CHIP":"TF ChIP"}
def mod(s): return MODS.get(s,s.replace("_"," ").capitalize())
def protein(s):
    s=str(s).replace("p.","").replace("(","").replace(")","")
    for k,v in AA.items():s=s.replace(k,v)
    return s.replace("X","*")
def numpos(s):
    m=re.search(r'\d+',str(s))
    if not m: raise ValueError(f"Could not map residue for {s}")
    return int(m.group())
def posnum(t):
    try:return float(t)
    except (TypeError,ValueError):return np.nan
def pval(p,prec=3):
    value=f"{float(p):.{prec}g}"
    if 'e' not in value:return value
    mantissa,exponent=value.split('e')
    return mantissa + ' × 10$^{' + str(int(exponent)) + '}$'
def pt(ax,title,letter=None):
    ax.set_title(title,loc='left',fontweight='bold',fontsize=9,pad=8)
    if letter:
        ax.annotate(letter,(0,1),xycoords='axes fraction',xytext=(-18,7),textcoords='offset points',fontsize=11,fontweight='bold',annotation_clip=False)
def fig(height): return plt.figure(figsize=(WIDTH,height/25.4),dpi=120,layout='constrained')
def axstyle(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(direction='out',length=3,width=.65,pad=3)
    ax.grid(axis='y',lw=.45,color=GRID)
    ax.set_axisbelow(True)
def section(w,si,start,second):
    a=w[si]; idx=next(i for i,row in enumerate(a) if len(row)>1 and str(row[0]).strip()==start and str(row[1]).strip()==second)
    cols=[str(x) if x is not None else f"_empty{k}" for k,x in enumerate(a[idx])]
    rows=[]
    for row in a[idx+1:]:
        if not row or row[0] is None or re.match(r'^[A-Z]\d?\. ',str(row[0])): break
        row = list(row) + [None]*max(0,len(cols)-len(row))
        rows.append(row[:len(cols)])
    return pd.DataFrame(rows,columns=cols).drop(columns=[c for c in cols if c.startswith('_empty')])
def box(ax,groups,names,colors,label,integer=False):
    for i,(vals,col) in enumerate(zip(groups,colors)):
        v=np.sort(np.asarray(vals,dtype=float))
        if not len(v):continue
        off=np.zeros(len(v))
        # Spread ties symmetrically; y coordinates always remain the observed values.
        if integer:
            for u in np.unique(v):
                sel=v==u; k=sum(sel); off[sel]=(np.arange(k)-(k-1)/2)*min(.07,.48/max(k-1,1))
        else:
            off=np.random.default_rng(418+i).uniform(-.19,.19,len(v))
        ax.boxplot(v,positions=[i],widths=.48,showfliers=False,patch_artist=True,manage_ticks=False,boxprops={'facecolor':col+'10','edgecolor':col,'linewidth':.8},medianprops={'color':DARK,'linewidth':1.2},whiskerprops={'color':col,'linewidth':.8},capprops={'color':col,'linewidth':.8},zorder=1)
        ax.scatter(i+off,v,s=11,c=col,edgecolor='white',lw=.3,zorder=2)
    ax.set_xticks(range(len(names)),names); ax.set_xlim(-.52,len(names)-.48);ax.set_ylabel(label);axstyle(ax)
def export(f,path,name):
    for ext in ['pdf','svg','png']+(['tif'] if name in ['Figure_1','Figure_2'] else []):
        extra={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tif' else {}
        target=path/f'{name}.{ext}'
        if ext == 'pdf':
            from datetime import datetime, timezone
            extra['metadata'] = {'CreationDate': datetime(2000, 1, 1, tzinfo=timezone.utc)}
        elif ext == 'svg':
            extra['metadata'] = {'Date': None}
        f.savefig(target,dpi=1000 if ext=='tif' else 400,facecolor='white',**extra)
        if ext == 'tif':
            with Image.open(target) as source:
                rgb = source.convert('RGB')
            rgb.save(target, compression='tiff_lzw', dpi=(1000, 1000))
        if ext=='svg' and plt.rcParams['font.family'][0]=='Liberation Sans':
            svg=target.read_text().replace("'Liberation Sans'","'Arial', 'Liberation Sans', sans-serif")
            target.write_text(svg)
    plt.close(f)

def load_data(args):
    wb=load_workbook(args.workbook,data_only=True,read_only=True)
    w=[list(s.values) for s in wb];wb.close()
    a=section(w,0,'Patient ID','Cohort')
    if 'Mixed reported age (years)' in a:
        a = a.rename(columns={'Mixed reported age (years)': 'Recurated mixed age (years)'})
    links=section(w,1,'Patient Master ID','Integrated variant ID'); ev=section(w,1,'Event ID','Patient ID')
    inv=section(w,1,'Origin','Variant ID'); germ=inv[(inv.Origin=='Germline')&(inv['Primary analysis']=='Yes')].copy()
    som=inv[inv.Origin.isin(['Acquired','Somatic'])].copy()
    gc=links[links['Integrated variant ID'].isin(germ['Variant ID'])].groupby('Integrated variant ID')['Patient Master ID'].nunique()
    a['position']=pd.to_numeric(a['Amino-acid position'],errors='coerce')
    primary=links[links['Integrated variant ID'].isin(germ['Variant ID'])].groupby('Patient Master ID')['Integrated variant ID'].agg(lambda x:x.iloc[0] if len(set(x))==1 else None)
    # Single-allele individuals can have position annotated only in S2.
    a['allele']=a['Patient ID'].map(primary)
    a['position']=a.position.fillna(a.allele.map(germ.set_index('Variant ID')['AA position']))
    labels=pd.read_csv(args.source/'Figure_1_display_coordinates.csv').set_index('Variant_ID')
    def normalize_cdna(x):
        return str(x).split(':')[-1].strip()
    def somkey(r):
        c=normalize_cdna(r['Acquired cDNA'])
        match=som[som['AlphaGenome cDNA'].map(normalize_cdna)==c] if c not in ['None','nan',''] else som.iloc[:0]
        if len(match)!=1:
            names=[protein(r['Normalized protein']),protein(r['Protein as reported'])]
            match=som[som['AlphaGenome protein'].map(protein).isin(names)]
        if len(match)==1:return match.iloc[0]['Variant ID']
        return None
    ev['vid']=ev.apply(somkey,axis=1)
    missing=ev[ev.vid.isna()][['Event ID','Acquired cDNA','Normalized protein']].to_dict('records')
    rows=[]
    for cat,frame,counts in [('Germline',germ,gc),('Acquired',som,ev.vid.value_counts())]:
        for _,r in frame.iterrows():
            vid=r['Variant ID']; count=counts.get(vid,0)
            if count<=0:continue
            cl=r['Functional class']
            con='Missense' if 'missense' in str(cl).lower() else 'In-frame indel' if bool(re.search('in[ -]?frame',str(cl),re.I)) else 'Truncating'
            # Reuse verified original descriptions for equivalent alleles with different reported protein labels.
            label=protein(r['AlphaGenome protein'])
            if 'Label' in labels.columns and vid in labels.index:
                label=re.sub(r'\s*\(n=\d+\)\s*$','',labels.loc[vid,'Label']).strip()
            protein_only = cat == 'Acquired' and not normalize_cdna(r['AlphaGenome cDNA']).startswith('c.')
            rows.append(dict(origin=cat,vid=vid,count=int(count),pos=float(r['AA position']),label=label,con=con,protein_only=protein_only))
    land=pd.DataFrame(rows)
    clinical=a.loc[a['Age-position eligible'].eq('Yes')].copy()
    mx=a.loc[a['Mixed age-position eligible'].eq('Yes')].copy()
    for data,key in [(clinical,'Clinical context age (years)'),(mx,'Recurated mixed age (years)')]:
        data['age']=pd.to_numeric(data[key],errors='raise')
        if data.age.isna().any() or data.position.isna().any():raise ValueError(f'Missing eligible numeric age/position for {key}')
    return w,a,ev,land,clinical,mx,missing

def spread(true,minimum_gap,lo=18,hi=1566):
    # Nearest nonoverlapping x slots while keeping the original genomic order.
    from scipy.optimize import isotonic_regression
    n=len(true); dx=min(minimum_gap,(hi-lo)/max(1,n-1)*.978)
    first=np.arange(n)*dx
    corrected=isotonic_regression(np.asarray(true)-first).x
    return np.clip(corrected,lo,hi-first[-1])+first

def landscape(ax,data,labels,size_title,top_title):
    lower,upper=labels
    ax.set_xlim(-30,1615);ax.set_ylim(-40,43);ax.set_yticks([])
    ax.xaxis.tick_top();ax.set_xticks([1,200,400,600,800,1000,1200,1400,1584]);ax.tick_params(axis='x',length=2,pad=1,labelsize=6.6)
    ax.spines[['bottom','left','right']].set_visible(False);ax.spines['top'].set_color(LIGHT)
    ax.add_patch(patches.Rectangle((1,-1.35),1583,2.7,facecolor='#edf0ec',edgecolor=GRAY,lw=.5,zorder=3))
    for name,lo,hi,color in DOMAINS:
        ax.add_patch(patches.Rectangle((lo,-1.35),hi-lo,2.7,facecolor=color,edgecolor='white',lw=.5,zorder=4))
        ax.text((lo+hi)/2,0,name,ha='center',va='center',fontsize=5.9,zorder=5)
    rows=[]
    for sign,group in [(1,upper),(-1,lower)]:
        sub=data[data.origin.eq(group)].sort_values(['pos','vid'])
        # Offset the labels only as far as their print size requires. The bases
        # always remain at the measured protein residues, including tied alleles.
        label_font = 6.0
        marker_scale = 6.2
        points_per_x = ax.bbox.width * 72 / ax.figure.dpi / np.diff(ax.get_xlim())[0]
        points_per_y = ax.bbox.height * 72 / ax.figure.dpi / np.diff(ax.get_ylim())[0]
        label_gap = (label_font + 1.4) / points_per_x
        xpos = spread(sub.pos, label_gap / 2 + 0.35)
        renderer = ax.figure.canvas.get_renderer()
        font = font_manager.FontProperties(size=label_font)
        previous_inner = -np.inf
        arranged = []
        for x, (_, r) in zip(xpos, sub.iterrows()):
            n = r['count']
            protein_only = bool(r.get('protein_only', False))
            label = r.label + (f' ({n})' if n > 1 else '') + ('§' if protein_only else '')
            extent = renderer.get_text_width_height_descent(label, font, ismath=False)[0]
            length_points = extent * 72 / ax.figure.dpi
            inner = (x - previous_inner >= label_gap)
            if inner:
                previous_inner = x
            arranged.append((x, r, label, length_points, inner))
        for x, r, label, length_points, inner in arranged:
            c = COLORS[r.con]
            n = r['count']
            size = marker_scale * n
            adjacent = [(length + np.sqrt(marker_scale * other['count']) / 2)
                        for ox, other, _, length, inner_neighbor in arranged
                        if inner_neighbor and abs(ox - x) < label_gap]
            extra = 0 if inner else max(adjacent, default=0) + 4
            base = 1.35 + (7.0 + extra) / points_per_y
            ax.plot([r.pos, r.pos, x, x], np.array([1.35, 2.45, 1.35 + 5.2 / points_per_y, base]) * sign,
                    color=c, lw=.55, zorder=1)
            protein_only = bool(r.get('protein_only', False))
            ax.scatter(x, sign * base, s=size, facecolor='white' if protein_only else c,
                       edgecolor=c if protein_only else 'white', lw=.85 if protein_only else .3, zorder=3)
            ax.annotate(label, (x, sign * base), xytext=(0, sign * (np.sqrt(size) / 2 + 1.9)),
                        textcoords='offset points', va='bottom' if sign > 0 else 'top', ha='center',
                        rotation=90 * sign, fontsize=label_font, color=DARK, zorder=4,
                        bbox=dict(facecolor='white', edgecolor='none', boxstyle='square,pad=0.08'))
            rows.append(dict(origin=group, allele=r.vid, position=r.pos, display_x=x,
                             count=n, area_pt2=size, label=label))
        ax.text(.002,.984 if sign>0 else .015,f'{group}  ·  {len(sub)} alleles',transform=ax.transAxes,ha='left',va='top' if sign>0 else 'bottom',fontweight='bold',fontsize=7.5)
    handles=[lines.Line2D([],[],ls='',marker='o',ms=4,mfc=c,mec=c,label=k) for k,c in COLORS.items() if k in set(data.con)]
    l1=ax.legend(handles=handles,loc='lower left',bbox_to_anchor=(0,-.057),ncol=len(handles),frameon=False,fontsize=6.7,handlelength=.65,handletextpad=.4,columnspacing=1)
    ax.add_artist(l1)
    if 'protein_only' in data and data.protein_only.any():
        ax.text(.0,-.080,'§  Protein-level report only',transform=ax.transAxes,fontsize=6.6,ha='left')
    handles=[lines.Line2D([],[],ls='',marker='o',ms=np.sqrt(6.2*n),mfc=GRAY,mec='white',mew=.27,label=str(n)) for n in [1,2,5,10]]
    ax.legend(handles=handles,title=size_title,loc='lower right',bbox_to_anchor=(1,-.076),ncol=4,frameon=False,fontsize=6.7,title_fontsize=6.7,handlelength=1.0,columnspacing=.8,handletextpad=.3)
    pt(ax,top_title)
    return rows

def s1_flags(a):
    return np.column_stack([a['Chromosome 7 loss/deletion'].eq('Reported'),pd.to_numeric(a['Intragenic event observations'],errors='raise').gt(0),a.UPD7q.eq('Reported')])

def roc_inputs(w,roc):
    rr=section(w,2,'Patient ID (S1)','Primary protein label')
    rr=rr[rr['ROC outcome (loss = 1)'].isin([0,1])].copy()
    x=rr['Position (aa)'].astype(float);y=rr['ROC outcome (loss = 1)'].astype(int)
    auc=roc_auc_score(y,x)
    assert len(rr)==roc['n'] and y.sum()==roc['loss'] and np.isclose(auc,roc['apparent_position_auc'])
    return x,y,int(roc['cutoff']['position_aa']),auc,*roc['apparent_position_auc_95_percentile_ci']

def scatter_age(ax,data,label):
    for selected,color,marker,name in [(data.age.le(18),TEAL,'o','Pediatric'),(data.age.gt(18),RUST,'^','Adult')]:
        ax.scatter(data.position[selected],data.age[selected],s=15,color=color,marker=marker,edgecolor='white',lw=.35,zorder=3,label=f'{name} ({selected.sum()})',clip_on=False)
    ax.axhline(18,color=GRAY,lw=.7,ls=(0,(3,3)))
    ax.set(xlim=(0,1600),ylim=(0,max(65,data.age.max()+2)),xticks=[0,400,800,1200,1600],xlabel='Germline residue position',ylabel=label)
    axstyle(ax);ax.legend(loc='upper right',fontsize=6.6,handlelength=1.1,borderaxespad=.25)
    r,p=spearmanr(data.age,data.position);return r,p

def figure1(w,land,clinical,roc,out):
    f=plt.figure(figsize=(WIDTH,172/25.4),dpi=120)
    ax=f.add_axes([.055,.449,.927,.465]); rows=landscape(ax,land,('Acquired','Germline'),'Carriers / observations per allele','')
    f.text(.014,.961,'A',fontsize=11,fontweight='bold');f.text(.067,.961,'Germline and acquired SAMD9L variation',fontsize=10,fontweight='bold')
    f.text(.982,.960,'NM_152703.5  ·  1,584 residues',ha='right',fontsize=7)
    a1=f.add_axes([.093,.101,.345,.206]);a2=f.add_axes([.668,.085,.275,.212])
    r,p=scatter_age(a1,clinical,'Clinical-context age (years)')
    a1.legend(loc='upper center',bbox_to_anchor=(.5,-.24),ncol=2,fontsize=6.6,handlelength=.8,columnspacing=1.2)
    f.text(.014,.355,'B',fontsize=11,fontweight='bold');f.text(.067,.355,'Germline position and clinical-context age',fontsize=9,fontweight='bold')
    f.text(.067,.331,f'Spearman ρ = {r:.3f}   P = {pval(p)}',fontsize=7.2)
    x,y,_,auc,lo,hi=roc_inputs(w,roc);fp,tp,_=roc_curve(y,x,drop_intermediate=False)
    a2.set_aspect('equal',adjustable='box');a2.plot([0,1],[0,1],ls=(0,(3,3)),color=GRAY,lw=.7)
    a2.plot(fp,tp,lw=1.5,color=TEAL)
    a2.set(xlim=(0,1),ylim=(0,1),xticks=np.linspace(0,1,6),yticks=np.linspace(0,1,6),xlabel='1 − specificity',ylabel='Sensitivity');axstyle(a2)
    f.text(.524,.355,'C',fontsize=11,fontweight='bold');f.text(.577,.355,'ROC for −7/del(7q)',fontsize=9,fontweight='bold')
    f.text(.577,.331,f'AUC {auc:.3f} (95% CI {lo:.3f}–{hi:.3f})',fontsize=7.2)
    export(f,out,'Figure_1');pd.DataFrame(rows).to_csv(out/'Figure_1_positions.csv',index=False)

def phase_groups(identifier):
    s=str(identifier).strip().lower()
    if 'both' in s:return 'Cis and alternative phase'
    if 'trans' in s:return 'Confirmed trans'
    if 'dual' in s:return 'Two germline alleles'
    if 'confirm' in s:return 'Confirmed cis'
    if 'report' in s:return 'Reported cis'
    if 'infer' in s:return 'Inferred cis'
    if s in ['unknown','unresolved','undetermined','not established']:return 'Unresolved'
    raise ValueError(f'Unmapped phase: {s}')
PHASE_STYLE={'Confirmed cis':('#276e97','o','-'),'Cis and alternative phase':('#524575','D','-'),'Reported cis':('#20776c','s','--'),'Inferred cis':(RUST,'^','-.'),'Confirmed trans':('#913d73','X',':'),'Two germline alleles':('#34393a','D','-'),'Unresolved':('#737b7e','o',':')}

def figure2(a,phase,scores,reg,out):
    flags=s1_flags(a);pat=np.array([(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1)])
    counts=[np.all(flags==l,axis=1).sum() for l in pat];marg=flags.sum(axis=0)
    f=plt.figure(figsize=(WIDTH,149/25.4),dpi=120)
    f.text(.014,.95,'A',fontsize=11,weight='bold');f.text(.067,.95,'Combinations of documented clonal findings',fontsize=9.5,weight='bold')
    ax=f.add_axes([.27,.72,.695,.18]); axstyle(ax)
    ax.bar(range(7),counts,color=[TEAL]*3+['#57837c']*3+[RUST],width=.51)
    for ix,n in enumerate(counts):ax.text(ix,n+.5,str(n),ha='center',va='bottom',fontsize=7)
    ax.set(xticks=[],xlim=(-.5,6.5),ylim=(0,max(counts)*1.17),ylabel='Individuals')
    mx=f.add_axes([.27,.56,.695,.125]);mx.set(xlim=(-.5,6.5),ylim=(2.55,-.55));mx.axis('off')
    for j,name in enumerate(['−7/del(7q)','Intragenic variant','UPD7q']):
        mx.text(-.6,j,f'{name}  ({marg[j]})',ha='right',va='center',fontsize=7.5)
        mx.plot([-.4,6.45],[j,j],color=GRID,lw=.55)
    for i,r in enumerate(pat):
        ids=np.flatnonzero(r);mx.plot([i,i],[ids.min(),ids.max()],color=TEAL if i<6 else RUST,lw=1.4)
        for j,v in enumerate(r):mx.scatter(i,j,s=23,facecolor=(TEAL if i<6 else RUST) if v else 'white',edgecolor=(TEAL if i<6 else RUST) if v else LIGHT,lw=.6,zorder=2)
    for xx,txt in [(1,'One finding'),(4,'Two findings'),(6,'Three')]:mx.text(xx,2.64,txt,va='top',ha='center',fontsize=7)
    ax=f.add_axes([.055,.085,.475,.336]);ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
    f.text(.014,.465,'B',fontsize=11,weight='bold');f.text(.067,.465,'p.Arg589Ter germline backgrounds',fontsize=9.5,weight='bold')
    stop=phase[(phase.somatic_pos==589)&(phase.somatic_label.map(protein)=='R589*')].sort_values('germline_pos')
    assert len(stop)==len(set(stop.patient))
    for x,l in [(0,'Acquired'),(.338,'Germline'),(.70,'Phase')]:ax.text(x,.97,l,weight='bold',fontsize=7)
    ax.text(0,.49,'p.Arg589Ter',va='center',fontsize=7.3,weight='bold')
    for (_,row),y in zip(stop.iterrows(),np.linspace(.81,.12,len(stop))):
        lab=phase_groups(row.phase);c,mark,ls=PHASE_STYLE[lab]
        ax.plot([.247,.317],[.49,y],c=c,ls=ls,lw=.9);ax.scatter(.317,y,c=c,s=9,marker=mark)
        ax.text(.338,y,row.germline_label,va='center',fontsize=7.15);ax.text(.70,y,lab,va='center',c=c,fontsize=6.9)
    filt=scores[scores.modality=='DNASE']; groups=[filt.loc[filt.dataset==k,'heme_median_signed_q'].dropna() for k in ['GERMLINE','SOMATIC']]
    model=pd.read_csv(reg/'origin_regressions.csv')
    m=model.loc[(model.specification=='all_linear') & (model.modality=='DNASE')].squeeze()
    f.text(.55,.465,'C',fontsize=11,weight='bold');f.text(.600,.465,'Hematopoietic DNase sensitivity',fontsize=9.2,weight='bold')
    f.text(.600,.425,f'Adjusted β = {m.Adjusted_beta:.3f}  (95% CI {m.CI_low:.3f} to {m.CI_high:.3f})',fontsize=7)
    f.text(.600,.396,f'BH-adjusted P = {pval(m.q_BH9)}',fontsize=7)
    ax=f.add_axes([.675,.115,.29,.245])
    box(ax,groups,[f'Germline\n(n = {len(groups[0])})',f'Acquired\n(n = {len(groups[1])})'],[TEAL,RUST],'Median signed quantile')
    ax.axhline(0,lw=.65,color=GRAY); ax.set_ylim(-1.08,1.08);ax.set_yticks([-1,-.5,0,.5,1])
    export(f,out,'Figure_2')

def figure_s1(a,land,cl,summary,clonal,out):
    f=fig(141);g=f.add_gridspec(2,2,height_ratios=[.63,1],wspace=.14,hspace=.18)
    ax=f.add_subplot(g[0,:]);ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off');pt(ax,'Source records and retained cohort','A')
    # The historical screening flow is fixed and is independent of clinical re-annotation.
    nodes=[(.016,'195','Source records'),(.272,'181','Eligible source records'),(.526,str(len(a)),'Retained individuals'),(.783,str(sum(land.origin=='Germline')),'Primary germline alleles')]
    for i,(x,n,txt) in enumerate(nodes):
        ax.add_patch(patches.FancyBboxPatch((x,.41),.192,.32,boxstyle='round,pad=.010,rounding_size=.02',facecolor='#f2f5f2',ec=LIGHT,lw=.65))
        ax.text(x+.096,.631,n,ha='center',va='center',fontsize=12,fontweight='bold',c=TEAL);ax.text(x+.096,.482,txt,ha='center',va='center',fontsize=6.6)
        if i<3:ax.annotate('',xy=(x+.247,.57),xytext=(x+.207,.57),arrowprops={'arrowstyle':'->','color':GRAY,'lw':.8})
    ax.text(.263,.15,'12 unresolved records\n2 unconfirmed relatives',ha='center',fontsize=6.7)
    ax.text(.514,.15,'2 duplicate records\n1 primary SAMD9 case',ha='center',fontsize=6.7)
    ax=f.add_subplot(g[1,0]);data=[cl.loc[cl.age<=18,'position'],cl.loc[cl.age>18,'position']]
    box(ax,data,[f'Pediatric ({len(data[0])})',f'Adult ({len(data[1])})'],[TEAL,RUST],'Germline residue position',True)
    pt(ax,'Position by clinical-context age','B')
    ax.text(.98,.97,f'P = {pval(summary["age"]["clinical"]["mw"])}',transform=ax.transAxes,ha='right',va='top')
    ax.set_ylim(0,1700)
    ax=f.add_subplot(g[1,1]);cc=clonal.groupby('clinical_group')[['chr7','intra','upd']].sum()
    for i,(key,lab,c) in enumerate([('chr7','−7/del(7q)',TEAL),('intra','Intragenic',RUST),('upd','UPD7q',PLUM)]):
        order=['Pediatric','Adult','Unknown'];v=cc.loc[order,key].astype(int);xs=np.arange(3)+(i-1)*.23
        ax.bar(xs,v,width=.21,color=c,label=lab)
        for x,n in zip(xs,v):ax.text(x,n+.5,str(n),fontsize=6.6,va='bottom',ha='center')
    ax.set(xticks=[0,1,2],xticklabels=['Pediatric','Adult','Age NA'],ylabel='Individuals',ylim=(0,cc[['chr7','intra','upd']].to_numpy().max()*1.25))
    pt(ax,'Clonal findings by age','C');axstyle(ax);ax.legend(loc='upper right',fontsize=6.5)
    export(f,out,'Figure_S1')

def figure_s2(reported,out):
    f = plt.figure(figsize=(5.27, 4.88), dpi=120)
    ax = f.add_axes([.13, .15, .83, .73])
    rho, p = scatter_age(ax, reported, 'Broader reported age (years)')
    f.text(.13, .95, 'Broader reported age and germline position', fontsize=9.5, weight='bold')
    f.text(.13, .90, f'Spearman ρ = {rho:.3f}   P = {pval(p)}', fontsize=7.2)
    export(f, out, 'Figure_S2')

def figure_s3(w,clinical_rows,summary,out):
    age=clinical_rows;f=fig(179);g=f.add_gridspec(3,2,wspace=.1,hspace=.15)
    for j,(name,letter,col,caption) in enumerate([('First event','A','First event age (years)','First clinical event'),('Hematologic event','B','Hematologic event age (years)','Hematologic event')]):
        ax=f.add_subplot(g[0,j]);part=age[age.endpoint.eq({'First event':'first','Hematologic event':'hematologic'}[name])];assert len(part)>0,(name,age.endpoint.unique())
        r,p=scatter_age(ax,part,col);ax.set_ylim(0, max(70,part.age.max()*1.23));pt(ax,caption,letter);ax.text(.02,.96,f'ρ = {r:.3f}  ·  P = {pval(p)}',va='top',transform=ax.transAxes,fontsize=6.5);ax.legend_.remove()
    ax=f.add_subplot(g[1,0]);sa=summary['age']
    label={'clinical':'Clinical context','first':'First event','hematologic':'Hematologic event','mixed':'Broader reported age'}
    def effects(ax,names,rho):
        ax.scatter(rho,range(len(names)),c=TEAL,s=22)
        for i,r in enumerate(rho):ax.annotate(f'{r:.3f}',(r,i),xytext=(5,0),textcoords='offset points',va='center',fontsize=6.5)
        ax.set(yticks=np.arange(len(names)),yticklabels=names,ylim=(len(names)-.5,-.5),xlim=(min(-.63,min(rho)-.09),.07),xlabel='Spearman ρ');ax.axvline(0,lw=.65,c=GRAY);axstyle(ax)
    effects(ax,[f'{lab} ({sa[key]["n"]})' for key,lab in label.items()],[sa[key]['rho'] for key in label]);pt(ax,'Sensitivity to age definition','C')
    ax=f.add_subplot(g[1,1]);sources=summary['source'];source_order=['Other','St Jude']
    effects(ax,[f'{("Other sources" if key=="Other" else "St Jude")} ({sources[key]["n"]})' for key in source_order],[sources[key]['rho'] for key in source_order]);pt(ax,'Source attribution','D')
    inv=section(w,1,'Origin','Variant ID');ger=inv[(inv.Origin=='Germline') & (inv['Primary analysis']=='Yes')]
    seen=ger['gnomAD status'].astype(str).str.startswith('Observed')
    unseen=ger['gnomAD status'].astype(str).str.startswith('Not observed')
    ax=f.add_subplot(g[2,0]); vs=[ger.loc[seen,'AA position'],ger.loc[unseen,'AA position']]
    box(ax,vs,[f'Observed ({len(vs[0])})',f'Unobserved ({len(vs[1])})'],[TEAL,RUST],'Germline residue position',True)
    pt(ax,'gnomAD v4.1.1 detection','E');ax.set_ylim(0,1860);ax.text(.99,.96,'P = '+pval(summary['gnomad']['p']),ha='right',va='top',transform=ax.transAxes,fontsize=7.6)
    ax=f.add_subplot(g[2,1]);ax.scatter(ger.loc[seen,'AA position'],ger.loc[seen,'AF'],s=14,c=TEAL,edgecolor='white',lw=.35);ax.set(yscale='log',xlabel='Germline residue position',ylabel='Joint allele frequency',xlim=(0,1600));pt(ax,'Frequency of observed alleles','F');axstyle(ax)
    ax.tick_params(axis='y', labelsize=7.6)
    export(f,out,'Figure_S3')

def figure_s4s5(w, summary, baseline, thresholds, out):
    f=fig(80);g=f.add_gridspec(1,2,wspace=.13)
    for i,co in enumerate(['A','B']):
        ax=f.add_subplot(g[0,i]); cohort = baseline[baseline.cohort.eq(co)]
        groups = [cohort.loc[cohort.loss, 'position'], cohort.loc[~cohort.loss, 'position']]
        rr=next(r for r in summary['sahoo'] if r['cohort']==co and not r['exclude_986'])
        assert [len(x) for x in groups]==[rr['loss'],rr['normal']],(co,groups)
        for j,(v,c) in enumerate(zip(groups,[TEAL,RUST])):
            v=v.to_numpy(float);off=np.zeros(len(v))
            for val in set(v):
                z=v==val;off[z]=np.linspace(-.11,.11,sum(z)) if sum(z)>1 else 0
            ax.scatter(j+off,v,s=20,c=c,edgecolor='white',linewidth=.4); ax.plot([j-.23,j+.23], [np.median(v)]*2,lw=1.25,c=DARK)
        ax.set(xlim=(-.52,1.52),ylim=(0,1730),ylabel='Germline residue position',yticks=[0,400,800,1200,1600],xticks=[0,1],xticklabels=[f'−7/del(7q)\n(n = {len(groups[0])})',f'Normal karyotype\n(n = {len(groups[1])})']);axstyle(ax);pt(ax,('Consecutive' if i==0 else 'Clinically selected')+' pediatric MDS',co)
        ax.text(.98,.97,'Exact P = '+pval(rr['p']),transform=ax.transAxes,ha='right',va='top',fontsize=7)
    export(f,out,'Figure_S4')
    x,y,cut,_,_,_=roc_inputs(w,summary['roc']);f=fig(140); ax=f.add_subplot(2,1,1)
    groups=[x[y==1],x[y==0]];box(ax,groups,[f'−7/del(7q)\n(n = {len(groups[0])})',f'Comparators\n(n = {len(groups[1])})'],[TEAL,RUST],'Germline residue position',True)
    ax.set(ylim=(0,1730),yticks=[0,400,800,1200,1600]);ax.axhline(cut,color=GRAY,dashes=(3,3),lw=.7);ax.text(.55,cut+24,f'{cut} aa',ha='center',fontsize=7)
    v=summary['roc']['exact_tied_rank_sum_two_sided_p']
    ax.text(.98,.97,f'Exact P = {pval(v)}',transform=ax.transAxes,ha='right',va='top');pt(ax,'Position by chromosome 7 status','A')
    ax=f.add_subplot(2,1,2)
    selected = thresholds.Residue.isin([840, 876, cut])
    ax.scatter(thresholds.loc[~selected, 'Residue'], thresholds.loc[~selected, 'Youden J'], s=9, color=TEAL)
    for residue, score in thresholds.loc[selected, ['Residue','Youden J']].itertuples(index=False, name=None):
        ax.scatter(residue, score, s=16, color=RUST)
        ax.annotate(str(residue), (residue, score), xytext=(0, 5), textcoords='offset points', ha='center', fontsize=7)
    ax.set(xlim=(0, 1600), ylim=(0, .42), xlabel='Germline variant position (amino acid)', ylabel='Empirical Youden index')
    axstyle(ax);pt(ax,'Exploratory cutoffs','B');export(f,out,'Figure_S5')

def figure_s6(phase,reg,out):
    structural=pd.read_csv(reg/'structural_pairs.csv')
    comparisons=pd.read_csv(reg/'structural_summary.csv')
    f=fig(141);g=f.add_gridspec(2,2,height_ratios=[1.15,.95],wspace=.18,hspace=.14)
    ax=f.add_subplot(g[0,0]);ax.plot([0,1600],[0,1600],color=LIGHT,dashes=[3,3],lw=.6)
    phase=phase.copy();phase['label']=phase.phase.map(phase_groups)
    for status,part in phase.groupby('label',sort=False):
        color,symbol,_=PHASE_STYLE[status]
        if status=='Two germline alleles':
            for _,r in part.iterrows():
                p=[numpos(s) for s in r.germline.split(';')];ax.plot(p,[r.somatic_pos]*len(p),color=color,lw=.7,marker='D',ms=3,mfc='white',label=status)
        else:
            ax.scatter(part.germline_pos,part.somatic_pos,s=16,color=color,marker=symbol,edgecolors='white',lw=.3,label=status,alpha=.9)
    ax.set(xlim=(0,1600),ylim=(0,1600),xticks=[0,400,800,1200,1600],yticks=[0,400,800,1200,1600],xlabel='Germline residue position',ylabel='Acquired residue position')
    axstyle(ax);pt(ax,'Phase of reported pairs','A')
    ax2=f.add_subplot(g[0,1]);ax2.axis('off')
    handles,labels=ax.get_legend_handles_labels();ax2.legend(handles,labels,loc='center left',title='Allelic phase',fontsize=8,title_fontsize=8,handlelength=1.2,labelspacing=.8)
    ax=f.add_subplot(g[1,0]);rows=comparisons[comparisons.Analysis.isin(['Pair-level 3D clustering','Background-collapsed'])&comparisons.Metric.isin(['Median Cα distance (Å)','Median heavy-atom distance (Å)'])].copy()
    names=[]
    for j,(_,r) in enumerate(rows.iterrows()):
        names.append(('Pair' if r.Analysis.startswith('Pair') else 'Background')+' · '+('Cα' if 'Cα' in r.Metric else 'heavy atom'))
        for v,dy,c,marker,lab in [('Observed',-.105,TEAL,'o','Observed'),('Null reference',.105,RUST,'s','Reassociated')]:
            ax.scatter(r[v],j+dy,c=c,marker=marker,s=19,label=lab if j==0 else None)
        ax.annotate('P = '+pval(r['Permutation p']),(max(r.Observed,r['Null reference']),j),xytext=(5,0),textcoords='offset points',va='center',fontsize=6.3)
    ax.set(yticks=range(len(names)),yticklabels=names,ylim=(len(rows)-.5,-.5),xlabel='Distance (Å)',xlim=(max(0,min(rows.Observed.min(),rows['Null reference'].min())-4),max(rows.Observed.max(),rows['Null reference'].max())+13));axstyle(ax);pt(ax,'Observed and reassociated distances','B')
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.28),ncol=2,fontsize=6.6)
    ax=f.add_subplot(g[1,1]);s=structural['Primary confidence'].eq('Yes')
    for sel,c,lab,mark in [(~s,GRAY,'Below criteria','o'),(s,TEAL,'Meets criteria','^')]:
        ax.scatter(structural.loc[sel,'Sequence distance aa'],structural.loc[sel,'Cα distance Å'],s=17,c=c,marker=mark,lw=.35,edgecolors='white',label=lab)
    ax.set(xlabel='Sequence separation (residues)',ylabel='Cα distance (Å)');pt(ax,'Sequence and structural separation','C');axstyle(ax)
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.28),ncol=2,fontsize=6.6);export(f,out,'Figure_S6')

def figure_s7(scores,reg,out):
    tissue=pd.read_csv(reg/'tissue_summary.csv');origin=pd.read_csv(reg/'origin_regressions.csv');junction=pd.read_csv(reg/'junction_position.csv')
    rt=pd.read_csv(reg/'rna_tests.csv');bg=pd.read_csv(reg/'rna_backgrounds.csv')
    f=fig(148);g=f.add_gridspec(2,2,wspace=.15,hspace=.19)
    ax=f.add_subplot(g[0,0]);order=['RNA_EXPRESSION','CAGE','SPLICE_JUNCTIONS','SPLICE_SITE_USAGE','POLYADENYLATION','ATAC','DNASE','CHIP_HISTONE','CHIP_TF']
    yy=np.arange(len(order))
    for idx,(key,c,label) in enumerate([('GERMLINE',TEAL,'Germline'),('SOMATIC',RUST,'Acquired')]):
        rr=tissue[tissue.dataset==key].set_index('modality').loc[order];ax.barh(yy+(idx-.5)*.35,rr.n_heme_higher/rr.n_alleles*100,height=.30,color=c,label=label)
    ax.set(yticks=yy,yticklabels=[mod(k) for k in order],ylim=(8.7,-.7),xlim=(0,104),xlabel='Alleles with higher hematopoietic P90 (%)')
    pt(ax,'Tissue specificity of isolated alleles','A');axstyle(ax);ax.grid(axis='y',visible=False);ax.grid(axis='x',lw=.4,color=GRID);ax.legend(loc='upper center',bbox_to_anchor=(.5,-.25),ncol=2,fontsize=6.3)
    ax=f.add_subplot(g[0,1]);d=scores[(scores.modality=='SPLICE_JUNCTIONS')&(scores.dataset=='GERMLINE')]
    ax.scatter(d.AA_position,d.delta_p90_q,s=14,c=TEAL,lw=.3,edgecolors='white');rr=junction[junction.Context.str.contains('delta|relative',case=False)].iloc[0]
    ax.set(xlim=(0,1600),xticks=[0,400,800,1200,1600],xlabel='Germline residue position',ylabel='Hematopoietic − other-tissue P90');pt(ax,'Splice-junction specificity','B');axstyle(ax)
    ax.set_ylim(top=d.delta_p90_q.max()+(d.delta_p90_q.max()-d.delta_p90_q.min())*.20)
    ax.text(.97,.96,f'ρ = {rr.rho:.3f}   BH q = {pval(rr.q_BH6)}',ha='right',va='top',transform=ax.transAxes,fontsize=6.7)
    ax=f.add_subplot(g[1,0]);rr=origin[origin.specification=='all_linear'].set_index('modality').loc[order]
    for y,(_,r) in enumerate(rr.iterrows()):
        c=RUST if r.q_BH9<.05 else TEAL
        ax.plot([r.CI_low,r.CI_high],[y,y],color=c,lw=1.15);ax.scatter(r.Adjusted_beta,y,s=17,c=c,zorder=2)
    ax.set(yticks=yy,yticklabels=[mod(k) for k in order],ylim=(8.7,-.7),xlabel='Adjusted acquired-origin coefficient (95% CI)')
    ax.axvline(0,color=GRAY,lw=.75,dashes=(3,3));pt(ax,'Origin effect by modality','C');axstyle(ax)
    ax=f.add_subplot(g[1,1]);d=bg[bg.subset.str.match(r'^confirmed\d+$')&bg.scope.eq('native_all')].copy()
    d['p']=d.background.map(numpos);d=d.sort_values('p');ax.bar(range(len(d)),d.advantage*1000,color=np.where(d.advantage>0,TEAL,RUST),width=.6)
    ax.axhline(0,color=GRAY,lw=.65);ax.set(xticks=range(len(d)),xticklabels=d.background,ylabel='RNA advantage (score × 1,000)');ax.tick_params(axis='x',labelrotation=35)
    axstyle(ax);pt(ax,'RNA advantage in confirmed cis pairs','D');ax.set_ylim(ax.get_ylim()[0],ax.get_ylim()[1]*1.3)
    rr=rt[(rt.subset==d.subset.iloc[0])&(rt.scope=='native_all')].iloc[0]
    ax.text(.99,.97,f'Sign P = {pval(rr.P_sign)}   BH q = {pval(rr.q_BH8)}',ha='right',va='top',transform=ax.transAxes,fontsize=6.7)
    export(f,out,'Figure_S7')

def figure_s8(reg,out):
    bg=pd.read_csv(reg/'gene_backgrounds.csv');bg=bg[bg.subset.str.startswith('all_modelled')]
    rna=pd.read_csv(reg/'rna_backgrounds.csv'); tests=pd.read_csv(reg/'rna_tests.csv')
    chrom=pd.read_csv(reg/'chromatin_tests.csv');chrom=chrom[chrom.analysis_family=='PRIMARY']
    f=fig(149);g=f.add_gridspec(2,2,wspace=.19,hspace=.23)
    for j,mo in enumerate(['RNA expression','Splice-site usage']):
        sg=g[0,j].subgridspec(1,2,wspace=.02)
        for k,(metric,lab,ref) in enumerate([('median_rescue_to_WT','Proximity score',0),('fraction_closer_to_WT','Fraction closer',.5)]):
            ax=f.add_subplot(sg[0,k]); ax.axvline(ref,ls=(0,(3,3)),color=GRAY,lw=.7);scopes=['native_all','native_hematopoietic'] if j==0 else ['native_all']
            for row,scope in enumerate(scopes):
                v=bg[(bg.scope==scope)&(bg.modality==mo)&(bg.endpoint==metric)].sort_values('background');ys=(row+np.random.default_rng(13+row).uniform(-.14,.14,len(v)))
                ax.scatter(v.value,ys,s=13,color=TEAL if row==0 else RUST,edgecolor='white',lw=.3)
                ax.plot(np.repeat(v.value.median(),2),row+np.array([-.22,.22]),lw=1.35,color=DARK)
            ax.set(yticks=range(len(scopes)),ylim=(len(scopes)-.5,-.5),yticklabels=(['Native','Hematopoietic'][:len(scopes)] if k==0 else []),xlabel=lab)
            if k:ax.set_xlim(-.035,1.035);ax.set_xticks([0,.5,1])
            axstyle(ax)
            if k==0:pt(ax,mo,'AB'[j])
    ax=f.add_subplot(g[1,0]);subsets=list(tests.subset.unique());labels=[]
    for j,sub in enumerate(subsets):
        base=('Confirmed cis' if j==0 else 'Confirmed/reported' if j==1 else 'Plus inferred' if j==2 else 'All modeled')
        n=int(tests.loc[tests.subset==sub,'n_events'].max());labels.append(f'{base} ({n})')
        for k,(scope,c,l) in enumerate([('native_all',TEAL,'Native'),('native_hematopoietic',RUST,'Hematopoietic')]):
            d=rna[(rna.subset==sub)&(rna.scope==scope)].advantage*1000
            yo=j+(k-.5)*.35;ax.scatter(d,np.repeat(yo,len(d)),s=11,color=c,alpha=.8,lw=.25,edgecolors='white',label=l if j==0 else None)
            ax.plot([d.median()]*2,[yo-.105,yo+.105],c=DARK,lw=1)
    ax.axvline(0,ls=(0,(3,3)),color=GRAY,lw=.7);ax.set(yticks=range(len(subsets)),yticklabels=labels,ylim=(len(subsets)-.5,-.5),xlabel='Background RNA advantage (score × 1,000)');axstyle(ax);pt(ax,'RNA advantage by phase','C')
    ax.legend(loc='upper center',bbox_to_anchor=(.48,-.24),ncol=2,fontsize=6.4)
    ax=f.add_subplot(g[1,1]);types=['DNASE','CHIP_HISTONE','CHIP_TF']
    for j,(m,c) in enumerate(zip(types,[TEAL,RUST,PLUM])):
        d=chrom[chrom.output_type==m].sort_values('test_id')
        ax.scatter(np.repeat(j,len(d))+np.linspace(-.19,.19,len(d)),d.q_family_global,c=c,s=15,lw=.3,edgecolor='white')
    ax.axhline(.05,ls=(0,(3,3)),color=GRAY,lw=.75);ax.annotate('0.05',(.99,.05),xycoords=('axes fraction','data'),xytext=(0,2),textcoords='offset points',ha='right',fontsize=6.4)
    ax.set(xticks=[0,1,2],xticklabels=[mod(m) for m in types],ylabel='BH-adjusted P',ylim=(.0,1.05),xlim=(-.5,2.5));axstyle(ax);pt(ax,'Matched chromatin comparisons','D')
    export(f,out,'Figure_S8')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workbook',required=True,type=Path)
    parser.add_argument('--source',required=True,type=Path,help='Directory containing S12_AG_isolated_scores.csv and Figure_1_display_coordinates.csv')
    parser.add_argument('--clinical',required=True,type=Path,help='Recomputed clinical output containing age_rows.csv')
    parser.add_argument('--regulatory',required=True,type=Path,help='Directory with regulatory result CSVs')
    parser.add_argument('--events',required=True,type=Path,help='Audited event-level phase and modeling manifest')
    parser.add_argument('--roc-summary',required=True,type=Path,help='Recomputed pooled ROC summary JSON')
    parser.add_argument('--output',required=True,type=Path);args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    font=font_manager.FontProperties(fname=font_manager.findfont(font_manager.FontProperties(family=['Arial','Liberation Sans','DejaVu Sans']))).get_name()
    plt.rcParams.update({'font.family':font,'mathtext.fontset':'custom','mathtext.rm':font,'mathtext.it':font + ':italic','mathtext.bf':font + ':bold','mathtext.cal':font,'font.size':7.2,'font.weight':'normal','axes.labelsize':7.3,'xtick.labelsize':6.7,'ytick.labelsize':6.7,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','svg.hashsalt':'SAMD9L-letter','axes.linewidth':.65,'axes.edgecolor':GRAY,'axes.labelcolor':DARK,'text.color':DARK,'xtick.color':DARK,'ytick.color':DARK,'figure.facecolor':'white','axes.facecolor':'white','legend.frameon':False,'savefig.facecolor':'white'})
    w,a,ev,land,cl,mx,missing=load_data(args)
    if missing: raise ValueError(f'Unmapped acquired event(s) {missing}')
    if len(ev)!=sum(land.loc[land.origin=='Acquired','count']):raise ValueError('Acquired observations do not reconcile with S2B')
    phase=pd.read_csv(args.events);assert set(ev['Event ID'])==set(phase.event),'Event manifest and S2B mismatch'
    scores=pd.read_csv(args.source/'S12_AG_isolated_scores.csv')
    normalization = pd.read_csv(Path(__file__).resolve().parents[2] / 'data/protein_normalizations.csv').set_index('variant_id')
    corrected = scores.variant_id_project.map(normalization.residue)
    scores.loc[corrected.notna(), 'AA_position'] = corrected.dropna().astype(int)
    summary=json.loads((args.clinical/'clinical_summary.json').read_text())
    summary['roc']=json.loads(args.roc_summary.read_text())
    clonal=pd.read_csv(args.clinical/'clonal.csv')
    ages=pd.read_csv(args.clinical/'age_rows.csv')
    for final,expected in [(cl,'clinical'),(mx,'mixed')]:
        d=ages[ages.endpoint==expected].set_index('patient').sort_index()
        a1=final.set_index('Patient ID').sort_index()
        if not d.index.equals(a1.index) or not np.allclose(d.age,a1.age) or not np.allclose(d.position,a1.position):
            raise ValueError(f'{expected} age roster in S1 differs from recomputed clinical output')
    figure1(w,land,cl,summary['roc'],args.output);figure2(a,phase,scores,args.regulatory,args.output)
    figure_s1(a,land,cl,summary,clonal,args.output);figure_s2(mx,args.output)
    figure_s3(w,ages,summary,args.output)
    baseline = pd.read_csv(args.clinical/'sahoo_baseline.csv')
    thresholds = pd.read_csv(args.clinical/'letter_roc_105/thresholds.csv')
    figure_s4s5(w, summary, baseline, thresholds, args.output)
    figure_s6(phase,args.regulatory,args.output);figure_s7(scores,args.regulatory,args.output);figure_s8(args.regulatory,args.output)
    stats={'font':font,'cohort':len(a),'clinical_age':len(cl),'broader_age':len(mx),'phase_n':len(phase),'structural_n':len(pd.read_csv(args.regulatory/'structural_pairs.csv')),'landscape_counts':land.groupby('origin')['count'].agg(['size','sum']).to_dict('index')}
    (args.output/'figure_manifest.json').write_text(json.dumps(stats,indent=2))
if __name__=='__main__':main()
