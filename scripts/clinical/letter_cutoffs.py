"""Rebuild the all-carrier validation and the Youden panel of Supplementary Figure S5."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pooled_roc import auc, design_matrix, held_out, groups, youden

def calculate(path, output):
    data = pd.read_csv(path)
    assert set(data.age_group) == {"adult", "pediatric"}
    adult = data.age_group.eq("adult").to_numpy()
    source = data.source_group.to_numpy()
    position = data.position_aa.to_numpy()
    y = data.chromosome_7_loss.to_numpy()
    blocks = groups(data.protein_variant.to_numpy())
    models = {
        "position": np.column_stack([np.ones(len(y)), (position - 986) / 100]),
        "position_age_source": design_matrix(adult, source, position),
        "age_source": design_matrix(adult, source),
    }
    output.mkdir(parents=True, exist_ok=True)
    results = {}
    predictions = data[['patient_id','protein_variant','chromosome_7_loss']].copy()
    for name, matrix in models.items():
        predictions[name] = held_out(y,matrix,blocks)
        results[name] = auc(y,predictions[name])
    predictions.to_csv(output/'held_out_105.csv',index=False)
    pd.Series(results,name='AUC').to_csv(output/'auc_105.csv')
    records=[]
    for threshold in sorted(set(position)):
        positive = position >= threshold
        if min(positive.sum(), (~positive).sum()) < 10: continue
        sensitivity=positive[y==1].mean()
        specificity=(~positive[y==0]).mean()
        records.append([threshold,sensitivity,specificity,sensitivity+specificity-1])
    scores=pd.DataFrame(records,columns=['Residue','Sensitivity','Specificity','Youden J'])
    scores.to_csv(output/'thresholds.csv',index=False)
    keep=position!=986
    omitted=[len(y[keep]),int(y[keep].sum()),auc(y[keep],position[keep]),youden(y[keep],position[keep])['position_aa']]
    pd.DataFrame([omitted],columns=['N','loss','AUC','optimal_residue']).to_csv(output/'omit_986.csv',index=False)
    assert np.allclose(list(results.values()),[.674,.632,.575],atol=.0005)
    assert omitted[3]==840 and round(omitted[2],3)==.736
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'Arial'],
                         'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none', 'svg.hashsalt': 'SAMD9L-S5B'})
    fig, ax = plt.subplots(figsize=(6.1, 3.31), layout='constrained')
    highlights = [840, 876, 986]
    ordinary = scores.loc[~scores.Residue.isin(highlights)]
    ax.scatter(ordinary.Residue, ordinary['Youden J'], s=8, color='#35586c', linewidths=0)
    for residue, offset in [(840, (-13, -9)), (876, (0, 6)), (986, (5, 4))]:
        value = scores.set_index('Residue').at[residue, 'Youden J']
        ax.scatter(residue, value, s=16, color='#bd6228', linewidths=0)
        ax.annotate(str(residue), (residue, value), xytext=offset, textcoords='offset points',
                    ha='center', fontsize=7)
    ax.set(xlabel='Germline variant position (amino acid)', ylabel='Empirical Youden index',
           xlim=(0, 1600), ylim=(0, 0.42))
    ax.spines[['right','top']].set_visible(False)
    for fmt in ('pdf', 'svg', 'png'):
        fig.savefig(output / f'Figure_S5B_Youden.{fmt}', dpi=600)
    plt.close(fig)
if __name__=='__main__':
    root=Path(__file__).resolve().parents[2]
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    calculate(root/'data/pooled_roc_input.csv',parser.parse_args().output)
