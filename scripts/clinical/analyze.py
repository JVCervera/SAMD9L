"""Exact tied rank-sum test by integer dynamic programming, retained from the archived Letter implementation."""
from pathlib import Path
import csv,json,math,argparse,hashlib
from collections import defaultdict,Counter
import numpy as np
from scipy.stats import rankdata,mannwhitneyu

def rank_distribution(values,n):
    ranks=np.rint(rankdata(values)*2).astype(int)
    # Integer DP coefficient of z^n x^rank_sum; identical observations remain
    # separate exchangeable individuals, so ties carry their multiplicity.
    cap=int(sum(sorted(ranks,reverse=True)[:n]))
    dp=np.zeros((n+1,cap+1),dtype=object);dp[0,0]=1
    for i,r in enumerate(ranks):
        for k in range(min(i+1,n),0,-1):dp[k,r:]+=dp[k-1,:cap+1-r]
    dist={int(i):int(v) for i,v in enumerate(dp[n]) if v}
    assert sum(dist.values())==math.comb(len(values),n)
    return ranks,dist

def exact_rank_test(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    if not len(a) or not len(b):return None
    vals=np.r_[a,b];n=len(a)
    # Smaller group makes DP cheaper; deviation test is invariant to complement.
    if len(a)>len(b): vals=np.r_[b,a];n=len(b)
    ranks,dist=rank_distribution(vals,n)
    obs=int(ranks[:n].sum());totalrank=int(ranks.sum());N=len(vals)
    dev=abs(obs*N-n*totalrank)
    extreme=sum(v for s,v in dist.items() if abs(s*N-n*totalrank)>=dev)
    total=math.comb(N,n)
    U=float(mannwhitneyu(a,b,method='asymptotic').statistic)
    return dict(n1=len(a),n0=len(b),median1=float(np.median(a)),median0=float(np.median(b)),q1_1=float(np.quantile(a,.25)),q3_1=float(np.quantile(a,.75)),q1_0=float(np.quantile(b,.25)),q3_0=float(np.quantile(b,.75)),median_difference=float(np.median(a)-np.median(b)),p_exact=extreme/total,extreme_assignments=str(extreme),total_assignments=str(total),U=U,probability_higher=(U/(len(a)*len(b))),cliffs_delta=(2*U/(len(a)*len(b))-1))
