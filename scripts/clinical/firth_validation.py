"""Firth and FLIC estimators, profile intervals and design matrices from the archived Letter implementation."""
from pathlib import Path
import argparse,csv,json,math,time,hashlib,sys
from collections import Counter,defaultdict
import numpy as np
import scipy
from scipy.special import expit
from scipy.optimize import minimize,brentq
from scipy.linalg import null_space
from scipy.stats import chi2,rankdata

SEED=20261110
SOURCE_ORDER=['Other St Jude','Sahoo A','Sahoo B','Nagata','Other literature']

def firth_parts(beta,X,y):
    eta=X@beta;p=expit(eta);w=p*(1-p)
    info=X.T@(w[:,None]*X)
    sign,logdet=np.linalg.slogdet(info)
    if sign<=0 or not np.isfinite(logdet):return -np.inf,None,None
    try:inv=np.linalg.inv(info)
    except np.linalg.LinAlgError:return -np.inf,None,None
    h=w*np.einsum('ij,jk,ik->i',X,inv,X)
    score=X.T@(y-p+h*(.5-p))
    pll=float(np.sum(y*eta-np.logaddexp(0,eta))+.5*logdet)
    return pll,score,inv

def fit_firth(X,y,start=None):
    beta=np.zeros(X.shape[1]) if start is None else start.copy();converged=False
    for it in range(150):
        ll,u,inv=firth_parts(beta,X,y)
        if u is None:break
        if np.max(np.abs(u))<1e-7:converged=True;break
        delta=inv@u;maxstep=np.max(np.abs(delta))
        if maxstep>4:delta*=4/maxstep
        scale=1.
        for step in range(35):
            trial=beta+scale*delta;newll=firth_parts(trial,X,y)[0]
            if np.isfinite(newll) and newll>=ll-1e-11:break
            scale*=.5
        else:break
        beta=trial
    ll,u,inv=firth_parts(beta,X,y)
    if not converged or u is None or np.max(np.abs(u))>2e-6:
        def fun(b):
            l,s,_=firth_parts(b,X,y)
            return (-l,-s) if s is not None else (1e100,np.zeros(X.shape[1]))
        opt=minimize(fun,beta,jac=True,method='BFGS',options={'gtol':1e-8,'maxiter':500})
        beta=opt.x;ll,u,inv=firth_parts(beta,X,y)
        converged=u is not None and np.max(np.abs(u))<2e-5
    if not converged:raise RuntimeError(f'Firth did not converge; score={None if u is None else np.max(abs(u))}')
    return dict(beta=beta,pll=ll,inv=inv,score_norm=float(np.max(abs(u))))

def design(rows,features,source_levels=None,check_rank=True):
    age=np.array([r['age_reported_mixed'] if r['age_reported_mixed'] is not None else np.nan for r in rows])
    adult=np.array([r['age_group']=='Adult' for r in rows],float)
    pos=np.array([(r['position']-1000)/100 for r in rows])
    cols=[np.ones(len(rows))];names=['intercept']
    if 'position' in features:cols.append(pos);names.append('position_per100')
    if 'cut986' in features:cols.append(np.array([r['position']>=986 for r in rows],float));names.append('position_ge986')
    if 'age_group' in features:cols.append(adult);names.append('adult')
    if 'age_continuous' in features:cols.append(np.log1p(age));names.append('log1p_age')
    if 'interaction' in features:cols.append(pos*adult);names.append('position_x_adult')
    levels=[]
    if 'source' in features:
        levels=source_levels if source_levels is not None else [s for s in SOURCE_ORDER if any(r['source_group']==s for r in rows)]
        for s in levels[1:]:cols.append(np.array([r['source_group']==s for r in rows],float));names.append('source_'+s)
    X=np.column_stack(cols)
    if not np.isfinite(X).all() or (check_rank and np.linalg.matrix_rank(X)!=X.shape[1]):raise ValueError('Nonfinite or rank-deficient design')
    return X,names,levels

def load_rows(path):
    rows=list(csv.DictReader(open(path,encoding='utf-8-sig')))
    for r in rows:
        for k in ['position','eligible','complex_karyotype']:
            r[k]=int(float(r[k])) if r[k] else None
        for k in ['age_reported_mixed','clinical_context_age','hematologic_age']:
            r[k]=float(r[k]) if r[k] else None
    return [r for r in rows if r['eligible']]

def auc(y,score):
    y=np.array(y,bool);n1=int(y.sum());n0=len(y)-n1
    return float((rankdata(score)[y].sum()-n1*(n1+1)/2)/(n1*n0)) if n1 and n0 else None

def performance(y,eta):
    p=expit(eta)
    return dict(auc=auc(y,eta),brier=float(np.mean((p-y)**2)),logloss=float(np.mean(np.logaddexp(0,eta)-y*eta)))

def calibration(y,eta):
    X=np.column_stack([np.ones(len(eta)),eta])
    def fun(b):
        z=X@b
        return float(np.sum(np.logaddexp(0,z)-y*z)),X.T@(expit(z)-y)
    opt=minimize(fun,np.array([0.,1.]),jac=True,method='BFGS',options={'gtol':1e-8,'maxiter':1000})
    assert np.max(abs(fun(opt.x)[1]))<1e-5
    citl=brentq(lambda a:np.mean(expit(eta+a))-np.mean(y),-80,80)
    return dict(calibration_in_large=float(citl),calibration_slope=float(opt.x[1]),free_intercept=float(opt.x[0]))

def recalibrate_intercept(eta,y):
    # Training-only FLIC intercept correction. Slopes remain the Firth estimates.
    shift=brentq(lambda d:np.mean(expit(eta+d))-np.mean(y),-80,80)
    return float(shift)

def fit_model(rows,features):
    y=np.array([r['cytogenetic_class']=='loss7' for r in rows],float)
    if len(np.unique(y))<2:raise ValueError('Only one outcome class')
    X,names,levels=design(rows,features);fit=fit_firth(X,y)
    b=fit['beta'];shift=recalibrate_intercept(X@b,y)
    return dict(**fit,X=X,y=y,names=names,source_levels=levels,shift=shift,n=len(rows),events=int(sum(y)),features=features,ids=[r['patient_id'] for r in rows],performance=performance(y,X@b+shift))

def constrained_pll(X,y,beta,c,value):
    c=np.array(c,float);Z=null_space(c.reshape(1,-1));base=c*value/np.dot(c,c)
    init=Z.T@(beta-base)
    def fun(g):
        ll,u,_=firth_parts(base+Z@g,X,y)
        return (-ll,-Z.T@u) if u is not None else (1e100,np.zeros(Z.shape[1]))
    opt=minimize(fun,init,jac=True,method='BFGS',options={'gtol':1e-8,'maxiter':1000})
    ll,u,_=firth_parts(base+Z@opt.x,X,y)
    if u is None or max(abs(Z.T@u),default=0)>1e-4:raise RuntimeError('Constrained profile did not converge')
    return float(ll)

def profile(fit,c):
    c=np.array(c,float);b=fit['beta'];est=float(c@b);se=float(np.sqrt(c@fit['inv']@c));limit=chi2.ppf(.95,1)
    def dev(v):return 2*(fit['pll']-constrained_pll(fit['X'],fit['y'],b,c,v))
    bounds=[]
    for sign in [-1,1]:
        width=max(.1,se)
        for k in range(20):
            v=est+sign*width
            if dev(v)>limit:break
            width*=1.7
        else:raise RuntimeError('Could not bracket profile interval')
        bounds.append(float(brentq(lambda v:dev(v)-limit,min(est,v),max(est,v),xtol=1e-7)))
    return dict(beta=est,OR=math.exp(est),profile95_OR=[math.exp(v) for v in bounds],nominal_profile_PLR_p=float(chi2.sf(max(0,dev(0)),1)),profile95_beta=bounds)

def contrast(fit,name):
    c=np.zeros(len(fit['names']))
    if name=='adult_position':
        c[fit['names'].index('position_per100')]=1;c[fit['names'].index('position_x_adult')]=1
    elif name=='pediatric_position':c[fit['names'].index('position_per100')]=1
    else:c[fit['names'].index(name)]=1
    return c
