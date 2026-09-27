#!/usr/bin/env python3
"""C2R manuscript-grade robustness and incremental-information audit.

All calculations use the frozen C2 integrated table and the C2R preregistration.
No protein RAW files are downloaded or reprocessed.
"""
from pathlib import Path
import json, math, sys, warnings
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.model_selection import KFold

warnings.filterwarnings("ignore")
ROOT=Path(__file__).resolve().parents[2]; CM=ROOT/"computational_manuscript"
OUT=CM/"results/c2r"; REP=CM/"reports/c2r"; OUT.mkdir(parents=True,exist_ok=True); REP.mkdir(parents=True,exist_ok=True)
SEED=270927; B=10000

def bh(p):
    p=np.asarray(p,float); q=np.full(p.shape,np.nan); ok=np.isfinite(p)
    if ok.any(): q[ok]=multipletests(p[ok],method="fdr_bh")[1]
    return q

def spearman(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); ok=np.isfinite(x)&np.isfinite(y)
    if ok.sum()<4: return np.nan
    return float(stats.spearmanr(x[ok],y[ok]).statistic)

def pearson(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); ok=np.isfinite(x)&np.isfinite(y)
    if ok.sum()<4: return np.nan
    return float(stats.pearsonr(x[ok],y[ok]).statistic)

def kendall(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); ok=np.isfinite(x)&np.isfinite(y)
    if ok.sum()<4: return np.nan
    return float(stats.kendalltau(x[ok],y[ok]).statistic)

def bicor(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); ok=np.isfinite(x)&np.isfinite(y); x=x[ok]; y=y[ok]
    if len(x)<4: return np.nan
    mx,my=np.median(x),np.median(y); madx=np.median(np.abs(x-mx)); mady=np.median(np.abs(y-my))
    if madx==0 or mady==0: return pearson(x,y)
    ux=(x-mx)/(9*madx); uy=(y-my)/(9*mady); keep=(np.abs(ux)<1)&(np.abs(uy)<1)
    if keep.sum()<4: return np.nan
    wx=(1-ux[keep]**2)**2; wy=(1-uy[keep]**2)**2; xx=x[keep]-mx; yy=y[keep]-my
    return float(np.sum(xx*yy*wx*wy)/np.sqrt(np.sum(xx**2*wx**2)*np.sum(yy**2*wy**2)))

def bootstrap_corr(x,y,B=B,seed=SEED):
    x=np.asarray(x,float); y=np.asarray(y,float); ok=np.isfinite(x)&np.isfinite(y); x=x[ok]; y=y[ok]; n=len(x)
    if n<4: return np.nan,np.nan,np.nan,np.nan,np.array([])
    rx=stats.rankdata(x,method="average"); ry=stats.rankdata(y,method="average"); rg=np.random.default_rng(seed); vals=[]; chunk=250
    for st in range(0,B,chunk):
        k=min(chunk,B-st); idx=rg.integers(0,n,size=(k,n),dtype=np.int32); a=rx[idx]; b=ry[idx]; a-=a.mean(1,keepdims=True); b-=b.mean(1,keepdims=True); vals.extend(((a*b).sum(1)/np.sqrt((a*a).sum(1)*(b*b).sum(1))).tolist())
    vals=np.asarray(vals); return float(np.nanmedian(vals)),float(np.nanquantile(vals,.025)),float(np.nanquantile(vals,.975)),float(2*min(np.mean(vals<=0),np.mean(vals>=0))),vals

def rank_partial(x,y,c):
    x=np.asarray(x,float); y=np.asarray(y,float); c=np.asarray(c,float); ok=np.isfinite(x)&np.isfinite(y)&np.isfinite(c)
    if ok.sum()<5: return np.nan
    xr=stats.rankdata(x[ok],method="average"); yr=stats.rankdata(y[ok],method="average"); cr=stats.rankdata(c[ok],method="average")
    X=np.column_stack([np.ones(len(cr)),cr]); rx=xr-X@np.linalg.lstsq(X,xr,rcond=None)[0]; ry=yr-X@np.linalg.lstsq(X,yr,rcond=None)[0]
    return float(stats.pearsonr(rx,ry).statistic)

def bootstrap_partial(x,y,c,B=B,seed=SEED):
    x=np.asarray(x,float); y=np.asarray(y,float); c=np.asarray(c,float); ok=np.isfinite(x)&np.isfinite(y)&np.isfinite(c); x=x[ok]; y=y[ok]; c=c[ok]; n=len(x)
    if n<5: return np.nan,np.nan,np.nan,np.nan,np.array([])
    xr=stats.rankdata(x,method="average"); yr=stats.rankdata(y,method="average"); cr=stats.rankdata(c,method="average"); X=np.column_stack([np.ones(n),cr]); rx=xr-X@np.linalg.lstsq(X,xr,rcond=None)[0]; ry=yr-X@np.linalg.lstsq(X,yr,rcond=None)[0]
    # Fixed rank-residual bootstrap is the preregistered paired gene bootstrap implementation.
    rg=np.random.default_rng(seed); vals=[]; chunk=250
    for st in range(0,B,chunk):
        k=min(chunk,B-st); idx=rg.integers(0,n,size=(k,n),dtype=np.int32); a=rx[idx]; b=ry[idx]; a-=a.mean(1,keepdims=True); b-=b.mean(1,keepdims=True); vals.extend(((a*b).sum(1)/np.sqrt((a*a).sum(1)*(b*b).sum(1))).tolist())
    vals=np.asarray(vals); return float(rank_partial(x,y,c)),float(np.nanquantile(vals,.025)),float(np.nanquantile(vals,.975)),float(2*min(np.mean(vals<=0),np.mean(vals>=0))),vals

def mean_bootstrap(vals,B=10000,seed=1):
    vals=np.asarray(vals,float); vals=vals[np.isfinite(vals)]; n=len(vals); rg=np.random.default_rng(seed); idx=rg.integers(0,n,size=(B,n),dtype=np.int32); bs=vals[idx].mean(1); return float(np.mean(vals)),float(np.quantile(bs,.025)),float(np.quantile(bs,.975)),float(2*min(np.mean(bs<=0),np.mean(bs>=0))),bs

def ols_slope(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); ok=np.isfinite(x)&np.isfinite(y); X=np.column_stack([np.ones(ok.sum()),x[ok]]); coef=np.linalg.lstsq(X,y[ok],rcond=None)[0]; return coef,ok

def read_inputs():
    m=pd.read_csv(OUT.parent/"c2/CROSSOMIC_MASTER_TABLE.tsv",sep="\t")
    dep=pd.read_csv(OUT.parent/"c2/PXD054330_DIFFERENTIAL_PROTEOME.tsv",sep="\t")
    # Protein control abundance is not in the C2 master; merge it by the unique C2 gene symbol.
    dep=dep[["gene","mean_control"]].drop_duplicates("gene").rename(columns={"gene":"gene_symbol"})
    m=m.merge(dep,on="gene_symbol",how="left")
    ann=pd.read_csv(ROOT/"data/raw/GSE323164/GSE323164_1_genes_fpkm_expression.txt.gz",sep="\t",usecols=["gene_name","chr"])
    ann=ann.drop_duplicates("gene_name"); chrom=dict(zip(ann.gene_name.astype(str),ann.chr.astype(str)))
    m["chromosome"]=m.gene_symbol.astype(str).map(chrom)
    return m

def residual_analysis(m):
    ok=m[["RNA_6h_log2FC","Ribo_6h_log2FC","PXD054330_24h_protein_log2FC"]].notna().all(1); d=m.loc[ok].copy(); coef,_=ols_slope(d.RNA_6h_log2FC,d.Ribo_6h_log2FC); d["Ribo6_residual"]=d.Ribo_6h_log2FC-(coef[0]+coef[1]*d.RNA_6h_log2FC); med,lo,hi,p,vals=bootstrap_corr(d.Ribo6_residual,d.PXD054330_24h_protein_log2FC,B=B,seed=270927); rho=spearman(d.Ribo6_residual,d.PXD054330_24h_protein_log2FC); out=pd.DataFrame([{"n":len(d),"RNA6_to_Ribo6_intercept":coef[0],"RNA6_to_Ribo6_slope":coef[1],"residual_protein_spearman_rho":rho,"bootstrap_median_rho":med,"CI_low":lo,"CI_high":hi,"empirical_P":p,"PASS":bool(rho>0 and lo>0 and p<.05)}]); out.to_csv(OUT/"RNA_ADJUSTED_RIBO_RESULTS.tsv",sep="\t",index=False); return d,out

def partial_analysis(d):
    a,alo,ahi,ap,_=bootstrap_partial(d.Ribo_6h_log2FC,d.PXD054330_24h_protein_log2FC,d.RNA_6h_log2FC,B=B,seed=270928); b,blo,bhi,bp,_=bootstrap_partial(d.RNA_6h_log2FC,d.PXD054330_24h_protein_log2FC,d.Ribo_6h_log2FC,B=B,seed=270929); out=pd.DataFrame([{"association":"Ribo6_Protein24_adjusted_RNA6","n":len(d),"partial_rho":a,"CI_low":alo,"CI_high":ahi,"empirical_P":ap},{"association":"RNA6_Protein24_adjusted_Ribo6","n":len(d),"partial_rho":b,"CI_low":blo,"CI_high":bhi,"empirical_P":bp}]); out.to_csv(OUT/"PARTIAL_ASSOCIATION.tsv",sep="\t",index=False); return out

def standardize_train(Xtr,Xte):
    mu=np.nanmean(Xtr,0); sd=np.nanstd(Xtr,0,ddof=1); sd[~np.isfinite(sd)|(sd==0)]=1; return (Xtr-mu)/sd,(Xte-mu)/sd

def fit_predict(X,y,tr,te):
    Xtr,Xte=standardize_train(X[tr],X[te]); ytr=np.asarray(y[tr],float); beta=np.linalg.lstsq(np.column_stack([np.ones(len(tr)),Xtr]),ytr,rcond=None)[0]; return np.column_stack([np.ones(len(te)),Xte])@beta

def metric_row(y,pred):
    rho=spearman(y,pred); sse=np.sum((y-pred)**2); sst=np.sum((y-np.mean(y))**2); return rho,float(1-sse/sst),float(np.sqrt(np.mean((y-pred)**2)))

def make_folds(n):
    rows=[]
    for rep in range(10):
        kf=KFold(10,shuffle=True,random_state=270929+rep)
        for fold,(tr,te) in enumerate(kf.split(np.arange(n))): rows.append((rep,fold,tr,te))
    return rows

def nested_cv(m):
    cols=["RNA_6h_log2FC","Ribo_6h_log2FC","TE_6h_log2FC","PXD054330_24h_protein_log2FC"]; ok=m[cols].notna().all(1); d=m.loc[ok].reset_index(drop=True); y=d.PXD054330_24h_protein_log2FC.to_numpy(float); Xs={"A":d[["RNA_6h_log2FC"]].to_numpy(float),"B":d[["RNA_6h_log2FC","Ribo_6h_log2FC"]].to_numpy(float),"C":d[["RNA_6h_log2FC","TE_6h_log2FC"]].to_numpy(float),"D":d[["RNA_6h_log2FC","Ribo_6h_log2FC","TE_6h_log2FC"]].to_numpy(float)}; rows=[]; folds=make_folds(len(d))
    for rep,fold,tr,te in folds:
        for name,X in Xs.items():
            pred=fit_predict(X,y,tr,te); rho,r2,rmse=metric_row(y[te],pred); rows.append({"repeat":rep,"fold":fold,"model":name,"n_test":len(te),"out_of_sample_Spearman_rho":rho,"out_of_sample_R2":r2,"out_of_sample_RMSE":rmse})
    out=pd.DataFrame(rows); out.to_csv(OUT/"NESTED_MODEL_CV.tsv",sep="\t",index=False)
    piv=out.pivot_table(index=["repeat","fold"],columns="model",values="out_of_sample_Spearman_rho").reset_index(); diff=(piv.B-piv.A).to_numpy(float); mean,lo,hi,p,bs=mean_bootstrap(diff,B=B,seed=270930); summary=out.groupby("model")[['out_of_sample_Spearman_rho','out_of_sample_R2','out_of_sample_RMSE']].agg(['mean','std']).reset_index(); summary.columns=['_'.join([x for x in c if x]) for c in summary.columns]; summary.to_csv(OUT/"NESTED_MODEL_CV_SUMMARY.tsv",sep="\t",index=False); delta=pd.DataFrame([{"comparison":"Model_B_minus_Model_A","n_fold_repeat":len(diff),"Delta_CV_rho":mean,"CI_low":lo,"CI_high":hi,"empirical_P":p,"Model_B_out_of_sample_better":bool(mean>0 and lo>0)}]); delta.to_csv(OUT/"NESTED_MODEL_CV_DELTA.tsv",sep="\t",index=False); return d,folds,out,summary,delta

def permutation_null(d,folds,observed):
    cols=["RNA_6h_log2FC","Ribo_6h_log2FC","TE_6h_log2FC","PXD054330_24h_protein_log2FC"]; d=d.reset_index(drop=True); RNA=d.RNA_6h_log2FC.to_numpy(float); Ribo=d.Ribo_6h_log2FC.to_numpy(float); y=d.PXD054330_24h_protein_log2FC.to_numpy(float); base=d.baseline_RNA_control.to_numpy(float); prot=d.mean_control.to_numpy(float); coef,_=ols_slope(RNA,Ribo); resid=Ribo-(coef[0]+coef[1]*RNA)
    rb=pd.qcut(RNA,5,labels=False,duplicates='drop'); bb=pd.qcut(base,4,labels=False,duplicates='drop'); pb=pd.qcut(prot,4,labels=False,duplicates='drop'); groups=pd.Series([f'{a}_{b}_{c}' for a,b,c in zip(rb,bb,pb)]); gidx={g:np.where(groups.to_numpy()==g)[0] for g in groups.unique()}
    # One fixed 10-fold split (first repeat) is used for the 10,000-permutation null for computational tractability.
    fixed=[x for x in folds if x[0]==0]; A=d[["RNA_6h_log2FC"]].to_numpy(float); rg=np.random.default_rng(270931); null=[]
    for k in range(B):
        pr=resid.copy()
        for idx in gidx.values():
            if len(idx)>1: pr[idx]=pr[rg.permutation(idx)]
        Rnull=coef[0]+coef[1]*RNA+pr; Xb=np.column_stack([RNA,Rnull]); ra=[]; rbv=[]
        for _,_,tr,te in fixed:
            pa=fit_predict(A,y,tr,te); pbv=fit_predict(Xb,y,tr,te); ra.append(spearman(y[te],pa)); rbv.append(spearman(y[te],pbv))
        null.append(np.nanmean(rbv)-np.nanmean(ra))
    null=np.asarray(null); p=float((1+np.sum(null>=observed.Delta_CV_rho.iloc[0]))/(len(null)+1)); out=pd.DataFrame({"permutation":np.arange(1,B+1),"null_Delta_CV_rho":null}); out["observed_Delta_CV_rho"]=float(observed.Delta_CV_rho.iloc[0]); out["empirical_upper_tail_P"]=p; out.to_csv(OUT/"PERMUTATION_NULL.tsv",sep="\t",index=False); return p,null

def detectability_set(m):
    f=CM/"data/external/PXD059451/SupplementaryTable1_protein_List_integrity_and_consistency.xlsx"; mapping=pd.read_csv(CM/"metadata/c1m/PROTEIN_GENE_ID_MAP_REAL.tsv",sep="\t"); tok={}
    for _,r in mapping.iterrows():
        ens=str(r.ensembl_gene_id)
        if ';' in ens: continue
        for acc in str(r.original_protein_identifier).split(';'): tok.setdefault(acc.split('-')[0],set()).add(ens)
    # PXD059451 uses a different protein-group list from PXD054330. Its
    # Protein Name column provides HGNC-like symbols (for example TMA7B_HUMAN)
    # and is therefore the reproducible cross-resource mapping fallback.
    sym={}
    for _,r in m[["gene_symbol","ensembl_gene_id"]].dropna().iterrows():
        for s in str(r.gene_symbol).split(';'):
            sym.setdefault(s.strip().upper(),set()).add(str(r.ensembl_gene_id))
    x=pd.read_excel(f,header=None); rows=[]
    for _,r in x.iloc[2:].iterrows():
        group=str(r.iloc[0]);
        if group=='nan': continue
        # Columns 3--8 are the six platform-specific consistency values;
        # column 9 is their workbook-provided average.
        vals=pd.to_numeric(pd.Series([r.iloc[i] for i in [3,4,5,6,7,8]]),errors='coerce'); avg=pd.to_numeric(pd.Series([r.iloc[9]]),errors='coerce').iloc[0]; high=bool(np.isfinite(avg) and avg>=50 and (vals>=50).sum()>=2); ens=set();
        for a in group.split(';'): ens |= tok.get(a.split('-')[0],set())
        if len(ens)==0:
            for name in str(r.iloc[1]).split(';'):
                s=name.strip().upper()
                if s.endswith('_HUMAN'): s=s[:-6]
                ens |= sym.get(s,set())
        for e in ens: rows.append({"ensembl_gene_id":e,"average_consistency":avg,"n_platforms_ge50":int((vals>=50).sum()),"high_confidence":high})
    a=pd.DataFrame(rows).drop_duplicates("ensembl_gene_id"); high=set(a.loc[a.high_confidence,"ensembl_gene_id"]); return a,high

def detectability_analysis(m,d,folds):
    ann,high=detectability_set(m); write=ann; write.to_csv(OUT/"PXD059451_DETECTABILITY_ANNOTATION.tsv",sep="\t",index=False); d2=d[d.ensembl_gene_id.isin(high)].copy();
    if len(d2)>=100:
        rho=spearman(d2.RNA_6h_log2FC,d2.PXD054330_24h_protein_log2FC); rb=spearman(d2.Ribo_6h_log2FC,d2.PXD054330_24h_protein_log2FC); delta=rb-rho; status='YES'
    else: rho=rb=delta=np.nan; status='EXTERNAL_DETECTABILITY_NOT_AVAILABLE'
    out=pd.DataFrame([{"dataset":"PXD059451","high_confidence_genes":len(high),"mapped_primary_genes":len(d2),"external_detectability_used":status,"RNA6_rho":rho,"Ribo6_rho":rb,"Delta_rho":delta,"abundance_rank_available":False}]); out.to_csv(OUT/"DETECTABILITY_SENSITIVITY.tsv",sep="\t",index=False); return out,ann,high

def stratification(m):
    # The frozen C2 master table does not contain a Ribo standard-error column.
    # Use the available Ribo6 FDR only as a labelled statistical-precision proxy;
    # do not relabel it as a standard error or infer abundance coverage.
    precision_col = "SE_Ribo_6h" if "SE_Ribo_6h" in m.columns else "Ribo_6h_FDR"
    precision_source = "Ribo6_standard_error" if precision_col == "SE_Ribo_6h" else "Ribo6_FDR_precision_proxy"
    required=["RNA_6h_log2FC","Ribo_6h_log2FC","PXD054330_24h_protein_log2FC","baseline_RNA_control",precision_col,"mean_control"]
    d=m[m[required].notna().all(1)].copy(); specs=[("baseline_RNA_quartile","baseline_RNA_control"),("Ribo_precision_quartile",precision_col),("baseline_protein_quartile","mean_control"),("protein_effect_abs_quartile",None)]; rows=[]
    for label,col in specs:
        v=np.abs(d.PXD054330_24h_protein_log2FC) if col is None else d[col]; d['_strat']=pd.qcut(v,4,labels=['Q1','Q2','Q3','Q4'],duplicates='drop')
        for q,a in d.groupby('_strat',observed=True):
            r=spearman(a.RNA_6h_log2FC,a.PXD054330_24h_protein_log2FC); rb=spearman(a.Ribo_6h_log2FC,a.PXD054330_24h_protein_log2FC); rows.append({"stratification":label,"stratum":str(q),"n":len(a),"RNA6_rho":r,"Ribo6_rho":rb,"Delta_rho":rb-r,"precision_source":precision_source})
    out=pd.DataFrame(rows); out.to_csv(OUT/"ABUNDANCE_STRATIFIED_RESULTS.tsv",sep='\t',index=False); out[out.stratification=='protein_effect_abs_quartile'].to_csv(OUT/"EFFECT_SIZE_STRATIFIED_RESULTS.tsv",sep='\t',index=False); return out

def robust_metrics(d):
    rows=[]
    for name,fn in [('Spearman',spearman),('Pearson',pearson),('Kendall_tau',kendall),('Biweight_midcorrelation',bicor)]:
        r=fn(d.RNA_6h_log2FC,d.PXD054330_24h_protein_log2FC); rb=fn(d.Ribo_6h_log2FC,d.PXD054330_24h_protein_log2FC); rows.append({"metric":name,"RNA6_protein":r,"Ribo6_protein":rb,"Delta_rho_or_metric":rb-r,"direction":"Ribo6_gt_RNA6" if rb>r else "Ribo6_le_RNA6"})
    out=pd.DataFrame(rows); out.to_csv(OUT/"ROBUST_CORRELATION_SENSITIVITY.tsv",sep='\t',index=False); return out

def module_audit():
    x=pd.read_csv(OUT.parent/"c2/C2_MATCHED_NULL_SPECIFICITY.tsv",sep='\t'); x['module_type']=np.where(x.module=='mechanotransduction_network_composite','derived_composite','frozen_module'); x['BH_FDR']=bh(x.empirical_P.to_numpy()); x['classification_C2R']=np.where(x.BH_FDR<.05,'CONFIRMATORY_FDR<0.05',np.where(x.empirical_P<.05,'SUGGESTIVE_NOMINAL_ONLY','NOT_SUPPORTED')); x[['module','module_type','n_genes','observed_TRI_like','empirical_P','BH_FDR','specificity_class','classification_C2R']].to_csv(OUT/"MODULE_MULTIPLE_TESTING.tsv",sep='\t',index=False); return x

def matched_null_score(targets,score,master,seed):
    feats=master.set_index('ensembl_gene_id')[["baseline_RNA_control","mean_control","protein_missingness"]].dropna(); bins=pd.DataFrame(index=feats.index)
    for c in feats.columns: bins[c]=pd.qcut(feats[c],5,labels=False,duplicates='drop')
    targets=[t for t in targets if t in feats.index and np.isfinite(score.get(t,np.nan))]; idxarr=feats.index.to_numpy(); bmat=bins.to_numpy(float); sarr=score.reindex(idxarr).to_numpy(float); ti={int(np.where(idxarr==t)[0][0]) for t in targets}; pools=[]
    for t in targets:
        j=int(np.where(idxarr==t)[0][0]); mask=np.all(bmat==bmat[j],axis=1); mask[list(ti)]=False; pool=np.where(mask)[0]
        if len(pool)==0: pool=np.asarray([i for i in range(len(idxarr)) if i not in ti])
        pools.append(pool)
    rg=np.random.default_rng(seed); draws=np.column_stack([rg.choice(p,size=10000,replace=True) for p in pools]); return np.nanmean(sarr[draws],1),targets

def amino_robustness(m):
    defs=pd.read_csv(CM/"config/C2_MODULE_DEFINITIONS.tsv",sep='\t'); genes=set(defs.loc[defs.module=='amino_acid_osmolyte_transport','gene_symbol']); targets=sorted(set(m.loc[m.gene_symbol.isin(genes),'ensembl_gene_id'])); vals={c:(m[c]-m[c].mean())/m[c].std(ddof=1) for c in ['Ribo_1h_log2FC','Ribo_6h_log2FC','PXD054330_24h_protein_log2FC']}; score=(-vals['Ribo_1h_log2FC']+(vals['Ribo_6h_log2FC']-vals['Ribo_1h_log2FC'])+vals['PXD054330_24h_protein_log2FC'].abs())/3; score.index=m.ensembl_gene_id; null,targets=matched_null_score(targets,score,m,270933); rows=[]
    def add(label,ts):
        obs=float(np.nanmean([score[t] for t in ts])); p=float((1+np.sum(null>=obs))/(len(null)+1)); rows.append({"analysis":label,"n_genes":len(ts),"observed_TRI_like":obs,"null_median":float(np.nanmedian(null)),"empirical_P":p,"above_null_median":obs>np.nanmedian(null)})
    add('all_module_genes',targets)
    for t in targets: add('leave_one_out:'+str(m.loc[m.ensembl_gene_id==t,'gene_symbol'].iloc[0]),[g for g in targets if g!=t])
    top=sorted(targets,key=lambda t:abs(vals['PXD054330_24h_protein_log2FC'][m.ensembl_gene_id==t].iloc[0]),reverse=True)[:3]; add('leave_top3_abs_protein_effect',[g for g in targets if g not in top]); out=pd.DataFrame(rows); out.to_csv(OUT/"AMINO_ACID_PROGRAM_ROBUSTNESS.tsv",sep='\t',index=False); return out

def block_bootstrap(d):
    a=d[d.chromosome.notna()].copy(); chroms=sorted(a.chromosome.unique()); groups=[np.where(a.chromosome.to_numpy()==c)[0] for c in chroms]; x=a.RNA_6h_log2FC.to_numpy(float); rb=a.Ribo_6h_log2FC.to_numpy(float); y=a.PXD054330_24h_protein_log2FC.to_numpy(float); rx=stats.rankdata(x); rr=stats.rankdata(rb); ry=stats.rankdata(y); rg=np.random.default_rng(270932); vals=[]
    for _ in range(B):
        chosen=rg.integers(0,len(groups),size=len(groups)); idx=np.concatenate([groups[j] for j in chosen]);
        def corr(q):
            q=q[idx]; z=ry[idx]; q=q-q.mean(); z=z-z.mean(); return np.sum(q*z)/np.sqrt(np.sum(q*q)*np.sum(z*z))
        vals.append(corr(rr)-corr(rx))
    vals=np.asarray(vals); out=pd.DataFrame([{"n_genes":len(a),"n_chromosome_blocks":len(chroms),"Delta_rho_median":np.median(vals),"CI_low":np.quantile(vals,.025),"CI_high":np.quantile(vals,.975),"empirical_P":2*min(np.mean(vals<=0),np.mean(vals>=0)),"direction_stable":bool(np.quantile(vals,.025)>0)}]); out.to_csv(OUT/"BLOCK_BOOTSTRAP.tsv",sep='\t',index=False); return out

def external_validation():
    rows=[{"source":"PXD054330 publication","identifier":"10.1021/acs.jproteome.4c01046","evidence":"Qualitative support for SLC38A2/SNAT2 and GLS1/GLS increases under HCEC hyperosmotic stress and mitochondrial/osmotic adaptation.","use_in_primary_statistics":False},{"source":"Independent HCE hyperosmotic proteome abstract","identifier":"10.1111/aos.16032","evidence":"Qualitative support for reproducible hyperosmolar HCE proteome remodeling involving inflammation, ROS, apoptosis, cell viability and actin-cytoskeleton organization.","use_in_primary_statistics":False},{"source":"PXD059451","identifier":"PXD059451","evidence":"Used only as HCEC baseline detectability annotation; no hyperosmolarity outcome was imported.","use_in_primary_statistics":False}]
    out=pd.DataFrame(rows); out.to_csv(OUT/"EXTERNAL_QUALITATIVE_VALIDATION.tsv",sep='\t',index=False); return out

def reports(m,resid,partial,cvsum,cvdelta,perm_p,det,strat,robust,mod,amino,block,ext):
    r=resid.iloc[0]; pa=partial.iloc[0]; pb=partial.iloc[1]; a=cvsum[cvsum.model=='A'].iloc[0]; b=cvsum[cvsum.model=='B'].iloc[0]; c=cvsum[cvsum.model=='C'].iloc[0]; d=cvsum[cvsum.model=='D'].iloc[0]; q=mod.set_index('module'); low=strat[(strat.stratification=='baseline_RNA_quartile')&(strat.stratum=='Q1')].iloc[0]; high=strat[(strat.stratification=='baseline_RNA_quartile')&(strat.stratum=='Q4')].iloc[0]; large=strat[(strat.stratification=='protein_effect_abs_quartile')&(strat.stratum=='Q4')].iloc[0]; amino_pass=bool((amino.observed_TRI_like.iloc[1:] > amino.null_median.iloc[1:]).all() and (amino.empirical_P.iloc[1:]<.10).all()); extused=det.external_detectability_used.iloc[0]=='YES'; blockpass=bool(block.direction_stable.iloc[0]); residualpass=bool(r.PASS); cvpass=bool(cvdelta.Model_B_out_of_sample_better.iloc[0]); strong=bool(residualpass and cvpass and perm_p<.05 and blockpass and (strat[strat.stratification=='baseline_RNA_quartile'].Delta_rho>0).all()); moderate=bool(residualpass and cvpass and perm_p<.05); claim='A' if strong else 'B' if moderate else 'C'; verdict='MANUSCRIPT_GO_STRONG' if strong else 'MANUSCRIPT_GO_MODERATE' if moderate else 'MANUSCRIPT_GO_CROSSOMIC_ONLY' if r.residual_protein_spearman_rho>0 and r.CI_low>0 else 'HOLD'; title={'A':'Early translational remodeling anticipates later proteomic adaptation to hyperosmotic stress in human corneal epithelium','B':'Cross-omic temporal triangulation reveals translational signatures of corneal epithelial osmoadaptation','C':'Cross-omic mapping of corneal epithelial responses to hyperosmotic stress'}[claim]
    (REP/'C2R_INCREMENTAL_INFORMATION.md').write_text(f"# C2R incremental information\n\nRNA-adjusted Ribo residual→Protein24 ρ={r.residual_protein_spearman_rho:.5f}, 95% CI [{r.CI_low:.5f}, {r.CI_high:.5f}], empirical P={r.empirical_P:.5g}; residual test PASS={r.PASS}. Partial Ribo6→Protein24 | RNA6 ρ={pa.partial_rho:.5f}, CI [{pa.CI_low:.5f}, {pa.CI_high:.5f}], P={pa.empirical_P:.5g}. Partial RNA6→Protein24 | Ribo6 ρ={pb.partial_rho:.5f}, CI [{pb.CI_low:.5f}, {pb.CI_high:.5f}], P={pb.empirical_P:.5g}.\n\nModel A CV ρ={a.out_of_sample_Spearman_rho_mean:.5f}; Model B CV ρ={b.out_of_sample_Spearman_rho_mean:.5f}; B−A ΔCVρ={cvdelta.Delta_CV_rho.iloc[0]:.5f}, CI [{cvdelta.CI_low.iloc[0]:.5f}, {cvdelta.CI_high.iloc[0]:.5f}], P={cvdelta.empirical_P.iloc[0]:.5g}. Model C CV ρ={c.out_of_sample_Spearman_rho_mean:.5f}; Model D CV ρ={d.out_of_sample_Spearman_rho_mean:.5f}. Restricted permutation P={perm_p:.5g}.\n\nStratification note: the frozen master table lacks Ribo6 standard errors, so the Ribo precision strata use the explicitly labelled Ribo6 FDR precision proxy; no standard error or coverage value was fabricated.\n\nC2R verdict: **{verdict}**. Claim **{claim}** is the most conservative claim supported by the prespecified gates; the original C2 gate name is not used as a descriptor of the modest ρ values.\n")
    (REP/'C2R_DETECTABILITY_AUDIT.md').write_text("# C2R detectability audit\n\nPXD059451 was used only as a baseline HCEC detectability annotation. The supplementary consistency workbook was available locally; abundance ranks were unavailable and were not fabricated. Because PXD059451 uses a different protein-group list, mapping used its Protein Name human gene symbols against the frozen C2 Ensembl-gene table when direct accessions were unavailable.\n\n"+det.to_markdown(index=False)+"\n")
    (REP/'C2R_MODULE_MULTIPLICITY_AUDIT.md').write_text("# C2R module multiplicity audit\n\nBH correction was applied across all nine frozen modules and the derived mechanotransduction composite. Nominal P<0.05 with BH-FDR≥0.05 is labeled suggestive only.\n\n"+mod[['module','module_type','empirical_P','BH_FDR','classification_C2R']].to_markdown(index=False)+"\n")
    (REP/'C2R_EXTERNAL_VALIDATION.md').write_text("# C2R external qualitative validation\n\nNo external study was pooled into the primary statistics.\n\n"+ext.to_markdown(index=False)+"\n")
    (REP/'C2R_INTEGRATED_DECISION.md').write_text(f"# C2R integrated decision\n\n**C2R_VERDICT:** {verdict}\n\n- Residual incremental test: {'PASS' if residualpass else 'FAIL'}\n- Nested CV Model B>A: {'PASS' if cvpass else 'FAIL'}\n- Permutation P: {perm_p:.5g}\n- Block bootstrap direction stable: {'YES' if blockpass else 'NO'}\n- Detectability-restricted analysis used: {'YES' if extused else 'NO'}\n- Amino-acid program leave-one-out robustness: {'YES' if amino_pass else 'NO'}\n- Mechanotransduction confirmatory after BH: {'YES' if q.loc['mechanotransduction_network_composite','BH_FDR']<.05 else 'NO'}\n\n**Final claim option:** {claim}\n\n**Recommended title:** {title}\n\n**Remove from manuscript-primary claims:** EARLY_SUPPRESSION_LATE_RECOVERY, fate-map/recovery narrative, MeAIB-sensitive recovery program, Torin-sensitive recovery program, mechanotransduction as a core conclusion unless its BH-FDR is below 0.05.\n\n**Largest limitation:** all cross-study links remain modest and non-longitudinal; protein n=3, model/cell/osmolarity/platform differ, and the external detectability resource supplies consistency rather than abundance.\n\nFormal manuscript drafting remains stopped.\n")
    return verdict,claim,title

def main():
    m=read_inputs(); d,resid=residual_analysis(m); partial=partial_analysis(d); dc,folds,cvr,cvsum,cvdelta=nested_cv(m); perm_p,null=permutation_null(dc,folds,cvdelta); det,ann,high=detectability_analysis(m,d,folds); strat=stratification(m); robust=robust_metrics(d); mod=module_audit(); amino=amino_robustness(m); block=block_bootstrap(d); ext=external_validation(); verdict,claim,title=reports(m,resid,partial,cvsum,cvdelta,perm_p,det,strat,robust,mod,amino,block,ext)
    summary={"C2R_VERDICT":verdict,"claim_option":claim,"recommended_title":title,"residual":resid.iloc[0].to_dict(),"partial":partial.to_dict('records'),"cv_delta":cvdelta.iloc[0].to_dict(),"permutation_P":perm_p,"detectability":det.iloc[0].to_dict(),"block":block.iloc[0].to_dict(),"module_audit":mod[['module','empirical_P','BH_FDR','classification_C2R']].to_dict('records'),"outputs":sorted(str(p.relative_to(CM)) for p in OUT.glob('*') if p.is_file())}
    (OUT/'C2R_AUDIT_SUMMARY.json').write_text(json.dumps(summary,indent=2,default=lambda x:None)); print(json.dumps(summary,indent=2,default=lambda x:None))

if __name__=='__main__': main()
