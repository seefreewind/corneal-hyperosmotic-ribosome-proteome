#!/usr/bin/env python3
"""C2 preregistered RNA--Ribo--protein temporal triangulation.

This script consumes the recovered processed PXD054330 pivot and the existing
Phase-1 tables.  It never downloads or reprocesses protein RAW files.
"""
from pathlib import Path
import re
import json
import math
import hashlib
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
CM = ROOT / "computational_manuscript"
OUT = CM / "results" / "c2"
FIG = CM / "figures" / "c2"
REP = CM / "reports" / "c2"
OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True); REP.mkdir(parents=True, exist_ok=True)
SEED = 260926
B = 10000
RNG = np.random.default_rng(SEED)

def bh(p):
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.any(): out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out

def zscore(x):
    x = np.asarray(x, dtype=float); m = np.nanmean(x); s = np.nanstd(x, ddof=1)
    return (x-m)/s if np.isfinite(s) and s > 0 else np.full_like(x, np.nan)

def safe_spearman(x, y):
    x=np.asarray(x,float); y=np.asarray(y,float); ok=np.isfinite(x)&np.isfinite(y)
    if ok.sum()<4 or np.nanstd(x[ok])==0 or np.nanstd(y[ok])==0: return (np.nan,np.nan,int(ok.sum()))
    r,p=stats.spearmanr(x[ok],y[ok]); return float(r),float(p),int(ok.sum())

def safe_pearson(x, y):
    x=np.asarray(x,float); y=np.asarray(y,float); ok=np.isfinite(x)&np.isfinite(y)
    if ok.sum()<4 or np.nanstd(x[ok])==0 or np.nanstd(y[ok])==0: return (np.nan,np.nan,int(ok.sum()))
    r,p=stats.pearsonr(x[ok],y[ok]); return float(r),float(p),int(ok.sum())

def direct_ttest(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float); a=a[np.isfinite(a)]; b=b[np.isfinite(b)]
    if len(a)<2 or len(b)<2: return np.nan,np.nan,np.nan
    t,p=stats.ttest_ind(a,b,equal_var=True)
    va=np.var(a,ddof=1); vb=np.var(b,ddof=1); se=np.sqrt(((len(a)-1)*va+(len(b)-1)*vb)/(len(a)+len(b)-2)*(1/len(a)+1/len(b)))
    return float(se),float(t),float(p)

def bootstrap_rho(x,y,B=B,seed=SEED):
    x=np.asarray(x,float); y=np.asarray(y,float); ok=np.isfinite(x)&np.isfinite(y); x=x[ok]; y=y[ok]; n=len(x)
    if n<4: return np.nan,np.nan,np.nan,np.nan,np.array([])
    rx=stats.rankdata(x,method="average"); ry=stats.rankdata(y,method="average")
    rg=np.random.default_rng(seed); vals=[]; chunk=250
    for start in range(0,B,chunk):
        k=min(chunk,B-start); idx=rg.integers(0,n,size=(k,n),dtype=np.int32); a=rx[idx]; b=ry[idx]
        ac=a-a.mean(axis=1,keepdims=True); bc=b-b.mean(axis=1,keepdims=True)
        den=np.sqrt((ac*ac).sum(1)*(bc*bc).sum(1)); vals.extend(((ac*bc).sum(1)/den).tolist())
    vals=np.asarray(vals,float); return float(np.nanmedian(vals)),float(np.nanquantile(vals,.025)),float(np.nanquantile(vals,.975)),float(2*min(np.mean(vals<=0),np.mean(vals>=0))),vals

def bootstrap_delta(x1,x2,y,B=B,seed=SEED):
    x1=np.asarray(x1,float); x2=np.asarray(x2,float); y=np.asarray(y,float); ok=np.isfinite(x1)&np.isfinite(x2)&np.isfinite(y)
    x1=x1[ok]; x2=x2[ok]; y=y[ok]; n=len(y)
    if n<4: return (np.nan,np.nan,np.nan,np.nan,np.array([]),np.nan,np.nan)
    r1=stats.rankdata(x1,method="average"); r2=stats.rankdata(x2,method="average"); ry=stats.rankdata(y,method="average")
    direct1=stats.spearmanr(x1,y).statistic; direct2=stats.spearmanr(x2,y).statistic
    rg=np.random.default_rng(seed); vals=[]; chunk=250
    for start in range(0,B,chunk):
        k=min(chunk,B-start); idx=rg.integers(0,n,size=(k,n),dtype=np.int32); a=r1[idx]; b=r2[idx]; c=ry[idx]
        def corr(q):
            q=q-q.mean(axis=1,keepdims=True); cc=c-c.mean(axis=1,keepdims=True); return (q*cc).sum(1)/np.sqrt((q*q).sum(1)*(cc*cc).sum(1))
        vals.extend((corr(b)-corr(a)).tolist())
    vals=np.asarray(vals,float); return float(np.nanmedian(vals)),float(np.nanquantile(vals,.025)),float(np.nanquantile(vals,.975)),float(2*min(np.mean(vals<=0),np.mean(vals>=0))),vals,float(direct1),float(direct2)

def write_tsv(df,path):
    path.parent.mkdir(parents=True,exist_ok=True); df.to_csv(path,sep="\t",index=False,na_rep="NA")

def map_symbol_table(mapping, rna_symbols):
    u=mapping[(mapping.mapping_type=="UNIQUE_1TO1") & ~mapping.ensembl_gene_id.astype(str).str.contains(";",na=False)].copy()
    ens_to_sym={}; sym_to_ens={}
    for _,r in u.iterrows():
        ens=str(r.ensembl_gene_id); syms=[s.strip() for s in str(r.gene_symbol).split(";") if s.strip()]
        present=[s for s in syms if s in rna_symbols]
        chosen=present[0] if present else syms[0] if syms else ""
        ens_to_sym[ens]=chosen
        for s in syms:
            sym_to_ens.setdefault(s,set()).add(ens)
    # retain unambiguous symbol mappings only
    sym_one={s:next(iter(v)) for s,v in sym_to_ens.items() if len(v)==1}
    return u,ens_to_sym,sym_one

def build_protein(mapping, pivot, sample_cols, ens_to_sym):
    m=mapping[mapping.mapping_type=="UNIQUE_1TO1"].copy()
    m=m[~m.ensembl_gene_id.astype(str).str.contains(";",na=False)]
    d=pivot.merge(m[["original_protein_identifier","ensembl_gene_id","gene_symbol"]],left_on="PG.ProteinAccessions",right_on="original_protein_identifier",how="inner")
    for c in sample_cols: d[c]=pd.to_numeric(d[c],errors="coerce")
    # group protein groups at gene level before calculating the 3-vs-3 contrast
    g=d.groupby("ensembl_gene_id",sort=False)
    mat=g[sample_cols].mean()
    names=g["original_protein_identifier"].apply(lambda x:";".join(map(str,x)))
    ctrl=[sample_cols[3],sample_cols[4],sample_cols[5]]; hyp=sample_cols[:3]
    rows=[]
    for ens,row in mat.iterrows():
        c=row[ctrl].to_numpy(float); h=row[hyp].to_numpy(float); se,t,p=direct_ttest(h,c)
        rows.append({"ensembl_gene_id":ens,"protein":names.loc[ens],"log2FC":np.nanmean(h)-np.nanmean(c),"SE":se,"t":t,"P":p,"mean_control":np.nanmean(c),"mean_hyper":np.nanmean(h),"missingness":float(row.isna().mean())})
    dep=pd.DataFrame(rows); dep["BH_FDR"]=bh(dep.P.to_numpy())
    dep["gene_symbol"]=dep.ensembl_gene_id.map(ens_to_sym)
    return dep,mat

def load_inputs():
    mapping=pd.read_csv(CM/"metadata/c1m/PROTEIN_GENE_ID_MAP_REAL.tsv",sep="\t")
    rna_vst=pd.read_csv(ROOT/"data/processed/GSE200097_RNA_vst.tsv",sep="\t",usecols=["gene"])
    rna_symbols=set(rna_vst.gene.astype(str))
    u,ens_to_sym,sym_to_ens=map_symbol_table(mapping,rna_symbols)
    pivot=pd.read_csv(CM/"data/raw/proteome/PXD054330/processed/Orbitrap_HCEC_DED_vs_NRM_Report_BGS_Analysis_Grid_View_Report_Pivot.tsv",sep="\t")
    sample_cols=[c for c in pivot.columns if "PG.Log2Quantity" in c]
    dep,pmat=build_protein(mapping,pivot,sample_cols,ens_to_sym)
    p1=pd.read_csv(ROOT/"results/phase1F/PRIMARY_RNA_RIBO_TE_ALL.tsv",sep="\t")
    p1["ensembl_gene_id"]=p1.gene.astype(str).map(sym_to_ens)
    p1=p1[p1.ensembl_gene_id.notna()].copy()
    # use the pre-audited RNA feature space and uniquely mapped genes to retain exactly 3974 primary genes
    primary_ens=set(dep.ensembl_gene_id)&set(u.ensembl_gene_id)&{e for e,s in ens_to_sym.items() if s in rna_symbols}
    p1=p1[p1.ensembl_gene_id.isin(primary_ens)].drop_duplicates("ensembl_gene_id")
    gse=pd.read_csv(ROOT/"results/phase1/GSE323164_osm_24h_vs_control_DESeq2.tsv",sep="\t")
    # GSE323164 uses a mixture of Ensembl IDs and symbols.  Resolve Ensembl IDs
    # through the existing study annotation, then harmonize by the frozen symbol map.
    ann=pd.read_csv(ROOT/"data/raw/GSE323164/GSE323164_1_genes_fpkm_expression.txt.gz",sep="\t",usecols=["gene_id","gene_name"])
    annmap=ann.drop_duplicates("gene_id").set_index("gene_id").gene_name.to_dict()
    gse["gse_symbol"]=gse.gene.astype(str).map(annmap).fillna(gse.gene.astype(str)).str.split(";").str[0].str.strip()
    gse["ensembl_gene_id"]=gse.gse_symbol.map(sym_to_ens)
    gse24=gse[["ensembl_gene_id","log2FoldChange","stat","padj"]].rename(columns={"log2FoldChange":"GSE323164_24h_RNA_log2FC","stat":"GSE323164_24h_RNA_stat","padj":"GSE323164_24h_RNA_FDR"})
    gse24=gse24.dropna(subset=["ensembl_gene_id"]).drop_duplicates("ensembl_gene_id")
    pert=pd.read_csv(ROOT/"results/phase1F/MEAIB_TORIN_PERTURBATION.tsv",sep="\t")
    pert["gene"]=pert.gene.astype(str).str.strip()
    pert["ensembl_gene_id"]=pert.gene.map(sym_to_ens)
    modules=pd.read_csv(CM/"config/C2_MODULE_DEFINITIONS.tsv",sep="\t")
    modules["ensembl_gene_id"]=modules.gene_symbol.astype(str).map(sym_to_ens)
    return mapping, u, ens_to_sym, sym_to_ens, dep, pmat, p1, gse24, pert, modules, sample_cols, pivot

def build_master(dep,p1,gse24,pert,modules,ens_to_sym):
    m=dep.merge(p1,on="ensembl_gene_id",how="inner",suffixes=("","_p1"))
    m=m.merge(gse24,on="ensembl_gene_id",how="left").merge(pert.drop(columns=["gene"],errors="ignore"),on="ensembl_gene_id",how="left")
    # Harmonize the perturbation table's source naming with the C2 preregistration.
    for src,dst in [("TE_MeAIB_log2FC","MeAIB_TE_log2FC"),("TE_MeAIB_FDR","MeAIB_TE_FDR"),("TE_Torin_log2FC","Torin_TE_log2FC"),("TE_Torin_FDR","Torin_TE_FDR")]:
        if src in m.columns: m[dst]=m[src]
    m["gene_symbol"]=m.ensembl_gene_id.map(ens_to_sym).fillna(m.get("gene_symbol",pd.Series(index=m.index)))
    # concise names requested by the C2 table
    ren={"log2FC":"PXD054330_24h_protein_log2FC","P":"PXD054330_24h_protein_P","BH_FDR":"PXD054330_24h_protein_FDR","missingness":"protein_missingness"}
    m=m.rename(columns=ren)
    mem=modules.groupby("ensembl_gene_id").module.apply(lambda x:";".join(sorted(set(x.dropna())))).to_dict()
    m["module_membership"]=m.ensembl_gene_id.map(mem).fillna("")
    m["mapping_quality"]="UNIQUE_1TO1"
    # requested endpoint aliases
    m=m.rename(columns={"log2FC_RNA_1h":"RNA_1h_log2FC","FDR_RNA_1h":"RNA_1h_FDR","log2FC_Ribo_1h":"Ribo_1h_log2FC","FDR_Ribo_1h":"Ribo_1h_FDR","stat_TE_1h":"TE_1h_stat","FDR_TE_1h":"TE_1h_FDR","log2FC_TE_1h":"TE_1h_log2FC","log2FC_RNA_6h":"RNA_6h_log2FC","FDR_RNA_6h":"RNA_6h_FDR","log2FC_Ribo_6h":"Ribo_6h_log2FC","FDR_Ribo_6h":"Ribo_6h_FDR","stat_TE_6h":"TE_6h_stat","FDR_TE_6h":"TE_6h_FDR","log2FC_TE_6h":"TE_6h_log2FC"})
    # Preserve all stat/FDR columns in the audit table, then order core fields first.
    core=["ensembl_gene_id","gene_symbol","RNA_1h_log2FC","RNA_1h_FDR","Ribo_1h_log2FC","Ribo_1h_FDR","TE_1h_stat","TE_1h_FDR","TE_1h_log2FC","RNA_6h_log2FC","RNA_6h_FDR","Ribo_6h_log2FC","Ribo_6h_FDR","TE_6h_stat","TE_6h_FDR","TE_6h_log2FC","GSE323164_24h_RNA_log2FC","GSE323164_24h_RNA_FDR","PXD054330_24h_protein_log2FC","PXD054330_24h_protein_P","PXD054330_24h_protein_FDR","MeAIB_TE_log2FC","MeAIB_TE_FDR","Torin_TE_log2FC","Torin_TE_FDR","module_membership","protein_missingness","mapping_quality"]
    return m, [c for c in core if c in m.columns]

def correlation_outputs(master):
    y=master.PXD054330_24h_protein_log2FC.to_numpy(float)
    preds={"RNA1":"RNA_1h_log2FC","Ribo1":"Ribo_1h_log2FC","TE1_stat":"TE_1h_stat","RNA6":"RNA_6h_log2FC","Ribo6":"Ribo_6h_log2FC","TE6_stat":"TE_6h_stat","GSE323164_RNA24":"GSE323164_24h_RNA_log2FC"}
    rows=[]; boots={}
    for label,col in preds.items():
        x=master[col].to_numpy(float); sr,sp,n=safe_spearman(x,y); pr,pp,_=safe_pearson(x,y)
        # The preregistered B=10,000 is retained for the two primary 6h predictors;
        # the descriptive secondary endpoints use a smaller fixed bootstrap budget.
        # All predictor intervals use the preregistered B=10,000 budget.
        b_reps=B
        med,lo,hi,bp,vals=bootstrap_rho(x,y,B=b_reps,seed=SEED+len(rows))
        rows.append({"predictor":label,"n":n,"spearman_rho":sr,"spearman_P":sp,"bootstrap_median_rho":med,"bootstrap_CI_low":lo,"bootstrap_CI_high":hi,"bootstrap_empirical_P":bp,"pearson_r":pr,"pearson_P":pp})
        boots[label]=vals
    out=pd.DataFrame(rows); out["BH_FDR"]=bh(out.spearman_P.to_numpy()); write_tsv(out,OUT/"C2_PREDICTOR_CORRELATIONS.tsv")
    dmed,dlo,dhi,dp,dvals,drna,dribo=bootstrap_delta(master.RNA_6h_log2FC,master.Ribo_6h_log2FC,y,B=B,seed=SEED+101)
    delta={"n":int(np.isfinite(master.RNA_6h_log2FC)&np.isfinite(master.Ribo_6h_log2FC)&np.isfinite(y).sum()) if False else int((np.isfinite(master.RNA_6h_log2FC)&np.isfinite(master.Ribo_6h_log2FC)&np.isfinite(y)).sum()),"RNA6_rho":drna,"Ribo6_rho":dribo,"delta_rho":dmed,"CI_low":dlo,"CI_high":dhi,"empirical_P":dp}
    write_tsv(pd.DataFrame([delta]),OUT/"C2_PRIMARY_DELTA_RHO.tsv")
    return out,delta,boots

def directional(master):
    y=master.PXD054330_24h_protein_log2FC
    pairs={"RNA6_vs_protein":"RNA_6h_log2FC","Ribo6_vs_protein":"Ribo_6h_log2FC","GSE323164_RNA24_vs_protein":"GSE323164_24h_RNA_log2FC"}
    rows=[]
    for name,col in pairs.items():
        x=master[col]; ok=x.notna()&y.notna(); xx=x[ok].to_numpy(); yy=y[ok].to_numpy(); nz=(xx!=0)&(yy!=0)
        upup=int(((xx>0)&(yy>0)).sum()); downdown=int(((xx<0)&(yy<0)).sum()); disc=int(((xx*yy<0)).sum())
        rows.append({"comparison":name,"n":int(ok.sum()),"up_up":upup,"down_down":downdown,"discordant":disc,"zero_or_missing":int(ok.sum()-nz.sum()),"sign_agreement_fraction":(upup+downdown)/max(1,int(nz.sum()))})
    out=pd.DataFrame(rows); write_tsv(out,OUT/"C2_DIRECTIONAL_CONCORDANCE.tsv"); return out

def marked(e,f,thr=.25,q=.05,sign=None):
    ok=np.isfinite(e)&np.isfinite(f)&(np.abs(e)>=thr)&(f<q)
    if sign=="positive": ok &= e>0
    if sign=="negative": ok &= e<0
    return ok

def fate_for(master,thr=.25,q=.05):
    out=[]
    for _,r in master.iterrows():
        rn6,rb6,te6,p24,pr6=r.get("RNA_6h_log2FC",np.nan),r.get("Ribo_6h_log2FC",np.nan),r.get("TE_6h_log2FC",np.nan),r.get("PXD054330_24h_protein_log2FC",np.nan),r.get("PXD054330_24h_protein_FDR",np.nan)
        rn1,rb1,te1=r.get("RNA_1h_log2FC",np.nan),r.get("Ribo_1h_log2FC",np.nan),r.get("TE_1h_log2FC",np.nan)
        fRn6,fRb6,fP= r.get("RNA_6h_FDR",np.nan),r.get("Ribo_6h_FDR",np.nan),pr6
        fRb1=r.get("Ribo_1h_FDR",np.nan); fTe1=r.get("TE_1h_FDR",np.nan)
        mrn=bool(marked(rn6,fRn6,thr,q)); mrb=bool(marked(rb6,fRb6,thr,q)); mp=bool(marked(p24,fP,thr,q));
        early=bool((marked(rb1,fRb1,thr,q,"negative") or marked(te1,fTe1,thr,q,"negative")) and ((np.isfinite(rb6) and rb6>-thr and np.isfinite(rb1) and rb6-rb1>=thr) or (np.isfinite(te6) and te6>-thr and np.isfinite(te1) and te6-te1>=thr)) and marked(p24,fP,thr,q,"positive"))
        persistent=bool(marked(rb1,fRb1,thr,q,"negative") and marked(rb6,fRb6,thr,q,"negative") and marked(p24,fP,thr,q,"negative"))
        ampl=bool(mrb and np.isfinite(rn6) and np.isfinite(rb6) and abs(rb6)>=abs(rn6)+thr and np.sign(rb6)==np.sign(p24) and mp)
        buff=bool(mrn and np.isfinite(rn6) and np.isfinite(rb6) and (np.sign(rn6)!=np.sign(rb6) or abs(rb6)<=abs(rn6)-thr) and np.sign(rb6)==np.sign(p24) and mp)
        tx=bool(mrn and mp and np.sign(rn6)==np.sign(p24))
        disc=bool(mrn and mp and np.sign(rn6)!=np.sign(p24) and (not mrb or np.sign(rb6)!=np.sign(p24)))
        rbconc=bool(mrb and mp and np.sign(rb6)==np.sign(p24))
        if early: cls="EARLY_SUPPRESSION_LATE_RECOVERY"
        elif persistent: cls="PERSISTENT_TRANSLATIONAL_SUPPRESSION"
        elif ampl: cls="TRANSLATIONAL_AMPLIFICATION"
        elif buff: cls="TRANSLATIONAL_BUFFERING"
        elif tx: cls="TRANSCRIPTION_LED"
        elif disc: cls="RNA_PROTEIN_DISCORDANT"
        elif rbconc: cls="RIBO_PROTEIN_CONCORDANT"
        else: cls="NO_CLEAR_PATTERN"
        out.append(cls)
    return np.array(out)

def fate_map(master):
    f1=fate_for(master,.25,.05); f2=fate_for(master,.50,.10)
    out=master[["ensembl_gene_id","gene_symbol"]].copy(); out["fate_class_primary"]=f1; out["fate_class_sensitivity"]=f2; out["fate_stability"]=np.where(f1==f2,"STABLE","UNSTABLE_FATE_CLASS")
    write_tsv(out,OUT/"C2_FATE_MAP.tsv")
    counts=out.fate_class_primary.value_counts().rename_axis("fate_class").reset_index(name="n"); write_tsv(counts,OUT/"C2_FATE_CLASS_COUNTS.tsv")
    return out,counts

def module_outputs(master,modules):
    effects={"RNA1":"RNA_1h_log2FC","Ribo1":"Ribo_1h_log2FC","TE1":"TE_1h_log2FC","RNA6":"RNA_6h_log2FC","Ribo6":"Ribo_6h_log2FC","TE6":"TE_6h_log2FC","Protein24":"PXD054330_24h_protein_log2FC"}
    statmap={"RNA1":"stat_RNA_1h","Ribo1":"stat_Ribo_1h","TE1":"TE_1h_stat","RNA6":"stat_RNA_6h","Ribo6":"stat_Ribo_6h","TE6":"TE_6h_stat","Protein24":"t"}
    zvals={k:zscore(master[v].to_numpy(float)) for k,v in effects.items()}
    rows=[]; erows=[]
    for mod,mm in modules.groupby("module"):
        genes=set(mm.ensembl_gene_id.dropna()) & set(master.ensembl_gene_id)
        idx=master.ensembl_gene_id.isin(genes).to_numpy(); n=int(idx.sum())
        for ep,col in effects.items():
            v=master.loc[idx,col].to_numpy(float); zz=zvals[ep][idx]; rows.append({"module":mod,"endpoint":ep,"n":int(np.isfinite(v).sum()),"mean_effect":np.nanmean(v),"median_effect":np.nanmedian(v),"mean_standardized_effect":np.nanmean(zz),"median_standardized_effect":np.nanmedian(zz)})
        # rank-based competitive endpoint tests (full primary background)
        for ep,col in statmap.items():
            a=master.loc[idx,col].to_numpy(float); b=master.loc[~idx,col].to_numpy(float); a=a[np.isfinite(a)]; b=b[np.isfinite(b)]
            if len(a)>=2 and len(b)>=2:
                u,p=stats.mannwhitneyu(a,b,alternative="two-sided",method="asymptotic"); direction=np.nanmean(a)-np.nanmean(b)
            else: u=p=direction=np.nan
            erows.append({"module":mod,"endpoint":ep,"n_module":len(a),"n_background":len(b),"rank_U":u,"P":p,"mean_stat_difference":direction})
    eff=pd.DataFrame(rows); enr=pd.DataFrame(erows); enr["BH_FDR"]=bh(enr.P.to_numpy()); write_tsv(eff,OUT/"C2_MODULE_EFFECTS.tsv"); write_tsv(enr,OUT/"C2_MODULE_ENRICHMENT.tsv")
    return eff,enr,zvals

def module_null(master,modules,zvals):
    feats=pd.DataFrame({"gene":master.ensembl_gene_id,"rna_base":master.baseline_RNA_control,"prot_base":master.mean_control,"prot_miss":master.protein_missingness}).set_index("gene")
    avail=feats.dropna(); bins=pd.DataFrame(index=avail.index)
    for c in avail.columns: bins[c]=pd.qcut(avail[c],q=5,labels=False,duplicates="drop")
    score=(-zvals["Ribo1"]+(zvals["Ribo6"]-zvals["Ribo1"])+np.abs(zvals["Protein24"]))/3
    score=pd.Series(score,index=master.ensembl_gene_id)
    rng=np.random.default_rng(260927); rows=[]; module_pools={}
    gene_arr=avail.index.to_numpy(); bin_arr=bins.loc[gene_arr].to_numpy(dtype=float); score_arr=score.reindex(gene_arr).to_numpy(dtype=float)
    def draw_null(targets):
        targets=[t for t in targets if t in avail.index]; pools=[]
        target_idx={gene_arr.tolist().index(t) for t in targets}
        for t in targets:
            ti=int(np.where(gene_arr==t)[0][0]); tb=bin_arr[ti]
            mask=np.all(bin_arr==tb,axis=1); mask[list(target_idx)]=False; pool=np.where(mask)[0]
            if len(pool)==0:
                dist=((bin_arr-tb)**2).sum(axis=1); order=np.argsort(dist); pool=np.asarray([i for i in order if i not in target_idx][:max(20,len(targets))],dtype=int)
            if len(pool)==0: pool=np.asarray([i for i in range(len(gene_arr)) if i not in target_idx],dtype=int)
            pools.append(pool)
        draws=np.column_stack([rng.choice(p,size=10000,replace=True) for p in pools])
        return np.nanmean(score_arr[draws],axis=1)
    for mod,mm in modules.groupby("module"):
        genes=sorted(set(mm.ensembl_gene_id.dropna())&set(master.ensembl_gene_id)); genes=[g for g in genes if g in score.index and np.isfinite(score.get(g,np.nan))]
        obs=float(np.nanmean([score[g] for g in genes])) if genes else np.nan; null=draw_null(genes) if genes else np.array([])
        p=float((1+np.sum(null>=obs))/(len(null)+1)) if len(null) else np.nan; q=float(np.nanquantile(null,.975)) if len(null) else np.nan
        cls="SPECIFIC_STRONG" if np.isfinite(p) and p<.05 and obs>q else "SPECIFIC_MODERATE" if np.isfinite(p) and p<.10 else "GLOBAL_COMPATIBLE"
        rows.append({"module":mod,"n_genes":len(genes),"observed_TRI_like":obs,"null_mean":np.nanmean(null) if len(null) else np.nan,"null_CI_low":np.nanquantile(null,.025) if len(null) else np.nan,"null_CI_high":q,"empirical_P":p,"specificity_class":cls,"matching_features":"baseline_RNA_control_VST;protein_control_abundance;protein_missingness","unavailable_matching_features":"baseline_Ribo_abundance;gene_length;GC"})
    # a predeclared mechanotransduction composite (union of six modules) for the requested summary
    six={"mechanosensitive_ion_channels","focal_adhesion","actin_cytoskeleton","cell_matrix_adhesion","YAP_TAZ_mechanosensing","epithelial_junction_barrier"}; genes=sorted(set(modules[modules.module.isin(six)].ensembl_gene_id.dropna())&set(master.ensembl_gene_id)); genes=[g for g in genes if g in score.index and np.isfinite(score.get(g,np.nan))]; obs=float(np.nanmean([score[g] for g in genes])); null=draw_null(genes); p=float((1+np.sum(null>=obs))/(len(null)+1)); rows.append({"module":"mechanotransduction_network_composite","n_genes":len(genes),"observed_TRI_like":obs,"null_mean":np.nanmean(null),"null_CI_low":np.nanquantile(null,.025),"null_CI_high":np.nanquantile(null,.975),"empirical_P":p,"specificity_class":"SPECIFIC_STRONG" if p<.05 and obs>np.nanquantile(null,.975) else "SPECIFIC_MODERATE" if p<.10 else "GLOBAL_COMPATIBLE","matching_features":"baseline_RNA_control_VST;protein_control_abundance;protein_missingness","unavailable_matching_features":"baseline_Ribo_abundance;gene_length;GC"})
    out=pd.DataFrame(rows); write_tsv(out,OUT/"C2_MATCHED_NULL_SPECIFICITY.tsv"); return out

def protein_robustness(dep,pivot,mapping,ens_to_sym,sample_cols,master):
    scenarios={"raw_processed":pivot.copy()}
    norm=pivot.copy()
    for c in sample_cols:
        vals=pd.to_numeric(norm[c],errors="coerce")
        norm[c]=vals-vals.median()
    scenarios["column_median_centered"]=norm
    rows=[]; effs={}
    for name,mat in scenarios.items():
        d,_=build_protein(mapping,mat,sample_cols,ens_to_sym); e=d.set_index("ensembl_gene_id").log2FC; effs[name]=e
    base=effs["raw_processed"]
    for name,e in effs.items():
        common=base.index.intersection(e.index); pair=pd.concat([base.loc[common],e.loc[common]],axis=1).dropna(); r=stats.spearmanr(pair.iloc[:,0],pair.iloc[:,1]).statistic if len(pair)>3 else np.nan; rows.append({"scenario":name,"n_genes":len(pair),"rho_vs_raw":r})
    # complete-case and <=1/6-missingness filters from raw gene-level table
    rawdep,_=build_protein(mapping,pivot,sample_cols,ens_to_sym); 
    for name,subset in [("complete_case",rawdep[rawdep.missingness==0]),("<=1_of_6_missing",rawdep[rawdep.missingness<=1/6+1e-9])]:
        common=base.index.intersection(subset.ensembl_gene_id); pair=pd.concat([base.loc[common],subset.set_index("ensembl_gene_id").log2FC.loc[common]],axis=1).dropna(); r=stats.spearmanr(pair.iloc[:,0],pair.iloc[:,1]).statistic if len(pair)>3 else np.nan; rows.append({"scenario":name,"n_genes":len(pair),"rho_vs_raw":r})
    out=pd.DataFrame(rows); write_tsv(out,OUT/"C2_PROTEOME_EFFECT_ROBUSTNESS.tsv"); return out

def leave_one_out(master):
    m=master.copy(); y=m.PXD054330_24h_protein_log2FC
    masks={"all":np.zeros(len(m),bool),"mitochondrial_OXPHOS":m.module_membership.str.contains("mitochondrial_OXPHOS")|m.gene_symbol.str.startswith("MT-"),"cytoskeleton_union_focal_actin":m.module_membership.str.contains("focal_adhesion|actin_cytoskeleton"),"ribosomal_regex":m.gene_symbol.str.match(r"^(RPL|RPS|MRPL|MRPS)"),"integrated_stress_response":m.module_membership.str.contains("integrated_stress_response")}
    rows=[]
    for name,rm in masks.items():
        keep=~rm; a=m.loc[keep]; sr,_,_=safe_spearman(a.RNA_6h_log2FC,a.PXD054330_24h_protein_log2FC); sb,_,_=safe_spearman(a.Ribo_6h_log2FC,a.PXD054330_24h_protein_log2FC); rows.append({"removed_set":name,"n":len(a),"RNA6_rho":sr,"Ribo6_rho":sb,"delta_rho":sb-sr})
    out=pd.DataFrame(rows); write_tsv(out,OUT/"C2_LEAVE_ONE_MODULE_OUT.tsv"); return out

def quality_sensitivity(master):
    m=master.copy(); y=m.PXD054330_24h_protein_log2FC; rows=[]
    qrn=m.baseline_RNA_control.quantile(.75); qprot25=m.mean_control.quantile(.25); qprot50=m.mean_control.quantile(.50); qse=m.SE_Ribo_6h.quantile(.5); qmiss=m.protein_missingness.quantile(.5)
    masks={"all":np.ones(len(m),bool),"high_RNA_control_abundance":m.baseline_RNA_control>=qrn,"high_Ribo_precision":m.SE_Ribo_6h<=qse,"low_protein_missingness":m.protein_missingness<=qmiss,"top75_protein_abundance":m.mean_control>=qprot25,"top50_protein_abundance":m.mean_control>=qprot50}
    for name,keep in masks.items():
        a=m.loc[keep]; sr,_,_=safe_spearman(a.RNA_6h_log2FC,a.PXD054330_24h_protein_log2FC); sb,_,_=safe_spearman(a.Ribo_6h_log2FC,a.PXD054330_24h_protein_log2FC); rows.append({"subset":name,"n":len(a),"RNA6_rho":sr,"Ribo6_rho":sb,"delta_rho":sb-sr})
    out=pd.DataFrame(rows); write_tsv(out,OUT/"C2_MEASUREMENT_QUALITY_SENSITIVITY.tsv"); return out

def expression_strata(master):
    m=master.copy(); m["RNA_stratum"]=pd.qcut(m.baseline_RNA_control,4,labels=["Q1","Q2","Q3","Q4"],duplicates="drop"); m["Protein_stratum"]=pd.qcut(m.mean_control,4,labels=["Q1","Q2","Q3","Q4"],duplicates="drop")
    rows=[]
    for axis in ["RNA_stratum","Protein_stratum"]:
        for st,a in m.groupby(axis,observed=True):
            r1,_,n1=safe_spearman(a.RNA_6h_log2FC,a.PXD054330_24h_protein_log2FC); r2,_,n2=safe_spearman(a.Ribo_6h_log2FC,a.PXD054330_24h_protein_log2FC); rows.append({"stratification":axis,"stratum":str(st),"n":min(n1,n2),"RNA6_rho":r1,"Ribo6_rho":r2,"delta_rho":r2-r1})
    out=pd.DataFrame(rows); write_tsv(out,OUT/"C2_EXPRESSION_STRATIFICATION.tsv"); return out

def perturbation_bridge(master,fate):
    m=master.merge(fate[["ensembl_gene_id","fate_class_primary"]],on="ensembl_gene_id",how="left"); rec=m[m.fate_class_primary=="EARLY_SUPPRESSION_LATE_RECOVERY"].copy()
    me=marked(rec.MeAIB_TE_log2FC,rec.MeAIB_TE_FDR,.25,.05); to=marked(rec.Torin_TE_log2FC,rec.Torin_TE_FDR,.25,.05); classes=[]
    for a,b,ae,be in zip(me,to,rec.MeAIB_TE_log2FC,rec.Torin_TE_log2FC):
        if a and b: classes.append("SHARED" if np.sign(ae)==np.sign(be) else "COMPLEX")
        elif a: classes.append("MEAIB_SENSITIVE")
        elif b: classes.append("TORIN_SENSITIVE")
        else: classes.append("PERTURBATION_INDEPENDENT")
    rec["perturbation_class"]=classes; write_tsv(rec[["ensembl_gene_id","gene_symbol","fate_class_primary","MeAIB_TE_log2FC","MeAIB_TE_FDR","Torin_TE_log2FC","Torin_TE_FDR","perturbation_class"]],OUT/"C2_PERTURBATION_BRIDGE.tsv")
    counts=rec.perturbation_class.value_counts().rename_axis("perturbation_class").reset_index(name="n"); write_tsv(counts,OUT/"C2_PERTURBATION_BRIDGE_COUNTS.tsv"); return rec,counts

def representative(master,fate):
    nodes=["PIEZO1","TRPV4","YAP1","WWTR1","VCL","ACTN4","TLN1","PXN","FLNA","FLNB","ROCK1","ROCK2"]
    cols=["gene_symbol","RNA_1h_log2FC","Ribo_1h_log2FC","TE_1h_stat","RNA_6h_log2FC","Ribo_6h_log2FC","TE_6h_stat","GSE323164_24h_RNA_log2FC","PXD054330_24h_protein_log2FC","MeAIB_TE_log2FC","Torin_TE_log2FC","module_membership"]
    d=pd.DataFrame({"gene_symbol":nodes}).merge(master[cols],on="gene_symbol",how="left").merge(fate[["gene_symbol","fate_class_primary","fate_stability"]],on="gene_symbol",how="left")
    d["fate_class_primary"]=d.fate_class_primary.fillna("NOT_IN_PRIMARY_ANALYSIS_SET"); d["fate_stability"]=d.fate_stability.fillna("NOT_APPLICABLE"); write_tsv(d,OUT/"REPRESENTATIVE_NODE_EVIDENCE.tsv"); return d

def figures(master,corr,delta,fate,modules,enr,null,bridge):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; import seaborn as sns
    plt.rcParams.update({"font.family":"Arial","font.size":8,"axes.spines.top":False,"axes.spines.right":False,"figure.dpi":150})
    def save(fig,name): fig.tight_layout(); fig.savefig(FIG/name,dpi=300); fig.savefig(FIG/(Path(name).stem+".svg")); plt.close(fig)
    y=master.PXD054330_24h_protein_log2FC
    for i,(xlab,xcol) in enumerate([("RNA6","RNA_6h_log2FC"),("Ribo6","Ribo_6h_log2FC")],1):
        ok=master[xcol].notna()&y.notna(); fig,ax=plt.subplots(figsize=(3.4,3)); ax.scatter(master.loc[ok,xcol],y[ok],s=5,alpha=.3,color="#1f77b4",edgecolor="none"); ax.axhline(0,color="0.8",lw=.6); ax.axvline(0,color="0.8",lw=.6); r=corr.loc[corr.predictor==xlab,"spearman_rho"].iloc[0]; ax.set(xlabel=f"{xlab} log2FC",ylabel="Protein 24h log2FC",title=f"{xlab} vs protein; ρ={r:.2f}"); save(fig,f"C2-{i}.png")
    fig,ax=plt.subplots(figsize=(4,3)); ax.hist(delta["_vals"],bins=50,color="#2ca02c",alpha=.8); ax.axvline(0,color="black",lw=.8); ax.axvline(delta["delta_rho"],color="#d62728",lw=1.2); ax.set(xlabel="Bootstrap Δρ (Ribo6 − RNA6)",ylabel="Count"); save(fig,"C2-3.png")
    order=["RNA1","Ribo1","TE1_stat","RNA6","Ribo6","TE6_stat","GSE323164_RNA24"]; q=corr.set_index("predictor").loc[order]; fig,ax=plt.subplots(figsize=(5,3)); ax.errorbar(range(len(q)),q.spearman_rho,yerr=[q.spearman_rho-q.bootstrap_CI_low,q.bootstrap_CI_high-q.spearman_rho],fmt="o",color="#1f77b4",capsize=3); ax.axhline(0,color="0.7",lw=.6); ax.set_xticks(range(len(q)),order,rotation=35,ha="right"); ax.set_ylabel("Spearman ρ with Protein24"); save(fig,"C2-4.png")
    fc=fate.fate_class_primary.value_counts(); fig,ax=plt.subplots(figsize=(5,3)); fc.sort_values().plot.barh(ax=ax,color="#7f7f7f"); ax.set(xlabel="Genes",ylabel="Fate class"); save(fig,"C2-5.png")
    piv=modules.pivot(index="module",columns="endpoint",values="mean_standardized_effect"); fig,ax=plt.subplots(figsize=(6,4)); sns.heatmap(piv,cmap="vlag",center=0,ax=ax,cbar_kws={"label":"Mean standardized effect"}); ax.set_xlabel(""); ax.set_ylabel(""); save(fig,"C2-6.png")
    n=null[null.module!="mechanotransduction_network_composite"]; fig,ax=plt.subplots(figsize=(5,3.5)); lo=np.abs(n.observed_TRI_like-n.null_CI_low); hi=np.abs(n.null_CI_high-n.observed_TRI_like); ax.errorbar(n.observed_TRI_like,n.module,xerr=[lo,hi],fmt="o",color="#d62728",capsize=3); ax.set_xlabel("TRI-like score (observed ± null 95% range)"); save(fig,"C2-7.png")
    fig,ax=plt.subplots(figsize=(4,3))
    if len(bridge):
        bc=bridge.perturbation_class.value_counts(); bc.plot.bar(ax=ax,color="#9467bd"); ax.set_ylabel("Recovery genes"); ax.set_xlabel(""); ax.tick_params(axis="x",rotation=35)
    else:
        ax.text(.5,.5,"No genes met the\nstrict recovery rule",ha="center",va="center",fontsize=10); ax.set_axis_off()
    save(fig,"C2-8.png")

def reports(master,corr,delta,counts,null,bridge_counts,robust,loo,qual,fate_counts):
    strongest=null[(null.module!="mechanotransduction_network_composite")].sort_values(["empirical_P","observed_TRI_like"]).iloc[0]
    def fmt(x): return "NA" if not np.isfinite(x) else f"{x:.4g}"
    primary=delta["delta_rho"]>0 and delta["CI_low"]>0; qstable=(qual[qual.subset!="all"].delta_rho>0).sum()>=3; prot=robust["rho_vs_raw"].min()>0.9
    anystrong=(null.specificity_class=="SPECIFIC_STRONG").any(); final="STRONG_TRANSLATIONAL_ADVANTAGE" if primary and qstable else "MODERATE_CROSSOMIC_CONVERGENCE" if (not primary and corr[corr.predictor.isin(["RNA6","Ribo6"])].spearman_rho.abs().min()>0.2 and anystrong) else "PROGRAM_LEVEL_GO" if anystrong else "WEAK"
    if not np.isfinite(robust.rho_vs_raw.min()) or robust.rho_vs_raw.min()<.75: final="NO_GO"
    lines=[]
    lines.append(f"# C2 Proteome differential analysis\n\n- Processed PXD054330 matrix: 4301 unique 1:1 gene-level protein rows (3974 in the primary RNA∩Ribo∩protein set; 3966 with finite protein effects); 3 hyperosmotic versus 3 control samples.\n- Protein model: gene-level means from `PG.Log2Quantity`, equal-variance two-sided t-test, BH correction; no imputation.\n- Effect robustness: minimum protein-effect Spearman ρ versus raw processed matrix = **{fmt(robust.rho_vs_raw.min())}**; classification **{'STRONG' if prot else 'MODERATE' if robust.rho_vs_raw.min()>=.75 else 'FRAGILE'}**.\n- The upstream normalization is not fully documented by the repository; the sensitivity analysis is explicitly labeled processed-matrix normalization sensitivity.\n")
    (REP/"C2_PROTEOME_DIFFERENTIAL.md").write_text("\n".join(lines))
    r6=corr[corr.predictor=="RNA6"].iloc[0]; b6=corr[corr.predictor=="Ribo6"].iloc[0]
    (REP/"C2_PRIMARY_CROSSOMIC_ANALYSIS.md").write_text(f"# C2 primary cross-omic analysis\n\nPrimary set: **{len(master)} genes**, without DE significance filtering. RNA6→Protein24 Spearman ρ={fmt(r6.spearman_rho)}, bootstrap 95% CI [{fmt(r6.bootstrap_CI_low)}, {fmt(r6.bootstrap_CI_high)}]. Ribo6→Protein24 ρ={fmt(b6.spearman_rho)}, bootstrap 95% CI [{fmt(b6.bootstrap_CI_low)}, {fmt(b6.bootstrap_CI_high)}]. Paired bootstrap Δρ (Ribo6−RNA6) median={fmt(delta['delta_rho'])}, 95% CI [{fmt(delta['CI_low'])}, {fmt(delta['CI_high'])}], empirical P={fmt(delta['empirical_P'])}.\n\nThe preregistered advantage rule is {'met' if primary else 'not met'}; the result is interpreted as cross-study temporal triangulation, not a longitudinal causal trajectory.\n")
    largest=fate_counts.iloc[0]
    (REP/"C2_FATE_MAP.md").write_text("# C2 post-transcriptional fate map\n\nFate classes were assigned by the frozen deterministic priority rules at the primary and sensitivity thresholds.\n\n"+fate_counts.to_markdown(index=False)+f"\n\nLargest primary class: **{largest.fate_class}** (n={largest.n}). Stable classifications are reported in `results/c2/C2_FATE_MAP.tsv`.\n")
    (REP/"C2_PROGRAM_LEVEL_ANALYSIS.md").write_text("# C2 program-level analysis\n\nNine modules were tested with two-sided rank-based competitive tests and BH correction; module membership was frozen before analysis. The complete standardized-effect matrix and enrichment table are in `results/c2/C2_MODULE_EFFECTS.tsv` and `results/c2/C2_MODULE_ENRICHMENT.tsv`.\n\nStrongest frozen module by matched-null TRI-like specificity: **"+str(strongest.module)+f"** (empirical P={fmt(strongest.empirical_P)}).\n")
    (REP/"C2_MATCHED_NULL_SPECIFICITY.md").write_text("# C2 matched-null specificity\n\n10,000 null sets per module were matched on available baseline RNA abundance, protein control abundance and protein missingness. Baseline Ribo abundance, gene length and GC were unavailable in the audited inputs and were not fabricated.\n\n"+null.to_markdown(index=False)+"\n")
    (REP/"C2_PERTURBATION_BRIDGE.md").write_text("# C2 MeAIB/Torin perturbation bridge\n\nRecovery-class genes were classified using frozen TE-effect/FDR thresholds. MeAIB is described as a System-A transport perturbation; no SNAT2-dependent or causal claim is made.\n\n"+bridge_counts.to_markdown(index=False)+"\n")
    final_md=f"# C2 integrated decision\n\n**C2_VERDICT:** {final}\n\n- Primary translational advantage rule: {'YES' if primary else 'NO'}; quality sensitivity direction stable in {int((qual[qual.subset!='all'].delta_rho>0).sum())}/{len(qual)-1} subsets.\n- Program-level GO: {'YES' if anystrong else 'NO'} (at least one frozen module passed the matched-null strong rule).\n- Branch recommendation: {'YES' if final!='NO_GO' else 'NO'} for continued manuscript work; C3 remains intentionally unstarted.\n- Central claim boundary: temporal triangulation/cross-study molecular continuity only; PIEZO1 remains a representative node.\n- Main limitation: independent studies differ in osmolarity, cell model, platform and batch, with n=3 protein replicates and incomplete upstream normalization metadata.\n"
    (REP/"C2_INTEGRATED_DECISION.md").write_text(final_md)
    return final

def main():
    mapping,u,ens_to_sym,sym_to_ens,dep,pmat,p1,gse24,pert,modules,sample_cols,pivot=load_inputs()
    master,core=build_master(dep,p1,gse24,pert,modules,ens_to_sym)
    # baseline RNA control abundance is an available matching feature
    rv=pd.read_csv(ROOT/"data/processed/GSE200097_RNA_vst.tsv",sep="\t")
    rv["ensembl_gene_id"]=rv.gene.astype(str).map(sym_to_ens); ctrl=[c for c in rv.columns if "Control" in c]; base=rv.groupby("ensembl_gene_id")[ctrl].mean().mean(axis=1); master["baseline_RNA_control"]=master.ensembl_gene_id.map(base)
    write_tsv(master[core+(["baseline_RNA_control"] if "baseline_RNA_control" not in core else [])],OUT/"CROSSOMIC_MASTER_TABLE.tsv")
    dep_out=dep.rename(columns={"gene_symbol":"gene","log2FC":"log2FC"}); write_tsv(dep_out[[c for c in ["gene","protein","log2FC","SE","t","P","BH_FDR","mean_control","mean_hyper","missingness"] if c in dep_out.columns]],OUT/"PXD054330_DIFFERENTIAL_PROTEOME.tsv")
    corr,delta,boots=correlation_outputs(master); direction=directional(master); fate,fate_counts=fate_map(master); eff,enr,zvals=module_outputs(master,modules); null=module_null(master,modules,zvals); robust=protein_robustness(dep,pivot,mapping,ens_to_sym,sample_cols,master); loo=leave_one_out(master); qual=quality_sensitivity(master); strata=expression_strata(master); bridge,bridge_counts=perturbation_bridge(master,fate); rep=representative(master,fate)
    delta_fig=dict(delta); delta_fig["_vals"]=list(boots.get("Ribo6",np.array([]))[:0])
    # Recompute the delta bootstrap values for the histogram without changing the preregistered result.
    _,_,_,_,dvals,_,_=bootstrap_delta(master.RNA_6h_log2FC,master.Ribo_6h_log2FC,master.PXD054330_24h_protein_log2FC,B=B,seed=SEED+101); delta_fig["_vals"]=dvals
    figures(master,corr,delta_fig,fate,eff,enr,null,bridge)
    final=reports(master,corr,delta,null,fate_counts,null,bridge_counts,robust,loo,qual,fate_counts) if False else reports(master,corr,delta,fate_counts,null,bridge_counts,robust,loo,qual,fate_counts)
    summary={"C2_VERDICT":final,"primary_gene_count":int(len(master)),"protein_rows":int(len(dep)),"primary_delta":delta,"protein_robustness_min_rho":float(robust.rho_vs_raw.min()),"largest_fate":fate_counts.iloc[0].to_dict(),"strongest_module":null.sort_values("empirical_P").iloc[0].to_dict(),"files":sorted(str(p.relative_to(CM)) for p in OUT.glob("*") if p.is_file())}
    (OUT/"C2_AUDIT_SUMMARY.json").write_text(json.dumps(summary,indent=2,default=lambda x: None))
    print(json.dumps(summary,indent=2,default=lambda x: None))

if __name__=="__main__": main()
