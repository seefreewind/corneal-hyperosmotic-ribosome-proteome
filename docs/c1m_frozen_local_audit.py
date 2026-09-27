#!/usr/bin/env python3
"""C4 clean-environment C1M audit using frozen processed inputs.

The original C1M script can refresh UniProt/MyGene mappings. The C4 rerun
uses the already frozen mapping tables so network availability cannot change
the manuscript-critical feature space.
"""
from pathlib import Path
import hashlib, json
import pandas as pd

HERE = Path(__file__).resolve()
CM = HERE.parents[2]
ROOT = CM.parent
PIVOT = CM / "data/raw/proteome/PXD054330/processed/Orbitrap_HCEC_DED_vs_NRM_Report_BGS_Analysis_Grid_View_Report_Pivot.tsv"
MAP = CM / "metadata/c1m/PROTEIN_GENE_ID_MAP_REAL.tsv"

def sha256(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()

pivot=pd.read_csv(PIVOT,sep="\t")
sample_cols=[c for c in pivot.columns if ".raw.PG.Log2Quantity" in c]
complete=int(pivot[sample_cols].notna().all(axis=1).sum())
missing=float(pivot[sample_cols].isna().mean().mean())
mapping=pd.read_csv(MAP,sep="\t")
unique=mapping[mapping.mapping_type.eq("UNIQUE_1TO1")].copy()
protein_symbols=set(unique.gene_symbol.dropna().astype(str).str.split(";").str[0])
rna=pd.read_csv(ROOT/"data/processed/GSE200097_RNA_vst.tsv",sep="\t",usecols=["gene"])
ribo=pd.read_csv(ROOT/"results/phase1/PHASE1E_RNA_RIBO_ALL.tsv",sep="\t",usecols=["gene"])
gse=pd.read_csv(ROOT/"data/processed/GSE323164_RNA_vst.tsv",sep="\t",usecols=["gene"])
rna_set=set(rna.gene.astype(str)); ribo_set=set(ribo.gene.astype(str)); gse_set=set(gse.gene.astype(str))
out={"mode":"FROZEN_PROCESSED_INPUTS_NO_NETWORK_REFRESH","pivot_sha256":sha256(PIVOT),"pivot_rows":int(len(pivot)),"pivot_columns":int(len(pivot.columns)),"sample_columns":int(len(sample_cols)),"complete_case_rows":complete,"overall_missingness":missing,"unique_1to1_groups":int((mapping.mapping_type=="UNIQUE_1TO1").sum()),"multi_gene_groups":int(mapping.mapping_type.isin(["ONE_PROTEIN_MULTI_GENE","PROTEIN_GROUP_MULTI_GENE"]).sum()),"unmapped_groups":int((mapping.mapping_type=="UNMAPPED").sum()),"protein_gene_symbols":len(protein_symbols),"rna_ribo_protein":len(rna_set&ribo_set&protein_symbols),"four_dataset":len(rna_set&ribo_set&gse_set&protein_symbols),"expected_primary_overlap":3974,"expected_four_dataset":3950,"PASS":bool(len(rna_set&ribo_set&protein_symbols)==3974 and len(rna_set&ribo_set&gse_set&protein_symbols)==3950 and complete==4518)}
(HERE.parent/"C1M_frozen_local_summary.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))
