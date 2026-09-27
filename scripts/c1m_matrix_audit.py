#!/usr/bin/env python3
import csv
import hashlib
import json
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PIVOT = ROOT / "data/raw/proteome/PXD054330/processed/Orbitrap_HCEC_DED_vs_NRM_Report_BGS_Analysis_Grid_View_Report_Pivot.tsv"
FASTA = ROOT / "data/raw/proteome/PXD054330/processed/20231227_uniprotkb_homo_sapiens_Review_ISO.fasta"
DATE = "2026-09-21"
for rel in ["metadata/c1m", "results/c1m", "logs/c1m"]:
    (ROOT / rel).mkdir(parents=True, exist_ok=True)

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def write_tsv(path, fields, rows):
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def val(x):
    if x is None or str(x).strip().lower() in {"", "nan", "na", "null"}:
        return None
    try:
        y = float(x)
        return y if math.isfinite(y) else None
    except ValueError:
        return None

def canonical(a):
    return re.sub(r"-\d+$", "", a.strip())

with PIVOT.open(newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    fields = reader.fieldnames or []
    rows = list(reader)
sample_cols = [f for f in fields if ".raw.PG.Log2Quantity" in f]
assert len(sample_cols) == 6, sample_cols

specs = [
    ("DED1_IS_DIA.raw", "UNKNOWN", 4032662296, "md5:bc6c6d63bee027bff73f8ded1d0786b5", "NO", "Thermo vendor RAW; hyperosmotic run."),
    ("DED2_IS_DIA.raw", "UNKNOWN", 4172725165, "md5:e27e01cbe50da097276bd89cbbf9e147", "NO", "Thermo vendor RAW; hyperosmotic run."),
    ("DED4_IS_DIA.raw", "UNKNOWN", 3898013318, "md5:e55d80b313ac70e50a59f2d62ba17091", "NO", "Thermo vendor RAW; hyperosmotic run."),
    ("NR1_IS_DIA_20220730051416.raw", "UNKNOWN", 3776299204, "md5:8ac68b1b7d586339ed653e79690b896c", "NO", "Thermo vendor RAW; isosmotic run."),
    ("NR3_IS_DIA.raw", "UNKNOWN", 4146914593, "md5:b38dc70f2807689183f5742b471344d9", "NO", "Thermo vendor RAW; isosmotic run."),
    ("NR4_IS_DIA.raw", "UNKNOWN", 4206950049, "md5:555ca45f8a1bafc3a13f7b84df95d980", "NO", "Thermo vendor RAW; isosmotic run."),
    ("Orbitrap_HCEC_DED_vs_NRM_Report_BGS Analysis Grid View Report (Pivot).tsv", "PROTEIN_MATRIX", 1194868, "sha256:" + sha256(PIVOT), "YES", "Tab-separated Spectronaut pivot export; downloaded unchanged."),
    ("20231227_uniprotkb_homo_sapiens_Review_ISO.fasta", "METADATA", 29232249, "sha256:" + sha256(FASTA), "YES", "Reviewed human reference FASTA; downloaded unchanged."),
    ("20240129_160649_Orbitrap_HCEC_DED_vs_NRM.sne", "SEARCH_OUTPUT", 3217281567, "md5:6122e8469d878c8432b0b6031fa1d6c3", "NO", "Spectronaut SNE/library output; not needed for matrix audit."),
]
inventory = []
for name, level, size, checksum, downloaded, note in specs:
    if name.endswith(".raw"):
        ftype, desc = "Thermo vendor RAW", "Thermo Fisher vendor RAW DIA file."
    elif name.endswith(".fasta"):
        ftype, desc = "reference FASTA", "Reviewed Homo sapiens UniProt reference FASTA."
    elif name.endswith(".sne"):
        ftype, desc = "Spectronaut search/library output", "Spectronaut SNE/library file."
    else:
        ftype, desc = "tabular protein-group abundance export", "Spectronaut pivot with six raw-run quantity columns."
    inventory.append({"file_name": name, "source_repository": "jPOST JPST003233 / OmicsDI PXD054330", "file_type": ftype, "reported_size": size, "download_url_or_identifier": "https://storage.jpostdb.org/JPST003233/" + quote(name, safe="._-"), "description": desc, "candidate_level": level, "downloaded": downloaded, "sha256": checksum.replace("sha256:", "") if checksum.startswith("sha256:") else "NOT_COMPUTED_PUBLIC_MD5_ONLY", "notes": note})
write_tsv(ROOT / "metadata/c1m/PXD054330_PROCESSED_FILE_INVENTORY.tsv", ["file_name", "source_repository", "file_type", "reported_size", "download_url_or_identifier", "description", "candidate_level", "downloaded", "sha256", "notes"], inventory)

def sample_info(col):
    raw = re.search(r"\[\d+\] (.+?\.raw)\.PG\.Log2Quantity$", col).group(1)
    if raw.startswith("DED"):
        condition = "hyperosmotic"
        rep = raw.split("_IS_DIA")[0].replace("DED", "")
        sid = "DED" + rep
    else:
        condition = "isosmotic"
        rep = raw.split("_IS_DIA")[0].replace("NR", "")
        sid = "NR" + rep
    return raw, sid, condition, rep

sample_map = []
for col in sample_cols:
    raw, sid, condition, rep = sample_info(col)
    sample_map.append({"matrix_column": col, "sample_id": sid, "raw_file": raw, "condition": condition, "biological_replicate": rep, "technical_replicate": "NOT_REPORTED", "mapping_source": "Exact matrix raw basename + jPOST file list; replicate index from filename; publication states biological triplicates.", "mapping_confidence": "INFERRED"})
write_tsv(ROOT / "metadata/c1m/PXD054330_SAMPLE_MAP.tsv", ["matrix_column", "sample_id", "raw_file", "condition", "biological_replicate", "technical_replicate", "mapping_source", "mapping_confidence"], sample_map)

values = {c: [val(r[c]) for r in rows] for c in sample_cols}
total_cells = len(rows) * len(sample_cols)
missing = sum(v is None for vs in values.values() for v in vs)
complete_idx = [i for i in range(len(rows)) if all(values[c][i] is not None for c in sample_cols)]
complete = np.array([[values[c][i] for c in sample_cols] for i in complete_idx], dtype=float)
cor = np.corrcoef(complete, rowvar=False) if len(complete_idx) >= 2 else np.full((6, 6), np.nan)
pca = np.full((6, 2), np.nan)
if len(complete_idx) >= 2:
    centered = complete - complete.mean(axis=0, keepdims=True)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    pca = centered @ vt[:2].T

def q(vals, p):
    return float(np.quantile(np.asarray(vals, dtype=float), p))

sample_qc = []
for j, c in enumerate(sample_cols):
    vals = [x for x in values[c] if x is not None]
    nmiss = len(rows) - len(vals)
    peers = [cor[j, k] for k in range(6) if k != j and math.isfinite(cor[j, k])]
    sample_qc.append({"sample_id": sample_map[j]["sample_id"], "matrix_column": c, "raw_file": sample_map[j]["raw_file"], "condition": sample_map[j]["condition"], "biological_replicate": sample_map[j]["biological_replicate"], "detected_protein_groups": len(vals), "missing_count": nmiss, "missing_fraction": f"{nmiss / len(rows):.8f}", "total_log2_quantity": f"{sum(vals):.8f}", "mean_log2_quantity": f"{statistics.mean(vals):.8f}", "median_log2_quantity": f"{statistics.median(vals):.8f}", "iqr_log2_quantity": f"{q(vals,.75)-q(vals,.25):.8f}", "mean_pairwise_complete_case_pearson": f"{statistics.mean(peers):.8f}", "pca1_complete_case": f"{pca[j,0]:.8f}", "pca2_complete_case": f"{pca[j,1]:.8f}"})
write_tsv(ROOT / "results/c1m/PXD054330_SAMPLE_QC.tsv", list(sample_qc[0]), sample_qc)

missing_rows = [{"level": "OVERALL", "feature": "ALL_SAMPLE_CELLS", "n_total": total_cells, "n_missing": missing, "missing_fraction": f"{missing / total_cells:.8f}", "category": "PARTIAL_MISSING" if missing else "COMPLETE_CASE"}]
for c in sample_cols:
    nmiss = len(rows) - sum(v is not None for v in values[c])
    frac = nmiss / len(rows)
    missing_rows.append({"level": "SAMPLE", "feature": c, "n_total": len(rows), "n_missing": nmiss, "missing_fraction": f"{frac:.8f}", "category": "COMPLETE_CASE" if not nmiss else ("SEVERE_MISSING" if frac >= .2 else "PARTIAL_MISSING")})
for i, r in enumerate(rows):
    nmiss = sum(values[c][i] is None for c in sample_cols)
    frac = nmiss / 6
    missing_rows.append({"level": "PROTEIN_GROUP", "feature": r["PG.ProteinAccessions"], "n_total": 6, "n_missing": nmiss, "missing_fraction": f"{frac:.8f}", "category": "COMPLETE_CASE" if not nmiss else ("SEVERE_MISSING" if frac >= .2 else "PARTIAL_MISSING")})
write_tsv(ROOT / "results/c1m/PXD054330_MISSINGNESS.tsv", ["level", "feature", "n_total", "n_missing", "missing_fraction", "category"], missing_rows)

matrix_qc = [
    {"metric": "matrix_classification", "value": "FULL_SAMPLE_MATRIX", "notes": "Sample-level protein-group abundance export; six individual raw-run columns; not DEP-only."},
    {"metric": "row_count", "value": len(rows), "notes": "Rows represent Spectronaut protein groups."},
    {"metric": "column_count", "value": len(fields), "notes": "Seven metadata/summary columns plus six individual sample columns."},
    {"metric": "sample_column_count", "value": len(sample_cols), "notes": "Six individual columns named by raw file."},
    {"metric": "protein_group_count", "value": len(rows), "notes": "One row per exported protein group."},
    {"metric": "unique_accession_token_count", "value": len(set(a for r in rows for a in r["PG.ProteinAccessions"].split(";") if a)), "notes": "Includes isoform tokens."},
    {"metric": "unique_canonical_accession_count", "value": len(set(canonical(a) for r in rows for a in r["PG.ProteinAccessions"].split(";") if a)), "notes": "Isoform suffixes collapsed for mapping."},
    {"metric": "candidate_flag_true_rows", "value": sum(r["PG.IsCandidate"] == "True" for r in rows), "notes": "Candidate flag is a subset column; the export contains 4,743 rows, not only 598 published DEPs."},
    {"metric": "complete_case_rows", "value": len(complete_idx), "notes": "No imputation; complete rows used for QC correlation/PCA only."},
    {"metric": "overall_missingness", "value": f"{missing / total_cells:.8f}", "notes": "Missing sample-level cells divided by all sample-level cells."},
    {"metric": "max_sample_missingness", "value": f"{max(len(rows) - sum(v is not None for v in values[c]) for c in sample_cols) / len(rows):.8f}", "notes": "Maximum among six runs."},
    {"metric": "quantification_field", "value": "PG.Log2Quantity", "notes": "Value type follows the deposited column name; no normalization inference."},
    {"metric": "normalization", "value": "UNKNOWN", "notes": "Exact Spectronaut normalization settings are not exposed in the pivot."},
    {"metric": "missing_value_handling", "value": "UNKNOWN; NO IMPUTATION", "notes": "Missingness audited; imputation deferred to C2."},
    {"metric": "differential_proteome_feasibility", "value": "PASS", "notes": "Three runs per condition permit condition means, log2FC and basic uncertainty; no cross-omic correlation run."},
]
write_tsv(ROOT / "results/c1m/PXD054330_MATRIX_QC.tsv", ["metric", "value", "notes"], matrix_qc)

# C1M.9-10: UniProt primary gene names followed by MyGene human symbol -> Ensembl.

def primary_gene(names):
    return names.split()[0] if names else ""

all_tokens = sorted(set(a for r in rows for a in r["PG.ProteinAccessions"].split(";") if a))
canonical_tokens = sorted(set(canonical(a) for a in all_tokens))
import time
import requests
http = requests.Session()
uniprot = {}
for start in range(0, len(canonical_tokens), 100):
    batch = canonical_tokens[start:start + 100]
    res = http.get("https://rest.uniprot.org/uniprotkb/search", params={"query": "accession:(" + " OR ".join(batch) + ")", "format": "tsv", "fields": "accession,gene_names", "size": 500}, timeout=90)
    res.raise_for_status()
    for line in res.text.splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) >= 2:
            uniprot[parts[0]] = (primary_gene(parts[1]), parts[1])
    time.sleep(.05)

symbols = sorted(set(v[0] for v in uniprot.values() if v[0]))
symbol_ens = defaultdict(set)
for start in range(0, len(symbols), 100):
    batch = symbols[start:start + 100]
    res = http.post("https://mygene.info/v3/query", data={"q": ",".join(batch), "scopes": "symbol", "species": "human", "fields": "ensembl.gene,symbol", "size": 100}, timeout=90)
    res.raise_for_status()
    payload = res.json() if isinstance(res.json(), list) else []
    for hit in payload:
        query = hit.get("query", "")
        symbol = hit.get("symbol", "")
        if query != symbol:
            continue
        ens = hit.get("ensembl", {}).get("gene", "") if isinstance(hit.get("ensembl"), dict) else ""
        if isinstance(ens, str) and ens:
            symbol_ens[query].add(ens)
        elif isinstance(ens, list):
            symbol_ens[query].update(x for x in ens if x)
    time.sleep(.05)

def decoy(s):
    return bool(re.search(r"(^|[;_])(DECOY|REV)[_;]", s, re.I)) or "DECOY" in s.upper() or "REVERSE" in s.upper()

def contaminant(s):
    return bool(re.search(r"(^|[;_])(CON__|CONTAM|CONT_)", s, re.I))

id_audit = []
mapping = []
for r in rows:
    group = r["PG.ProteinAccessions"]
    tok = [x for x in group.split(";") if x]
    canon = [canonical(x) for x in tok]
    syms = sorted(set(uniprot.get(x, ("", ""))[0] for x in canon if uniprot.get(x, ("", ""))[0]))
    ens = sorted(set(e for s in syms for e in symbol_ens.get(s, set())))
    is_decoy = decoy(group) or decoy(r["PG.ProteinDescriptions"])
    is_contam = contaminant(group) or contaminant(r["PG.ProteinDescriptions"])
    if is_decoy or is_contam or not syms or not ens:
        mtype = "UNMAPPED"
    elif len(ens) == 1:
        mtype = "UNIQUE_1TO1"
    elif len(tok) == 1:
        mtype = "ONE_PROTEIN_MULTI_GENE"
    else:
        mtype = "PROTEIN_GROUP_MULTI_GENE"
    id_audit.append({"protein_group": group, "n_accession_tokens": len(tok), "n_unique_canonical_accessions": len(set(canon)), "n_primary_gene_symbols": len(syms), "n_ensembl_gene_ids": len(ens), "decoy_flag": "YES" if is_decoy else "NO", "contaminant_flag": "YES" if is_contam else "NO", "missing_gene_annotation": "YES" if not syms or not ens else "NO", "mapping_type": mtype})
    mapping.append({"original_protein_identifier": group, "uniprot_accession": ";".join(sorted(set(canon))), "gene_symbol": ";".join(syms), "ensembl_gene_id": ";".join(ens), "mapping_type": mtype, "n_group_accessions": len(tok), "decoy_flag": "YES" if is_decoy else "NO", "contaminant_flag": "YES" if is_contam else "NO", "mapping_source": "UniProt REST gene_names + MyGene.info human symbol-to-Ensembl; retrieved 2026-09-21.", "mapping_rule": "Primary set retains only one Ensembl gene ID per protein group; multi-gene groups are sensitivity-only."})
write_tsv(ROOT / "metadata/c1m/PXD054330_PROTEIN_ID_AUDIT.tsv", ["protein_group", "n_accession_tokens", "n_unique_canonical_accessions", "n_primary_gene_symbols", "n_ensembl_gene_ids", "decoy_flag", "contaminant_flag", "missing_gene_annotation", "mapping_type"], id_audit)
write_tsv(ROOT / "metadata/c1m/PROTEIN_GENE_ID_MAP_REAL.tsv", ["original_protein_identifier", "uniprot_accession", "gene_symbol", "ensembl_gene_id", "mapping_type", "n_group_accessions", "decoy_flag", "contaminant_flag", "mapping_source", "mapping_rule"], mapping)

def genes(path):
    with path.open(newline="") as fh:
        return {r[0] for r in csv.reader(fh, delimiter="\t") if r and r[0] not in {"gene", ""}}

rna = genes(ROOT.parent / "data/processed/GSE200097_RNA_vst.tsv")
ribo = genes(ROOT.parent / "results/phase1/PHASE1E_RNA_RIBO_ALL.tsv")
gse323 = genes(ROOT.parent / "data/processed/GSE323164_RNA_vst.tsv")
protein_symbols = {m["gene_symbol"].split(";")[0] for m in mapping if m["mapping_type"] == "UNIQUE_1TO1" and m["gene_symbol"]}
overlap = [
    ("RNA_Ribo", "GSE200097_RNA", "GSE200097_Ribo", "", "", len(rna & ribo), "KNOWN", "Native gene-symbol intersection."),
    ("RNA_Protein", "GSE200097_RNA", "PXD054330_protein", "", "", len(rna & protein_symbols), "KNOWN", "Protein side restricted to UNIQUE_1TO1."),
    ("Ribo_Protein", "GSE200097_Ribo", "PXD054330_protein", "", "", len(ribo & protein_symbols), "KNOWN", "Protein side restricted to UNIQUE_1TO1."),
    ("RNA_Ribo_Protein", "GSE200097_RNA", "GSE200097_Ribo", "PXD054330_protein", "", len(rna & ribo & protein_symbols), "KNOWN", "Coverage only; no C2 inference."),
    ("RNA_Ribo_GSE323164_Protein", "GSE200097_RNA", "GSE200097_Ribo", "GSE323164_RNA", "PXD054330_protein", len(rna & ribo & gse323 & protein_symbols), "KNOWN", "Cross-study feature coverage; not longitudinal."),
]
overlap_rows = [{"overlap_label": a, "dataset_a": b, "dataset_b": c, "dataset_c": d, "dataset_d": e, "n_features": n, "identifier_rule": "Native transcriptomic gene symbols; protein UNIQUE_1TO1 gene symbols audited to Ensembl", "status": s, "notes": note} for a,b,c,d,e,n,s,note in overlap]
write_tsv(ROOT / "results/c1m/CROSSOMIC_FEATURE_OVERLAP_REAL.tsv", ["overlap_label", "dataset_a", "dataset_b", "dataset_c", "dataset_d", "n_features", "identifier_rule", "status", "notes"], overlap_rows)

summary = {"pivot_sha256": sha256(PIVOT), "fasta_sha256": sha256(FASTA), "matrix_rows": len(rows), "matrix_columns": len(fields), "sample_columns": sample_cols, "missing_cells": missing, "missing_fraction": missing / total_cells, "complete_case_rows": len(complete_idx), "unique_1to1_groups": sum(x["mapping_type"] == "UNIQUE_1TO1" for x in mapping), "multi_gene_groups": sum(x["mapping_type"] == "PROTEIN_GROUP_MULTI_GENE" for x in mapping), "one_protein_multi_gene": sum(x["mapping_type"] == "ONE_PROTEIN_MULTI_GENE" for x in mapping), "unmapped_groups": sum(x["mapping_type"] == "UNMAPPED" for x in mapping), "unique_primary_gene_symbols": len(protein_symbols), "rna_ribo_protein": len(rna & ribo & protein_symbols), "four_way": len(rna & ribo & gse323 & protein_symbols), "queried_canonical_accessions": len(canonical_tokens), "uniprot_accessions_with_gene": len(uniprot), "queried_gene_symbols": len(symbols), "gene_symbols_with_ensembl": sum(bool(symbol_ens.get(s)) for s in symbols)}
(ROOT / "logs/c1m/C1M_AUDIT_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
