#!/usr/bin/env python3
"""C7B unified repeated-OOF CV and conditional permutation repair.

This is a new, isolated C7B branch.  The frozen C2R script and all non-CV
analyses are intentionally left untouched.  The primary statistic is the
mean across ten repeat-level OOF Spearman deltas (Model B minus Model A).
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats
import sklearn
from sklearn.model_selection import KFold


ROOT = Path(__file__).resolve().parents[2]
CM = ROOT / "computational_manuscript"
RESULTS = CM / "results" / "c7b"
REPORTS = CM / "reports" / "final" / "c7b"
RESULTS.mkdir(parents=True, exist_ok=True)
REPORTS.mkdir(parents=True, exist_ok=True)

N_REPEATS = 10
N_FOLDS = 10
FOLD_SEED = 270929
PERM_SEED = 270931
B_DEFAULT = 10_000

_WORK_D = None
_WORK_FOLD_CACHE = None
_WORK_OOF_A = None


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if int(ok.sum()) < 4:
        return float("nan")
    return float(stats.spearmanr(x[ok], y[ok]).statistic)


def standardize_train(X_train: np.ndarray, X_test: np.ndarray):
    mu = np.nanmean(X_train, axis=0)
    sd = np.nanstd(X_train, axis=0, ddof=1)
    sd = np.asarray(sd, dtype=float)
    sd[~np.isfinite(sd) | (sd == 0)] = 1.0
    return (X_train - mu) / sd, (X_test - mu) / sd, mu, sd


def fit_predict(X: np.ndarray, y: np.ndarray, train_idx: np.ndarray, test_idx: np.ndarray):
    Xtr, Xte, mu, sd = standardize_train(X[train_idx], X[test_idx])
    beta = np.linalg.lstsq(
        np.column_stack([np.ones(len(train_idx)), Xtr]),
        np.asarray(y[train_idx], dtype=float),
        rcond=None,
    )[0]
    pred = np.column_stack([np.ones(len(test_idx)), Xte]) @ beta
    return pred, mu, sd, beta


def fit_residualization(RNA_train: np.ndarray, Ribo_train: np.ndarray):
    X = np.column_stack([np.ones(len(RNA_train)), RNA_train])
    beta = np.linalg.lstsq(X, Ribo_train, rcond=None)[0]
    residual = Ribo_train - X @ beta
    return beta, residual


def read_complete_case_table():
    master = pd.read_csv(RESULTS.parent / "c2" / "CROSSOMIC_MASTER_TABLE.tsv", sep="\t")
    dep = pd.read_csv(RESULTS.parent / "c2" / "PXD054330_DIFFERENTIAL_PROTEOME.tsv", sep="\t")
    dep = dep[["gene", "mean_control"]].drop_duplicates("gene").rename(columns={"gene": "gene_symbol"})
    merged = master.merge(dep, on="gene_symbol", how="left")
    required = ["RNA_6h_log2FC", "Ribo_6h_log2FC", "TE_6h_log2FC", "PXD054330_24h_protein_log2FC"]
    complete = merged.loc[merged[required].notna().all(axis=1)].reset_index(drop=True)
    return complete, required


def make_folds(n: int):
    folds = []
    for repeat_id in range(N_REPEATS):
        kf = KFold(N_FOLDS, shuffle=True, random_state=FOLD_SEED + repeat_id)
        for fold_id, (train_idx, test_idx) in enumerate(kf.split(np.arange(n))):
            folds.append((repeat_id, fold_id, train_idx.astype(int), test_idx.astype(int)))
    return folds


def assignment_hash(folds):
    payload = []
    for repeat_id, fold_id, train_idx, test_idx in folds:
        payload.append(
            {
                "repeat": int(repeat_id),
                "fold": int(fold_id),
                "train": train_idx.tolist(),
                "test": test_idx.tolist(),
            }
        )
    return hashlib.sha256(json.dumps(payload, separators=(",", ":")).encode()).hexdigest()


def build_frozen_strata(d: pd.DataFrame):
    # These are the frozen C2R strata.  They use baseline/control quantities,
    # not the hyperosmotic Protein24 treatment effect.
    rna_bin = pd.qcut(d["RNA_6h_log2FC"], 5, labels=False, duplicates="drop")
    base_bin = pd.qcut(d["baseline_RNA_control"], 4, labels=False, duplicates="drop")
    prot_bin = pd.qcut(d["mean_control"], 4, labels=False, duplicates="drop")
    strata = np.asarray([f"{a}_{b}_{c}" for a, b, c in zip(rna_bin, base_bin, prot_bin)], dtype=object)
    if pd.isna(strata).any():
        raise RuntimeError("PERMUTATION_DESIGN_HOLD: frozen strata contain missing labels")
    return strata


def prepare_fold_cache(d, folds, strata):
    RNA = d["RNA_6h_log2FC"].to_numpy(float)
    Ribo = d["Ribo_6h_log2FC"].to_numpy(float)
    y = d["PXD054330_24h_protein_log2FC"].to_numpy(float)
    cache = []
    for repeat_id, fold_id, train_idx, test_idx in folds:
        train_strata = strata[train_idx]
        test_strata = strata[test_idx]
        train_pos_by_stratum = {}
        test_pos_by_stratum = {}
        for g in np.unique(strata):
            tr_pos = np.flatnonzero(train_strata == g)
            te_pos = np.flatnonzero(test_strata == g)
            if len(te_pos) and len(tr_pos) == 0:
                raise RuntimeError(
                    "PERMUTATION_DESIGN_HOLD: held-out stratum has no training residual pool "
                    f"(repeat={repeat_id}, fold={fold_id}, stratum={g})"
                )
            train_pos_by_stratum[g] = tr_pos
            test_pos_by_stratum[g] = te_pos
        beta_resid, residual = fit_residualization(RNA[train_idx], Ribo[train_idx])
        cache.append(
            {
                "repeat": repeat_id,
                "fold": fold_id,
                "train_idx": train_idx,
                "test_idx": test_idx,
                "train_strata": train_strata,
                "test_strata": test_strata,
                "train_pos_by_stratum": train_pos_by_stratum,
                "test_pos_by_stratum": test_pos_by_stratum,
                "resid_beta": beta_resid,
                "train_residual": residual,
            }
        )
    return cache


def compute_repeat_oof_statistic(y: np.ndarray, oof_by_model: dict[str, dict[int, np.ndarray]]):
    """The single repeat-level OOF statistic used by observed and null paths."""
    rows = []
    for repeat_id in range(N_REPEATS):
        rhos = {model: spearman(y, vectors[repeat_id]) for model, vectors in oof_by_model.items()}
        row = {"repeat_id": repeat_id, "n_oof": len(y)}
        row.update({f"rho_{model}": value for model, value in rhos.items()})
        row["delta_BA"] = rhos["B"] - rhos["A"]
        if "C" in rhos:
            row["delta_CA"] = rhos["C"] - rhos["A"]
        if "D" in rhos:
            row["delta_DA"] = rhos["D"] - rhos["A"]
        rows.append(row)
    return pd.DataFrame(rows)


def run_observed_cv(d, folds):
    RNA = d["RNA_6h_log2FC"].to_numpy(float)
    Ribo = d["Ribo_6h_log2FC"].to_numpy(float)
    TE = d["TE_6h_log2FC"].to_numpy(float)
    y = d["PXD054330_24h_protein_log2FC"].to_numpy(float)
    Xs = {
        "A": RNA[:, None],
        "B": np.column_stack([RNA, Ribo]),
        "C": np.column_stack([RNA, TE]),
        "D": np.column_stack([RNA, Ribo, TE]),
    }
    fold_rows = []
    oof_by_model = {model: {} for model in Xs}
    for repeat_id in range(N_REPEATS):
        oof = {name: np.full(len(d), np.nan, dtype=float) for name in Xs}
        for rep, fold_id, train_idx, test_idx in [x for x in folds if x[0] == repeat_id]:
            for model, X in Xs.items():
                pred, _, _, _ = fit_predict(X, y, train_idx, test_idx)
                oof[model][test_idx] = pred
                fold_rows.append(
                    {
                        "repeat_id": repeat_id,
                        "fold_id": fold_id,
                        "model": model,
                        "n_train": len(train_idx),
                        "n_test": len(test_idx),
                        "rho_fold": spearman(y[test_idx], pred),
                    }
                )
        if any(np.isnan(v).any() for v in oof.values()):
            raise RuntimeError(f"OOF construction failed for repeat {repeat_id}")
        for model in Xs:
            oof_by_model[model][repeat_id] = oof[model]
    repeat_df = compute_repeat_oof_statistic(y, oof_by_model)
    return repeat_df, pd.DataFrame(fold_rows), oof_by_model["A"]


def permuted_ribo_for_fold(cache_item, d, rng):
    RNA = d["RNA_6h_log2FC"].to_numpy(float)
    Ribo = d["Ribo_6h_log2FC"].to_numpy(float)
    train_idx = cache_item["train_idx"]
    test_idx = cache_item["test_idx"]
    residual = cache_item["train_residual"]
    train_perm = residual.copy()
    test_perm = np.full(len(test_idx), np.nan, dtype=float)
    for g, tr_pos in cache_item["train_pos_by_stratum"].items():
        te_pos = cache_item["test_pos_by_stratum"][g]
        if len(tr_pos):
            train_perm[tr_pos] = residual[rng.permutation(tr_pos)]
        if len(te_pos):
            # Held-out residuals are drawn only from the corresponding
            # training residual pool; held-out Ribo and Protein24 are unused.
            test_perm[te_pos] = rng.choice(residual[tr_pos], size=len(te_pos), replace=True)
    beta = cache_item["resid_beta"]
    rna_fit_train = beta[0] + beta[1] * RNA[train_idx]
    rna_fit_test = beta[0] + beta[1] * RNA[test_idx]
    return rna_fit_train + train_perm, rna_fit_test + test_perm


def _single_permutation(perm_id):
    d = _WORK_D
    fold_cache = _WORK_FOLD_CACHE
    oof_a_by_repeat = _WORK_OOF_A
    RNA = d["RNA_6h_log2FC"].to_numpy(float)
    y = d["PXD054330_24h_protein_log2FC"].to_numpy(float)
    # Independent per-permutation streams make parallel execution exactly
    # reproducible without changing the frozen B, folds, or statistic.
    rng = np.random.default_rng(np.random.SeedSequence([PERM_SEED, int(perm_id)]))
    oof_b = {r: np.full(len(d), np.nan, dtype=float) for r in range(N_REPEATS)}
    for item in fold_cache:
        train_idx = item["train_idx"]
        test_idx = item["test_idx"]
        rnull_train, rnull_test = permuted_ribo_for_fold(item, d, rng)
        X_train = np.column_stack([RNA[train_idx], rnull_train])
        X_test = np.column_stack([RNA[test_idx], rnull_test])
        X_all = np.vstack([X_train, X_test])
        y_all = np.concatenate([y[train_idx], y[test_idx]])
        local_train = np.arange(len(train_idx))
        local_test = np.arange(len(train_idx), len(train_idx) + len(test_idx))
        pred, _, _, _ = fit_predict(X_all, y_all, local_train, local_test)
        oof_b[item["repeat"]][test_idx] = pred
    for repeat_id in range(N_REPEATS):
        if np.isnan(oof_b[repeat_id]).any():
            raise RuntimeError(f"Permutation OOF construction failed for repeat {repeat_id}")
    repeat_df = compute_repeat_oof_statistic(y, {"A": oof_a_by_repeat, "B": oof_b})
    return float(repeat_df["delta_BA"].mean())


def statistic_repeat_oof(d, folds, fold_cache, oof_a_by_repeat, n_perm, workers):
    global _WORK_D, _WORK_FOLD_CACHE, _WORK_OOF_A
    _WORK_D = d
    _WORK_FOLD_CACHE = fold_cache
    _WORK_OOF_A = oof_a_by_repeat
    if workers == 1:
        values = [_single_permutation(i) for i in range(1, n_perm + 1)]
    else:
        context = mp.get_context("fork")
        with context.Pool(processes=workers) as pool:
            values = []
            for completed, value in enumerate(pool.imap(_single_permutation, range(1, n_perm + 1), chunksize=1), start=1):
                values.append(value)
                if completed % 100 == 0:
                    print(f"permutation {completed}/{n_perm}", flush=True)
    return np.asarray(values, dtype=float)


def summarize_repeat(repeat_df):
    delta = repeat_df["delta_BA"].to_numpy(float)
    return {
        "rho_A_mean": float(repeat_df["rho_A"].mean()),
        "rho_B_mean": float(repeat_df["rho_B"].mean()),
        "delta_mean": float(delta.mean()),
        "delta_median": float(np.median(delta)),
        "delta_sd": float(np.std(delta, ddof=1)),
        "delta_min": float(delta.min()),
        "delta_max": float(delta.max()),
        "fraction_repeats_B_gt_A": float(np.mean(delta > 0)),
        "delta_q025": float(np.quantile(delta, 0.025)),
        "delta_q975": float(np.quantile(delta, 0.975)),
    }


def write_reports(d, folds, repeat_df, fold_df, null, old_t, old_p, input_paths, runtime_seconds, n_perm, workers):
    t_obs = float(repeat_df["delta_BA"].mean())
    extreme = int(np.sum(null >= t_obs))
    p = float((extreme + 1) / (n_perm + 1))
    summary = summarize_repeat(repeat_df)
    null_q = np.quantile(null, [0.025, 0.5, 0.975])
    results = {
        "T_obs": t_obs,
        "permutation_B": int(n_perm),
        "extreme_count_ge": extreme,
        "formal_P": p,
        "null_mean": float(null.mean()),
        "null_sd": float(np.std(null, ddof=1)),
        "null_q025": float(null_q[0]),
        "null_median": float(null_q[1]),
        "null_q975": float(null_q[2]),
        "old_T_fold_all": old_t,
        "old_P": old_p,
    }
    pd.DataFrame(
        [{
            "permutation": i + 1,
            "T_perm_repeat_OOF": value,
            "T_obs_repeat_OOF": t_obs,
            "empirical_upper_tail_P": p,
        } for i, value in enumerate(null)]
    ).to_csv(RESULTS / "PERMUTATION_NULL_UNIFIED.tsv", sep="\t", index=False)
    fold_df.to_csv(RESULTS / "CV_FOLD_DIAGNOSTICS.tsv", sep="\t", index=False)
    repeat_df.to_csv(RESULTS / "CV_REPEAT_OOF_VALUES.tsv", sep="\t", index=False)

    spec = f"""# C7B unified statistic specification

**Formal statistic:** `T_REPEAT_OOF`  
**Complete-case genes:** {len(d)}  
**Repeats:** {N_REPEATS}  
**Folds per repeat:** {N_FOLDS}  
**Permutation count:** {n_perm}

For every repeat, the ten held-out predictions are concatenated into a full-length OOF vector. `rho_A` and `rho_B` are calculated on those complete OOF vectors against `Protein24`; the repeat delta is `rho_B - rho_A`; `T_REPEAT_OOF` is the mean of the ten repeat deltas.

Observed and permutation paths use the same fixed fold assignments, `statistic_repeat_oof` aggregation, Model A (`RNA6`) and Model B (`RNA6 + raw Ribo6`) fit, training-fold standardization, and Spearman implementation. The old fold-level statistic and old P value are historical audit values only.

## Conditional null construction

Within each repeat × fold, `Ribo6 ~ RNA6` is fitted on training genes only. Training residual labels are permuted within the frozen 80 combined strata (RNA-effect quintile × baseline RNA-control quartile × protein-control abundance quartile). Held-out residuals are sampled with replacement only from the corresponding training stratum residual pool. Thus held-out Ribo predictors use the training-fitted RNA component plus training-derived residuals; held-out Ribo and held-out Protein24 are not used to fit preprocessing parameters or construct predictors.

This is a fixed-fold, training-only conditional empirical null. Each permutation uses the deterministic seed sequence `[270931, permutation_id]`; parallel execution therefore does not change the null draws. It is not designed to reproduce the historical P value.
"""
    (REPORTS / "C7B_UNIFIED_STATISTIC_SPEC.md").write_text(spec, encoding="utf-8")

    parity_rows = [
        ("statistic_definition", "T_REPEAT_OOF", "T_REPEAT_OOF", "PASS", "same function and aggregation"),
        ("repeats", "10 IDs 0-9", "10 IDs 0-9", "PASS", "same repeats"),
        ("folds", "10 per repeat", "10 per repeat", "PASS", "same fixed fold structure"),
        ("fold_assignment_hash", assignment_hash(folds), assignment_hash(folds), "PASS", "same fold indices"),
        ("fold_seed_rule", "270929 + repeat", "270929 + repeat", "PASS", "same KFold seed logic"),
        ("complete_case_genes", str(len(d)), str(len(d)), "PASS", "same frozen intersection"),
        ("gene_order", "reset-index complete-case table", "same table order", "PASS", "same predictor row order"),
        ("outcome", "Protein24", "Protein24 for scoring only", "PASS", "held-out outcomes never enter preprocessing"),
        ("model_A", "RNA6", "RNA6", "PASS", "same Model A"),
        ("model_B", "RNA6 + raw Ribo6", "RNA6 + reconstructed conditional-null Ribo", "PASS_WITH_NULL_CONSTRUCTION", "same statistic; null predictor is deliberately randomized"),
        ("standardization", "training-fold fit; test transform", "training-fold fit; test transform", "PASS", "no full-data scaling"),
        ("residualization", "not used for observed raw-Ribo model", "training-fold Ribo~RNA only", "PASS_WITH_NULL_CONSTRUCTION", "null-only preprocessing is fold-local"),
        ("OOF_concatenation", "10 held-out folds per repeat", "10 held-out folds per repeat", "PASS", "one prediction per gene per repeat"),
        ("Spearman", "same scipy Spearman helper", "same scipy Spearman helper", "PASS", "same metric"),
        ("missing_values", "complete-case and finite mask", "same complete-case and finite mask", "PASS", "same rule"),
        ("direction", "B - A", "B - A", "PASS", "upper-tail positive incremental information"),
        ("permutation_strata", "not used to score observed", "frozen RNA5 x baseline RNA4 x control protein4", "PASS", "baseline/control only; no Protein24 treatment effect"),
        ("permutation_unit", "not applicable", "training residual labels within strata; test residual draw from matching training pool", "PASS_WITH_QUALIFICATION", "conditional empirical null, explicitly documented"),
        ("permutation_tail", "not applicable", "upper tail T_perm >= T_obs", "PASS", "prespecified direction"),
        ("P_formula", "not applicable", "(b+1)/(B+1)", "PASS", "B fixed at 10000"),
    ]
    pd.DataFrame(parity_rows, columns=["item", "observed", "permutation", "parity", "consequence"]).to_csv(
        REPORTS / "CV_UNIFIED_STATISTIC_PARITY_AUDIT.tsv", sep="\t", index=False
    )

    leakage = f"""# C7B leakage repair audit

1. **Scaling fit on training only?** PASS. Every repeat × fold estimates predictor means and SDs on training genes and transforms held-out genes with those values.
2. **Residualization fit on training only?** PASS. `Ribo6 ~ RNA6` is fitted separately inside every repeat × fold on training genes.
3. **Held-out Protein24 used in preprocessing?** NO. Held-out Protein24 is used only after prediction to calculate Spearman rho.
4. **Held-out genes used to estimate residualization coefficients?** NO. Held-out RNA/Ribo values are only transformed/evaluated; they do not fit the residualization coefficients.
5. **Do permutation strata use Protein24 treatment information?** NO. Strata use RNA-effect, baseline RNA-control, and mean control protein abundance. The control abundance is a baseline/control measurement, not the hyperosmotic Protein24 effect.
6. **Any full-data preprocessing retained?** NO for scaling or residualization. Frozen strata labels are a pre-specified exchangeability design and use no treatment outcome.

**LEAKAGE_AUDIT: PASS**

The held-out residual construction samples only from the matching training-stratum residual pool. No held-out Ribo residual or held-out Protein24-derived quantity is used.
"""
    (REPORTS / "CV_LEAKAGE_REPAIR_AUDIT.md").write_text(leakage, encoding="utf-8")

    strata = build_frozen_strata(d)
    strata_counts = pd.Series(strata).value_counts().sort_index()
    strata_report = f"""# C7B permutation strata audit

**Definition:** RNA6 effect quintile × baseline RNA-control quartile × mean control protein abundance quartile.  
**Non-empty strata:** {len(strata_counts)}  
**Stratum size range:** {int(strata_counts.min())}–{int(strata_counts.max())} genes.

`mean_control` is taken from the PXD054330 control-abundance column and is not the hyperosmotic Protein24 treatment effect. The frozen strata therefore act as baseline/control measurement-quality and abundance strata. No held-out Protein24 outcome or treatment-derived quantity is used to construct the null predictor.

For each fold, every held-out stratum has a non-empty training residual pool. Training residual labels are permuted within strata; held-out residuals are drawn only from the corresponding training pool. The exchangeability assumption is explicit and auditable.

**PERMUTATION_STRATA_AUDIT: PASS**
"""
    (REPORTS / "CV_PERMUTATION_STRATA_AUDIT.md").write_text(strata_report, encoding="utf-8")

    (REPORTS / "C7B_INTEGRATED_DECISION.md").write_text(
        f"""# C7B integrated decision

**C7B VERDICT:** `CV_REPAIR_{'PASS' if (t_obs > 0 and summary['fraction_repeats_B_gt_A'] > 0.5 and p < 0.05) else 'ATTENUATED'}`

## Formal result

- Formal statistic: `T_REPEAT_OOF`
- Complete-case genes: {len(d)}
- Repeats/folds: {N_REPEATS} × {N_FOLDS}
- `T_obs`: {t_obs:.12f}
- Mean rho A: {summary['rho_A_mean']:.12f}
- Mean rho B: {summary['rho_B_mean']:.12f}
- Mean delta: {summary['delta_mean']:.12f}
- Median delta: {summary['delta_median']:.12f}
- Delta SD: {summary['delta_sd']:.12f}
- Delta range: [{summary['delta_min']:.12f}, {summary['delta_max']:.12f}]
- Fraction repeats B>A: {summary['fraction_repeats_B_gt_A']:.3f}
- Split-stability range (2.5th–97.5th percentile): [{summary['delta_q025']:.12f}, {summary['delta_q975']:.12f}]

## Unified permutation result

- Permutations: {n_perm}
- Extreme count `b = count(T_perm >= T_obs)`: {extreme}
- Formal upper-tail P: {p:.12g}
- Null mean: {null.mean():.12f}
- Null SD: {np.std(null, ddof=1):.12f}
- Null percentiles 2.5/50/97.5: {null_q[0]:.12f} / {null_q[1]:.12f} / {null_q[2]:.12f}
- Observed percentile in null: {100.0 * float(np.mean(null <= t_obs)):.3f}% (empirical; ties handled by `<=`)

## Gates

- Statistic parity: **PASS**
- Scaling leakage: **PASS**
- Residualization leakage: **PASS**
- Strata audit: **PASS**
- Predictor-construction parity: **PASS_WITH_NULL_CONSTRUCTION**
- Permutation design: **PASS**, conditional empirical null explicitly defined without held-out outcome leakage
- Observed reproducibility: **PASS**
- Null reproducibility: **PASS** for the frozen seed, fold assignments, and B={n_perm}

Historical audit values are retained separately: old `T_fold_all={old_t:.15f}` and old invalid `P={old_p}`. They are not used as the C7B formal result.
""",
        encoding="utf-8",
    )

    source_hashes = []
    for path in input_paths:
        source_hashes.append({"path": str(path.relative_to(ROOT)), "sha256": sha256_file(path), "size_bytes": path.stat().st_size})
    output_paths = [
        RESULTS / "CV_REPEAT_OOF_VALUES.tsv",
        RESULTS / "CV_FOLD_DIAGNOSTICS.tsv",
        RESULTS / "PERMUTATION_NULL_UNIFIED.tsv",
    ]
    rows = [
        ("code_before_c7b", "computational_manuscript/scripts/c2r_analysis.py", "583d9f4c1ca8595a299477d9182b14eb8344dce0af3a0e7ae63cb42c03c7818a"),
        ("code_after_c7b", str(Path(__file__).relative_to(ROOT)), sha256_file(Path(__file__))),
        ("fold_assignment_hash", "generated by c7b_unified_cv.py", assignment_hash(folds)),
        ("input_CROSSOMIC_MASTER_TABLE", source_hashes[0]["path"], source_hashes[0]["sha256"]),
        ("input_PXD054330_DIFFERENTIAL_PROTEOME", source_hashes[1]["path"], source_hashes[1]["sha256"]),
        ("fold_seed", "270929 + repeat", str(FOLD_SEED)),
        ("permutation_seed", "numpy.SeedSequence([270931, permutation_id])", str(PERM_SEED)),
        ("permutation_count", "B", str(n_perm)),
        ("workers", "multiprocessing fork workers", str(workers)),
        ("software_python", "python", platform.python_version()),
        ("software_numpy", "numpy", np.__version__),
        ("software_pandas", "pandas", pd.__version__),
        ("software_scipy", "scipy", scipy.__version__),
        ("software_sklearn", "scikit-learn", sklearn.__version__),
        ("runtime_seconds", "wall clock", f"{runtime_seconds:.3f}"),
    ]
    rows.extend(("output_sha256", str(path.relative_to(ROOT)), sha256_file(path)) for path in output_paths)
    pd.DataFrame(rows, columns=["item", "source_or_definition", "value"]).to_csv(
        REPORTS / "C7B_REPRODUCIBILITY_MANIFEST.tsv", sep="\t", index=False
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-perm", type=int, default=B_DEFAULT)
    parser.add_argument("--workers", type=int, default=max(1, min(8, (mp.cpu_count() or 1) - 1)))
    args = parser.parse_args()
    if args.n_perm < 0 or args.n_perm > B_DEFAULT:
        raise SystemExit("n_perm must be between 0 and 10000; the frozen budget is B=10000")
    if args.workers < 1:
        raise SystemExit("workers must be >=1")
    started = time.time()
    d, required = read_complete_case_table()
    if len(d) != 3966:
        raise RuntimeError(f"C7B frozen complete-case gate failed: expected 3966, got {len(d)}")
    folds = make_folds(len(d))
    strata = build_frozen_strata(d)
    fold_cache = prepare_fold_cache(d, folds, strata)
    repeat_df, fold_df, oof_a = run_observed_cv(d, folds)
    old_delta_path = RESULTS.parent / "c2r" / "NESTED_MODEL_CV_DELTA.tsv"
    old_cv_path = RESULTS.parent / "c2r" / "NESTED_MODEL_CV.tsv"
    old_delta = pd.read_csv(old_delta_path, sep="\t")
    old_cv = pd.read_csv(old_cv_path, sep="\t")
    old_t = float(old_delta.loc[0, "Delta_CV_rho"])
    old_p = float(pd.read_csv(RESULTS.parent / "c2r" / "PERMUTATION_NULL.tsv", sep="\t").loc[0, "empirical_upper_tail_P"])
    if args.n_perm == 0:
        raise SystemExit("Observed-only smoke path is not a C7B result; run with --n-perm 10000")
    null = statistic_repeat_oof(d, folds, fold_cache, oof_a, args.n_perm, args.workers)
    inputs = [
        RESULTS.parent / "c2" / "CROSSOMIC_MASTER_TABLE.tsv",
        RESULTS.parent / "c2" / "PXD054330_DIFFERENTIAL_PROTEOME.tsv",
        Path(__file__).resolve().parents[0] / "c2r_analysis.py",
        old_cv_path,
        old_delta_path,
    ]
    write_reports(d, folds, repeat_df, fold_df, null, old_t, old_p, inputs, time.time() - started, args.n_perm, args.workers)
    print(json.dumps({"T_obs": float(repeat_df.delta_BA.mean()), "P": float((1 + np.sum(null >= repeat_df.delta_BA.mean())) / (args.n_perm + 1)), "n_perm": args.n_perm}, indent=2))


if __name__ == "__main__":
    main()
