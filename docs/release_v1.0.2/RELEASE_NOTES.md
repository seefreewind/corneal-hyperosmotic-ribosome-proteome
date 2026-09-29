# v1.0.2 release notes

This release adds the corrected C7B repeat-level complete-out-of-fold analysis and its reproducibility package.

## Correction scope

During the pre-submission statistical audit, the earlier cross-validation/permutation implementation used different aggregation rules for the observed and permuted statistics. The v1.0.2 package uses a unified repeat-level complete-OOF statistic, with training-only preprocessing applied identically to observed and permuted data. The permutation design, data sources, genes, fold assignments, seeds, thresholds, and primary paired-gene results are unchanged.

## Canonical outputs

- Model A mean complete-OOF Spearman ρ: 0.10969.
- Model B mean complete-OOF Spearman ρ: 0.15012.
- Mean repeat-level Δρ: 0.04042; 10/10 repeats positive.
- 95% split-stability interval: 0.03712–0.04408.
- Conditional residual-permutation null: B = 10,000, b = 0, formal P = 0.00009999000099990002, displayed as P < 1 × 10⁻⁴.

Raw third-party files are not redistributed. Accession identifiers in `docs/DATA_ACCESSION_AUDIT.tsv` provide the source-data route.

## Zenodo

The v1.0.2 Zenodo DOI will be inserted after the deposit is created and independently verified.
