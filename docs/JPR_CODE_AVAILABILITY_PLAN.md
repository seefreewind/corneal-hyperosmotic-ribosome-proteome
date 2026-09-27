# JPR Code Availability Plan

## Release contents

- `scripts/c2_analysis.py`: frozen C2 primary cross-omic analysis and report generation.
- `scripts/c2r_analysis.py`: frozen C2R residual, partial, nested-CV, permutation, detectability, stratification, robust-correlation, module, amino-program, and block-bootstrap analyses.
- `scripts/build_jpr_package.py`: deterministic packaging and draft-figure script.
- `config/C2_ANALYSIS_PREREGISTRATION.yaml` and `config/C2R_ROBUSTNESS_PREREGISTRATION.yaml`.
- Nonrestricted derived TSV/JSON outputs under `results/c2/` and `results/c2r/`.
- Mapping, sample-map, and matrix-QC files under `metadata/` and `results/c1m/`.
- A plain-text environment manifest generated at release time.

## Reproducibility requirements

The original frozen analysis used Python 3.9.6 and R 4.4.3. The C4 clean rerun used Python 3.11.15 on macOS with NumPy 1.26.4, pandas 2.3.3, SciPy 1.13.1, statsmodels 0.14.6, scikit-learn 1.6.1, matplotlib 3.9.4, seaborn 0.13.2, openpyxl 3.1.5, requests 2.32.3, and tabulate 0.9.0. Seeds, thresholds, model formulas, fold construction, and resampling budgets are recorded in the scripts and `config/FINAL_RANDOM_SEEDS.yaml`. The frozen repository package is publicly available at https://github.com/seefreewind/corneal-hyperosmotic-ribosome-proteome and archived at https://doi.org/10.5281/zenodo.22985908 (GitHub release `v1.0.0`).

## Data-use boundary

Large third-party raw files are not copied into the code repository. The accession list and the processed matrix provenance are retained. The derived outputs are nonrestricted unless repository terms or contributor agreements indicate otherwise.

## Required author action

The public repository, `v1.0.0` release tag, Zenodo archive, DOI, manuscript availability statements, and clean-environment rerun are complete. Before submission, retain the accession-based data-use boundary and perform the final journal-format and deposited-table checks.
