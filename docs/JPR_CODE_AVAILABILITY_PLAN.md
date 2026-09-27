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

The release should record Python 3.9.6, R 4.4.3, NumPy 1.26.4, pandas 2.3.3, SciPy 1.13.1, statsmodels 0.14.6, scikit-learn 1.6.1, matplotlib 3.9.4, seaborn 0.13.2, and openpyxl 3.1.5, subject to final environment verification. Seeds, thresholds, model formulas, fold construction, and resampling budgets are already recorded in the scripts/configuration. The workspace has no Git repository or commit identifier; a new version-controlled public repository must be created before submission.

## Data-use boundary

Large third-party raw files are not copied into the code repository. The accession list and the processed matrix provenance are retained. The derived outputs are nonrestricted unless repository terms or contributor agreements indicate otherwise.

## Required author action

Create the public repository, add a release tag matching the submitted manuscript, deposit an archive in Zenodo or Figshare, update the DOI and URL in the manuscript/availability statements, and run the scripts from a clean environment before submission.
