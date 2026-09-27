# Final Reproducibility Snapshot

**Date:** 2026-09-27  
**Clean environment:** PASS; Python 3.11.15 venv with pinned requirements; original project environment recorded separately.  
**C1M:** PASS using frozen processed inputs and frozen mapping; 4,301 mapped rows, 3,974 primary overlap, 3,950 four-way overlap, 4,518 complete protein rows.  
**C2/C2R:** PASS; all manuscript-critical targets reproduced within tolerance.  
**Comparison:** `reproducibility/c4_clean_run/C4_REPRODUCIBILITY_COMPARISON.tsv`.  
**Input checksums:** `reproducibility/C4_INPUT_SHA256.txt`.  
**Environment capture:** `reproducibility/C4_ENVIRONMENT_CAPTURE.txt`.  
**Random seeds:** `config/FINAL_RANDOM_SEEDS.yaml`.  
**Git commit:** `ef469a5` (frozen public-repository package commit); no remote publication was performed. Subsequent documentation-only metadata may be committed separately.

## Key outputs

- RNA6→Protein24 ρ=0.11347795924096821.
- Ribo6→Protein24 ρ=0.16326415613882245.
- Δρ=0.04971168011039095.
- RNA-adjusted Ribo residual ρ=0.11058808019012686.
- Partial Ribo|RNA ρ=0.11845983274673053.
- Model B−A CV Δρ=0.039280648992946095.
- Structured permutation P=0.000099990001.
- Chromosome-block bootstrap Δρ=0.049833798404433.

No exploratory analysis, threshold change, new pathway search, PIEZO1 analysis, or fate-class redefinition was added.
