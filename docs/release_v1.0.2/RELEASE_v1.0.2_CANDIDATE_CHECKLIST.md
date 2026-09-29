# v1.0.2 release-candidate checklist

Status: **CANDIDATE PREPARED — NOT PUBLISHED**

- [x] Scientific freeze remains C7C/C7B; no analysis rerun in C8.
- [x] Corrected repeat-level complete-OOF statistic and matched permutation outputs identified.
- [x] Final Figure 3, Table 2, source-data, SI, and C8 audit artifacts prepared.
- [x] Release notes include the observed/permutation aggregation correction and its scope.
- [ ] Copy candidate artifacts into the public repository working tree.
- [ ] Run repository tests and checksum manifest.
- [ ] Update README, CITATION.cff, and data/code availability metadata.
- [ ] Commit and tag `v1.0.2` only after author approval.
- [ ] Push GitHub release.
- [ ] Verify GitHub release contents and DOI metadata after publication.

Proposed candidate payload:

- `scripts/c7b_unified_cv.py`
- `scripts/c7c/build_figure3_c7b.py`
- `results/c7b/CV_REPEAT_OOF_VALUES.tsv`
- `results/c7b/PERMUTATION_NULL_UNIFIED.tsv`
- `figures/FIGURE3_C7B_FINAL.*`
- `figures/source_data/final/FIGURE3_SOURCE_C8.tsv`
- `figures/source_data/final/FIGURE3_PERMUTATION_NULL_C7B.tsv`
- `figures/source_data/final/TABLE2_SOURCE_C8.tsv`
- `reports/final/c7b/`
- `reports/final/c7c/`
- `reports/final/c8/`
- `docs/JPR_DATA_AVAILABILITY.md`

No push or release operation was performed by C8.
