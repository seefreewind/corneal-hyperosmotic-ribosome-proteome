# Simulated Reviewer 1 — Proteomics / Multi-omics

1. **Processed-matrix normalization:** Please provide the exact upstream normalization and Spectronaut settings for PXD054330. **Status:** UNRESOLVABLE_LIMITATION; the manuscript now reports this explicitly and supplies matrix checksums/QC.
2. **Protein mapping:** Clarify why multi-gene groups were excluded and whether aggregation preceded the contrast. **Status:** ALREADY_ADDRESSED; Methods and SI specify unique 1:1 retention and gene-level means.
3. **Proteome sample size:** Three runs per condition limit protein-effect precision. **Status:** ALREADY_ADDRESSED; n=3 is stated and limitation is explicit.
4. **Cross-study mismatch:** Explain why 500 mOsm/6 h can be compared with 450 mOsm/24 h. **Status:** ALREADY_ADDRESSED; the claim is cross-study continuity, not longitudinal transfer.
5. **Ribo versus RNA estimand:** Define Ribo6 separately from TE6. **Status:** ALREADY_ADDRESSED.
6. **Residualization:** Explain why OLS residualization does not establish causal independence. **Status:** TEXT_FIX completed; residual is called RNA-adjusted information, not causal TE.
7. **Partial rank association:** Provide the exact rank-residual procedure and reciprocal result. **Status:** ALREADY_ADDRESSED.
8. **Cross-validation leakage:** Confirm fold construction and fold-specific scaling. **Status:** ALREADY_ADDRESSED; Methods and preregistration specify shared folds and training-fold scaling.
9. **Permutation design:** Justify the matched strata and fixed folds. **Status:** ALREADY_ADDRESSED in Methods/SI.
10. **Gene dependence:** Explain chromosome-block resampling and its scope. **Status:** ALREADY_ADDRESSED; 23 blocks and limitation stated.
11. **Detectability bias:** Clarify that PXD059451 has no abundance outcome. **Status:** ALREADY_ADDRESSED.
12. **Public-data novelty:** Distinguish the present test from the source papers. **Status:** TEXT_FIX completed in v2 Introduction/Cover Letter.
13. **Title “anticipates”:** Could imply prospective prediction. **Status:** RESOLVED by selecting conservative T2.
14. **Figure scaling:** Do not visually inflate ρ=0.11–0.16. **Status:** FIGURE_FIX; full scatter axes and explicit values retained.

No new analysis is required to address these comments.
