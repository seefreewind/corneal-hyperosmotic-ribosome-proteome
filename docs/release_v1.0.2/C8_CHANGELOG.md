# C8 changelog

C8 is a pre-release packaging and consistency gate after the C7C scientific freeze.

- No new scientific analysis, rerun, resampling, fold assignment, permutation design, seed, threshold, gene set, pathway, or claim was introduced.
- Main manuscript, Table 2, SI, source-data notes, and cover letter were normalized to the C7B canonical values and JPR-friendly P typography.
- Mixed PMID/PMCID fields were removed from the C8 candidate reference list; DOI punctuation was normalized.
- Figure 3 uses the C7B final asset and source data; its displayed permutation result is `P < 1 × 10⁻⁴`, with formal `P = 0.00009999000099990002` retained in Methods/SI.
- Table 2 now points to `TABLE2_SOURCE_C8.tsv`; its CV row is labeled as a split-stability interval.
- DOCX metadata cleanup removed comments/people parts and tracked-change tags from generated C8 DOCX files.
- Data/code wording distinguishes the historical v1.0.1 package from the corrected public v1.0.2 release.
- Cover letter, submission checklist, author metadata, release checklist, Zenodo metadata, and TOC graphic specification were prepared.

Transparent correction wording used in the cover letter: “During pre-submission statistical audit, an earlier CV/permutation implementation was found to use different aggregation rules for the observed and permuted statistics. Before submission, this branch was replaced with a unified repeat-level complete out-of-fold statistic with training-only preprocessing applied identically to observed and permuted data.”

No GitHub release, Zenodo upload, or journal submission was performed.
