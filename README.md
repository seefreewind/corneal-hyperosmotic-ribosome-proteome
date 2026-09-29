# Cross-study translational information in corneal epithelial hyperosmolarity

This repository package contains the frozen computational workflow supporting the manuscript **“Ribosome-Level Responses Provide Incremental Information about Later Proteomic Remodeling under Hyperosmotic Stress in Human Corneal Epithelial Models.”** It tests whether early ribosome-level changes contain incremental information about an independent later proteomic response beyond RNA-level changes.

## Study overview

The workflow integrates GSE200097 (early RNA-seq and ribosome profiling), GSE323164 (independent 24 h transcriptome), PXD054330/JPST003233 (independent 24 h DIA proteome), and PXD059451 (baseline HCEC detectability sensitivity). The analysis is cross-study and non-longitudinal. It does not establish causality, PIEZO1 dependence, mechanotransduction, or clinical prediction.

## Public inputs

- GSE200097 — early RNA/Ribo responses.
- GSE323164 — independent HCE-2 transcriptome.
- PXD054330 / JPST003233 — HCEC DIA proteome.
- PXD059451 — HCEC detectability resource.

Large third-party raw files are not included. Accession identifiers and download instructions are recorded in `docs/DATA_ACCESSION_AUDIT.tsv`. The processed PXD054330 pivot and source tables may only be redistributed after author and repository-term review; see `DATA_USE_NOTICE.md`.

## Required software

The locked C4 clean run used Python 3.11.15 on macOS with packages in `environment/requirements.txt`. R 4.4.3 is recorded for the broader project environment; the C2/C2R manuscript-critical scripts are Python-based. Recreate the environment with `python -m venv .venv` and `pip install -r environment/requirements.txt`, or use `environment/environment.yml` with conda.

## Exact analysis order

1. Recover or place the permitted processed inputs according to the project path map.
2. Run the C1M matrix and identifier audit using the frozen mapping rules.
3. Run `scripts/c2_analysis.py`.
4. Run `scripts/c2r_analysis.py`.
5. Compare outputs with `docs/EXPECTED_KEY_OUTPUTS.tsv` and the frozen preregistrations.

The C4 rerun used frozen processed inputs for C1M because network refresh of UniProt/MyGene mappings is not deterministic. The release includes the mapping provenance and input checksums; `docs/c1m_frozen_local_audit.py` records the frozen local audit logic.

## Reproduction commands

From the repository root after restoring the permitted input paths:

```bash
python scripts/c2_analysis.py
python scripts/c2r_analysis.py
```

Run the C1M audit before these commands. No raw DIA/RNA third-party files should be copied into this repository without permission.

## Corrected v1.0.2 outputs

The primary set is 3,974 genes, with 3,966 complete cases. The corrected C7B release uses ten repeated 10-fold gene-level splits, complete out-of-fold aggregation within each repeat, and the formal statistic `T_REPEAT_OOF = mean(Δρr)`. Model A mean ρ is 0.10969, Model B mean ρ is 0.15012, and mean Δρ is 0.04042; all ten repeat-level increments are positive. The matched conditional residual-permutation null uses B = 10,000 and b = 0, giving formal P = 0.00009999000099990002, displayed as P < 1 × 10⁻⁴. The corrected outputs and source tables are in `results/c7b/`, `figures/`, and `figures_source/`.

The observed and permuted statistics use identical repeat-level complete-OOF aggregation with training-only preprocessing. The correction harmonizes statistic definition only; it does not change the permutation design, datasets, genes, folds, seeds, thresholds, or primary paired-gene results. The detailed release audit is in `docs/release_v1.0.2/`.

The C4 release audit is recorded in `docs/FINAL_REPRODUCIBILITY_SNAPSHOT.md`, `docs/FINAL_RESULT_CONSISTENCY_AUDIT.tsv`, `docs/FINAL_REFERENCE_AUDIT.tsv`, `docs/FINAL_ACCESSION_AUDIT.tsv`, `docs/FIGURE_SOURCE_DATA_AUDIT.tsv`, `docs/SI_FINAL_QC.tsv`, and `docs/PUBLIC_DATA_NOVELTY_PARAGRAPH.md`.

## License and citation

Original code and documentation are released under the MIT License. Source accession data remain governed by their originating repositories; see `DATA_USE_NOTICE.md`. Cite the manuscript and the original accession records; see `CITATION.cff`.

## AI assistance disclosure

ChatGPT and DeepSeek were used to assist with code writing, language polishing, and formatting adjustments. They were not authors. The authors reviewed and validated the analyses, numerical results, interpretations, figures, source-data links, and references. No generative-AI manuscript figure or TOC imagery was used.

## Release status

The public GitHub repository is https://github.com/seefreewind/corneal-hyperosmotic-ribosome-proteome. Release `v1.0.2` contains the corrected C7B outputs and is archived in Zenodo at https://doi.org/10.5281/zenodo.23048255; the all-versions DOI remains https://doi.org/10.5281/zenodo.22985907. The historical `v1.0.1` record is https://doi.org/10.5281/zenodo.22985978. No manuscript submission has been made.
