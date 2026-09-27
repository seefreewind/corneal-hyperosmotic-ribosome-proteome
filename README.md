# Cross-study translational information in corneal epithelial hyperosmolarity

This repository package contains the frozen computational workflow supporting the Journal of Proteome Research manuscript **“Ribosome-Level Responses Provide Incremental Information about Later Proteomic Remodeling during Corneal Epithelial Hyperosmotic Stress.”** It tests whether early ribosome-level changes contain incremental information about an independent later proteomic response beyond RNA-level changes.

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

## Expected key outputs

The primary set is 3,974 genes, with 3,966 complete cases. Expected C2/C2R values are in `docs/EXPECTED_KEY_OUTPUTS.tsv`; the frozen random seeds are in `config/FINAL_RANDOM_SEEDS.yaml`.

The C4 release audit is recorded in `docs/FINAL_REPRODUCIBILITY_SNAPSHOT.md`, `docs/FINAL_RESULT_CONSISTENCY_AUDIT.tsv`, `docs/FINAL_REFERENCE_AUDIT.tsv`, `docs/FINAL_ACCESSION_AUDIT.tsv`, `docs/FIGURE_SOURCE_DATA_AUDIT.tsv`, `docs/SI_FINAL_QC.tsv`, and `docs/PUBLIC_DATA_NOVELTY_PARAGRAPH.md`.

## License and citation

Original code and documentation are released under the MIT License. Source accession data remain governed by their originating repositories; see `DATA_USE_NOTICE.md`. Cite the manuscript and the original accession records; see `CITATION.cff`.

## AI assistance disclosure

ChatGPT and DeepSeek were used to assist with code writing, language polishing, and formatting adjustments. They were not authors. The authors reviewed and validated the analyses, numerical results, interpretations, figures, source-data links, and references. No generative-AI manuscript figure or TOC imagery was used.

## Release status

The public GitHub repository is https://github.com/seefreewind/corneal-hyperosmotic-ribosome-proteome. Release `v1.0.0` is archived in Zenodo at https://doi.org/10.5281/zenodo.22985908. No manuscript submission has been made.
