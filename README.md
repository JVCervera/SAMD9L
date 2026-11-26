# Germline variation and somatic adaptation in SAMD9L-associated disease

Data, analysis code and figures for the *Leukemia* Letter by Alessandro Liquori and José Cervera (revision v20).

The GitHub archive contains the reviewed workbook, code, artwork and archived model summaries. The complete archive accompanying the *Leukemia* submission additionally contains the full track exports.

## Reproduce the analyses

Use 64-bit Python 3.12. From the extracted directory, create a virtual environment and install the pinned dependencies.

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python reproduce.py --skip-track-checks
```

On Windows, activate with `.venv\Scripts\activate`. All analyses run locally after dependency installation. The runner writes results, logs and `verification.json` to `reproduced/`. A successful run reports `"status": "PASS"` after it checks reference statistics, participant membership, workbook cells and revised figure assets. Use `--output PATH` to select a new output directory or `--workbook PATH` to use a separate copy of the submission workbook.

The compact run starts from the frozen exported summaries and checks the resulting statistics.

`python reproduce.py --copy-figures-only --output submission_artwork` copies the exact editorial PDFs, SVGs and main TIFFs from `reference/figures/` into `submission_artwork/figures/` without running analyses. The same assets accompany the default analysis output. To make fresh PDFs from the SVGs, install Inkscape and add `--render-svg`. Those PDFs are written to `figures_rendered/`. Arial is the first choice for the SVG text, with Liberation Sans as a fallback.

`--redraw-figures` generates independent analytical redraws in `figures_redrawn/`. Layout, fonts and styling may differ from the editorial figures in `reference/figures/`. The Figure 1B input in `data/Figure_1B_points.csv` contains one clinical-context age and germline residue per eligible individual. The default run checks each point against the recomputed clinical roster. Figure S2 uses the broader reported-age endpoint.

## Analysis inputs and generated outputs

| Analysis | Inputs | Outputs relative to the run directory |
| --- | --- | --- |
| Clinical phenotypes, variants, clonal observations and age–position contrasts | Supplementary Tables S1–S3 and the source and family support checks in `data/` | `clinical/clinical_summary.json`, `age_rows.csv`, `age_exclusions.csv`, `family.csv`, source strata, inflammation, gnomAD and clonal tables |
| Any-time chromosome 7 loss | `data/pooled_roc_input.csv`, linked to S3I | `clinical/roc/` contains the case-level input, held-out predictions, ROC coordinates, bootstrap draws and summary |
| Consecutive and selected pediatric MDS comparisons | Source-specific baseline karyotypes in S3 | `clinical/sahoo_baseline.csv` and four exact comparisons in the clinical summary |
| Regulatory predictions and alternate configurations | AlphaGenome exports and frozen control assignments in `archive/analysis/` and event assignments from S2 | `regulatory/` contains adjusted origin analyses, junction results, paired gene-centered and RNA results and matched chromatin null distributions |
| Structural comparisons | AlphaFold DB reference monomer and predicted aligned errors, plus mapped residues in S2 | `regulatory/structural_pairs.csv`, geometric summaries and reassociation results |
| Figures | Adjudicated workbook, Figure 1B point data and reference artwork | `figures/` contains the exact submission assets. Optional numerical redraws appear in `figures_redrawn/`. |

The active clinical workbook is `data/SAMD9L_Supplementary_Tables.xlsx`. The revision uses 178 individuals, 169 with documented manifestations. Clinical, first-event, hematologic-event and broader reported ages have distinct eligibility rules. The source-coded diagnosis comparison and exploratory broad phenotype comparison use the same 126 eligible source records, with separate definitions for each outcome. A missing diagnosis entry is treated as absence of documentation in that source field.

The any-time analysis uses the reviewed series of 105 individuals with documented chromosome 7 data and one resolvable germline allele position per individual. It includes 57 losses and 48 comparators. The primary Firth and held-out models use all 105 individuals with documented age groups, including P102, documented as pediatric (PMID 34417303). P102 enters those models by age group. Its clinical-context age is recorded separately for Figure 1B. S3J gives adjusted sensitivities. S3M reports the 105-person validation and omission of substitutions at residue 986. The runner saves the standalone Supplementary Figure S5B panel and the cutoff data under `clinical/letter_roc_105/`. The roster records the age and evidence sources and the observed status. The historical protein labels in `resampling_group` preserve the order of the same bootstrap groups across nomenclature corrections. Evidence for P064 is limited to a 46,XX lymphocyte karyotype, which is documented in the comparator roster. P050 contributes its observed marrow monosomy 7 despite having normal blood counts and no documented disease manifestations at that evaluation. The pooled analysis and the pediatric baseline comparison address distinct clinical questions.

## Molecular inputs

The offline analysis starts from frozen AlphaGenome 1 track exports generated with client 0.9.0. It reruns the statistical comparisons from the archived summaries. A new model-side generation of predictions requires external service access. Source references, version identifiers and hashes for the archived AlphaFold model and predicted aligned error matrix are in `archive/analysis/model_inputs/manifest.json`. Track provenance is listed in `archive/analysis/track_inputs/sources.csv`. See [Molecular analyses](scripts/regulatory/README.md) for the exact modules and analytical units.

The acquired variant inventory comprises 50 distinct nucleotide-defined alleles across 57 patient–event observations in 38 individuals. Isolated predictions are available for all 50 acquired alleles. The combined-sequence analyses use the 54 configurations with paired predictions. The structural analysis includes 54 configurations. Reviewed protein annotations and residue anchors are listed in `data/protein_normalizations.csv`. The archived model exports retain their source labels and nucleotide-level inputs.

## Release integrity

Use `sha256sum -c SHA256SUMS` in the repository root to verify the distributed files. The full validation record is in `validation/verification.json`. See [sources and attribution](RIGHTS.md) and [citation metadata](CITATION.cff).

The included Colab notebook documents the paired chromatin prediction run. The fixed-window RNA and chromatin designs have 54 observed configurations and five assigned alternatives per configuration (175 distinct alternatives). The input configuration and control manifests in `archive/analysis/portable_inputs/`, the frozen track exports, and the runnable Python scripts document and reproduce the reported analyses offline. Generating new predictions uses the external model service. The final chromatin statistics use the completed 229-configuration export.
