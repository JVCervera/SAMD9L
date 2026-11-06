# SAMD9L germline context and somatic adaptation

Reproducibility code for **Germline variation and somatic adaptation in SAMD9L-associated disease**, Alessandro Liquori and José Cervera.

The repository contains curated analytical inputs, saved prediction summaries, statistical algorithms, model coordinates, reference results and publication figure code. The cohort comprises 178 individuals (170 affected), 88 primary germline alleles and 48 exact acquired alleles. Combined-sequence and structural analyses use 53 eligible configurations: seven confirmed cis, 11 additional source-reported cis, two inferred cis and 33 unresolved.

## Quick start

Use Python 3.12. Create a virtual environment and install the dependencies:

    python -m venv .venv
    # Linux/macOS:
    source .venv/bin/activate
    # Windows PowerShell instead:
    # .venv\Scripts\Activate.ps1
    python -m pip install -r requirements.txt
    python reproduce.py

Results and a machine-readable verification.json are written to results/. The script runs from any working directory and does not modify included inputs or reference results.

    # Also generate Figures 1–4 and S1–S2 in PDF, PNG, SVG and HTML:
    python reproduce.py --figures

    # Choose an output directory:
    python reproduce.py --output /path/to/results

    # Verify the 96 saved chromatin null distributions without resampling:
    python reproduce.py --verify-only

The full offline workflow includes the 10,000-draw chromatin comparisons, 100,000 structural reassociations and 10,000-replicate documented-family bootstrap. It makes **no AlphaGenome requests** and needs no API key. Reference outputs are compared automatically.

For publication typography, install Nimbus Sans (fonts-urw-base35 on Debian/Ubuntu). Without this font, Matplotlib uses its installed fallback; numerical data and layouts are retained but typography may differ. PDF/PNG exports do not require a browser.

## What is reproduced

| Manuscript component | Entry point | Output / supplementary table |
| --- | --- | --- |
| Age, inflammation, source strata, chromosome 7 and population comparisons | analysis/clinical_statistics/clinical_and_regulatory.py | clinical/; S2–S3 |
| Documented-family and source-support sensitivities | analysis/run_family.py | family/; S3F |
| Isolated-allele origin models and sensitivities | analysis/clinical_statistics/clinical_and_regulatory.py | clinical/origin_regression_sensitivity.csv; S4B,J |
| Six position/splice-junction tests | analysis/run_position_splicing.py | splicing/; S4C |
| Gene-centered tests by phase and scope | analysis/clinical_statistics/clinical_and_regulatory.py | clinical/gene_centered_phase_statistics.csv; S4G,I |
| Fixed-window RNA matched comparisons | analysis/run_match_summary.py | rna/; S4D,E |
| Chromatin PRIMARY and EFFECTIVE_MIN3 families | analysis/run_chromatin.py | chromatin/; S4F (48 + 48) |
| Reference-monomer geometry and reassociations | analysis/run_structure.py | structure/; S5 |
| Publication figures, including detailed lollipop | figure_code/ | figures/; Figures 1–4 and S1–S2 |

Figure generation uses the fixed publication source tables after statistical verification. It does not substitute a different analysis or infer missing patient information. reference/Supplementary_Tables.xlsx is the reader-facing S1–S5 workbook. Computational workbook sheet names such as S14_Age_position and S24_Structural_pairs are internal dataset identifiers, not supplementary table numbers.

## Files

- analysis/: statistical code and curated numerical inputs.
- analysis/portable_inputs/: fixed phase/control assignments and exported prediction metrics.
- analysis/model_inputs/: versioned AlphaFold PDB/PAE files and source manifest.
- analysis/clinical_statistics/: source records and family annotations.
- figure_code/: final publication plotting code and source tables.
- reference/: frozen results, null distributions, supplementary workbook and publication PDFs.
- validation/: result of the delivered package's offline execution.
- MANIFEST_SHA256.csv: content hashes for delivered files.
- docs/: methods, input provenance and reproducibility scope.

Missing predictions, uncertain phase and unknown clinical observations are not assigned zero values or negative findings. Related individuals remain in individual-level analyses; the documented-family sensitivity uses only resolved family/de novo units.

## Regulatory predictions and external inference

The offline starting point is the exported prediction summaries and fixed matching assignments. These reproduce the statistical results without a new model call. This repository does not claim to recreate every original per-track prediction from DNA sequence.

An optional fixed-window RNA inference entry point is supplied:

    python -m pip install -r analysis/requirements_inference.txt
    # Set ALPHAGENOME_API_KEY in your environment; do not commit it.
    python analysis/run_rna_inference.py --execute --output results/new_rna_predictions

This optional path requires AlphaGenome service access and external GRCh38 sequence retrieval, and pins the client to 0.9.0. It is not part of offline verification and has not been rerun for this release. Its inference input retains the 54-event design for alignment; downstream study statistics exclude demonstrated-trans event SGR025 and use the 53-event phase manifest. The dual-germline observation is excluded. Backend changes may change new predictions even with the same client. Do not automatically replace publication reference results with new service outputs.

## Uploading to GitHub

Create an empty repository, extract this folder and upload its contents, retaining the directory structure. Include the data inputs and reference/ files: scripts alone are insufficient. Do not upload generated results/, a virtual environment, API keys or .env files.

CITATION.cff supplies the title and authors; add the repository URL and article DOI when available. No DOI or release URL is invented here. See docs/RIGHTS.md for code and source-data rights.

These model outputs are not measurements of cellular rescue or individual prognosis.

## Final submission materials
Figure S3 is distributed as the published image in figure_code/Figure_S3_published.png. The figure generator covers Figures 1–4 and S1–S2. The age-position dataset S14 supplies the individual-level observations for S3. The junctions directory contains all eligible variant–junction–track rows for Figure 4C. Run python analysis/junctions/verify.py to reproduce the six position contrasts from these rows. SGR053 remains analytically unresolved because the authors proposed probable cis without direct or clonal phasing evidence.

## Final submission materials
Figure S3 is distributed as the published image in figure_code/Figure_S3_published.png. The figure generator covers Figures 1–4 and S1–S2. The age-position dataset S14 supplies the individual-level observations for S3. The junctions directory contains all eligible variant–junction–track rows for Figure 4C. Run python analysis/junctions/verify.py to reproduce the six position contrasts from these rows. SGR053 remains analytically unresolved because the authors proposed probable cis without direct or clonal phasing evidence.

Sahoo cohorts A and B (Supplementary Figure S4 and Table S3G): run python analysis/sahoo/compare.py. The case-level input is transcribed from Sahoo et al. Tables 3 and 6 (PMID 34621053). The script calculates four exact tests and renders the figure.
