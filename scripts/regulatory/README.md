# Molecular analyses

Run from the root of the distributed project with the released supplementary workbook.

```sh
python scripts/regulatory/recalculate.py --archive archive --tables data/SAMD9L_Supplementary_Tables.xlsx --output results/regulatory
```

The script calculates the molecule-level statistics for Supplementary Tables S4–S5 and their figure panels. The run builds the patient–event roster from S1A and S2B. The outputs include the individuals and phase evidence contributing to each analysis. The original per-track exports came from AlphaGenome 1 using client 0.9.0. All computations below are offline. New sequence predictions require an external model run.

## Inputs and units

- **Isolated alleles.** The primary fit uses 88 primary germline and 50 acquired origin-labelled entries in each of nine modalities. Normalized inputs, query metadata and the archived track exports cover all origin-labelled entries. Predictions used 1-Mb GRCh38 windows. The primary regression fits the signed median hematopoietic quantile as a function of origin, scaled linear position and truncating consequence with HC3 uncertainty, followed by BH correction across nine modalities. Two DNA alleles carry both origin labels. Four germline or context modifier entries lie outside the primary allele set. Supplementary Table S4J reports alternative specifications starting from all 88 germline and 50 acquired entries. The run recalculates those 27 estimates from the corresponding rows of the final score matrix and compares them with S4J. The germline-only position analysis uses six splice-junction comparisons with a separate adjustment.
- **Paired predictions.** S2B has 57 patient–event observations. SGR016 and SGR025 were demonstrated in trans. SGR044 has two germline variants. Of the remaining 54 with mappable residues, all 54 have archived combined-sequence predictions. The nested phase sets contain 12, 23, 25 and 54 observations. One event (SGR007) has a cis configuration among several clones and was also observed on the reference background at the germline site. The observed cis pairing contributes once. The files with prefix `omit_dual_configuration_` give a sensitivity excluding this observation, with original subset labels retained for cross-reference. The subset sizes in these files are one lower. Per-modality effective counts appear in each result.
- **Gene-centered predictions.** The program uses native and native hematopoietic strong-track summaries, with Wilcoxon tests across germline backgrounds, followed by adjustment across four modalities in each scope, subset and endpoint.
- **Fixed-window RNA.** Five frozen alternatives per scored observation yield 270 assignments and 175 distinct alternative sequences. Common finite native RNA tracks were matched in the original exports. The analysis takes successive medians over tracks, assigned alternatives and observed events for each germline background. Exact two-sided sign tests omit background differences of zero and use Benjamini–Hochberg adjustment across eight phase-by-scope comparisons.
- **Paired chromatin.** The analysis preserves the archived track thresholds defined from observed germline effects. For each of the two 48-comparison families, 10,000 draws sample alternatives within each event. The program returns all fitted parameters, endpoint-specific event counts, available control identities, adjusted probabilities and the simulated distributions.
- **Protein geometry.** All 54 eligible pairs contribute, using their mapped protein residues. The archived AlphaFold v6 monomer and PAE matrix supply all distances. The 100,000 permutations use the original generator, seed 20261004 and 5,000-row batches. The program verifies the recomputed geometry of pairs with previously tabulated distances against the original individual distance table and reports results with and without SGR056.

## Rebuilding the archived track summaries

Use the complete review archive to rebuild allele and configuration summaries from the track-level records. The compact archive starts from the frozen summaries supplied under `archive/analysis/`.
