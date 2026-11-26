"""Compare the computed patient selections and table entries with the publication."""
import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from openpyxl.utils.cell import coordinate_to_tuple

ROOT = Path(__file__).resolve().parents[1]


def verify_membership(output, membership):
    people = pd.read_csv(membership)
    assert len(people) == 178 and people["Patient ID"].is_unique
    endpoints = {"Clinical age": "clinical", "First clinical event": "first",
                 "Hematologic event": "hematologic", "Broad reported age": "mixed"}
    reference = people.melt(id_vars="Patient ID", value_vars=endpoints,
                            var_name="endpoint", value_name="age").dropna(subset=["age"])
    reference["endpoint"] = reference.endpoint.map(endpoints)
    reference = reference.rename(columns={"Patient ID": "patient"}).set_index(["patient", "endpoint"]).age
    rows = pd.read_csv(output / "clinical/age_rows.csv")
    pd.testing.assert_series_equal(rows.set_index(["patient", "endpoint"]).age.sort_index(),
                                   reference.sort_index(), rtol=1e-10, atol=1e-10)
    for field, filename in [("Fixed ROC", "roc_input.csv"), ("Age-adjusted ROC", "held_out_predictions.csv")]:
        observed = pd.read_csv(output / "clinical/roc" / filename).patient_id
        expected = people.loc[people[field].eq("Yes"), "Patient ID"]
        assert observed.is_unique and set(observed) == set(expected), field
    return len(rows)


def verify_sahoo_immune_phenotypes(cohort_sheet):
    positives = pd.read_csv(ROOT / "data/sahoo_immune_positive.csv")
    assert len(positives) == 9 and positives.patient_id.is_unique and positives.sahoo_upn.is_unique
    headers = cohort_sheet[4]
    individuals = pd.DataFrame(cohort_sheet[5:183], columns=headers).set_index("Patient ID")
    assert individuals.index.is_unique
    for item in positives.itertuples(index=False):
        record = individuals.loc[item.patient_id]
        assert item.sahoo_upn in record["Aliases"].split("; "), (item.patient_id, item.sahoo_upn)
        assert "Immune/inflammatory" in record["Phenotype group"], item.patient_id
    return len(positives)


def compare_numeric_block(sheet, first_row, frame):
    count = 0
    for row_number, record in enumerate(frame.itertuples(index=False, name=None), first_row):
        for column_number, value in enumerate(record, 1):
            if not isinstance(value, (float, int, np.number)):
                continue
            cell = sheet[row_number - 1][column_number - 1]
            np.testing.assert_allclose(cell, value, rtol=1e-9, atol=1e-11,
                                       err_msg=f"Row {row_number}, column {column_number}")
            count += 1
    return count


def molecular_tables(output, sheets):
    regulatory = output / "regulatory"

    def load(name):
        return pd.read_csv(regulatory / (name + ".csv"))

    s4, s5 = sheets["S4 Regulatory predictions"], sheets["S5 Structural context"]
    gene = load("gene_tests")
    levels = gene["subset"].drop_duplicates().tolist()
    chromatin = load("chromatin_tests")
    origin = load("origin_regressions")
    sensitivity_origin = load("S4J_sensitivity")[origin.columns]
    isolated_counts = load("tissue_summary").pivot(index="modality", columns="dataset", values="n_alleles")
    isolated_counts = isolated_counts.assign(total=isolated_counts.sum(axis=1), distinct=isolated_counts.sum(axis=1)-2)
    blocks = [
        (s4, 6, isolated_counts.reset_index()),
        (s4, 20, origin),
        (s4, 34, load("junction_position")),
        (s4, 45, load("rna_tests")),
        (s4, 57, load("rna_backgrounds").sort_values(["subset", "scope", "background"])),
        (s4, 164, chromatin[chromatin.analysis_family.eq("PRIMARY")]),
        (s4, 217, chromatin[chromatin.analysis_family.eq("EFFECTIVE_MIN3")]),
        (s4, 270, gene[gene["subset"].eq(levels[2])].drop(columns="subset")),
        (s4, 301, gene[gene["subset"].eq(levels[0])]),
        (s4, 317, gene[gene["subset"].eq(levels[1])]),
        (s5, 6, load("structural_summary")),
        (s5, 22, load("structural_pairs")),
    ]
    for start, spec in [(338, "exclude_shared_linear"), (347, "exclude_shared_quadratic"),
                        (356, "all_quadratic")]:
        blocks.append((s4, start, sensitivity_origin[sensitivity_origin.specification.eq(spec)]))
    return sum(compare_numeric_block(*block) for block in blocks)


def arg589ter(output, s4):
    events = pd.read_csv(output / "regulatory/events.csv")
    hits = events[events.somatic_cdna.eq("c.1765C>T")].set_index("event")
    summaries = pd.read_csv(ROOT / "archive/analysis/source_data/Recovered_gene_centred_pair_summaries.csv")
    endpoints = [("native_all", "RNA expression"), ("native_hematopoietic", "RNA expression"),
                 ("native_all", "Splice-site usage")]
    profiles = []
    seen = set()
    for row in s4[290:295]:
        event = row[0].split("/")[-1].strip()
        seen.add(event)
        phase = hits.loc[event, "phase"]
        assert phase == row[1], event
        for text, (scope, modality) in zip(row[3:6], endpoints):
            matched = summaries[(summaries.pair_id == event) & (summaries.scope == scope) &
                                (summaries.modality == modality)]
            assert len(matched) <= 1
            if matched.empty:
                assert text == "ND (0 eligible tracks)", (event, text)
                value, tracks = None, 0
            else:
                value = float(matched.iloc[0].median_rescue_to_WT)
                tracks = int(matched.iloc[0].n_strong_rows)
                number = re.fullmatch(r"([-+]?\d+\.\d+) \(n=(\d+)\)", text)
                assert number and int(number[2]) == tracks, (event, text)
                np.testing.assert_allclose(float(number[1]), value, atol=.00005000001, rtol=0)
            profiles.append([event, hits.loc[event, "patient"], phase, scope, modality, value, tracks])
    assert seen == set(hits.index)
    labels = ["event", "patient", "phase", "scope", "modality", "proximity", "strong_tracks"]
    pd.DataFrame(profiles, columns=labels).to_csv(output / "regulatory/Arg589Ter_profiles.csv", index=False)
    return len(seen)


def verify_acquired_inventory(output, sheet):
    event, = [row for row in sheet if row[:1] == ("SGR057",)]
    allele, = [row for row in sheet if row[:2] == ("Somatic", "SOMVAR050")]
    cdna = "c.2053_2055del"
    protein = "p.Glu685del"
    genomic = "chr7:93133916:GTTC>G"
    assert event[3:6] == (cdna, protein, protein)
    assert allele[3:8] == (cdna, protein, 685, "In-frame deletion", genomic)
    missense, = [row for row in sheet if row[:2] == ("Somatic", "SOMVAR049")]
    assert missense[3:8] == ("c.3795A>T", "p.Lys1265Asn", 1265, "Missense", "chr7:93132177:T>A")

    events = pd.read_csv(output / "regulatory/events.csv").set_index("event")
    row = events.loc["SGR057"]
    assert row.somatic_cdna == cdna and row.somatic_label == protein
    assert row.structural and row.somatic_pos == 685 and row.sequence
    assert pd.isna(row.reason_without_sequence)
    assert events.loc[["SGR056", "SGR057"], "sequence"].all()

    summary = json.loads((output / "clinical/clinical_summary.json").read_text())["acquired"]
    counts = {key: summary[key] for key in ("events", "individuals", "exact_events", "exact_alleles")}
    assert tuple(counts.values()) == (57, 38, 57, 50)
    assert summary["observations_per_exact_allele"] == {"1": 46, "2": 3, "5": 1}
    assert summary["exact_types"] == {"Missense": 24, "Truncating / stop-gain": 12,
                                      "Frameshift": 11, "In-frame deletion": 3}
    return {**counts, "exact_singletons": 46, "isolated_acquired_predictions": 50,
            "combined_sequence_predictions": int(events.sequence.sum())}


def verify_revised_clinical(output, sheets):
    people = pd.DataFrame(sheets["S1 Cohort"][5:183], columns=sheets["S1 Cohort"][4])
    assert len(people) == 178 and people["Patient ID"].is_unique
    affected = people[people["Affected"].eq("Yes")]
    categories = {
        "Hematologic predominant": 80,
        "Mixed: Hematologic + Neurologic": 26,
        "Neurologic predominant": 21,
        "Mixed: Hematologic + Immune/inflammatory": 21,
        "Mixed: Hematologic + Neurologic + Immune/inflammatory": 19,
        "Immune/inflammatory predominant": 1,
        "Mixed: Neurologic + Immune/inflammatory": 1,
    }
    assert len(affected) == 169 and affected["Phenotype group"].value_counts().to_dict() == categories
    summary = json.loads((output / "clinical/clinical_summary.json").read_text())
    for key, expected in {"inflammation_individuals": (12, 14, 0, 112),
                          "inflammation_broad": (12, 14, 22, 112)}.items():
        counts = summary[key]
        fields = ("positive_trunc", "all_trunc", "positive_missense", "all_missense")
        assert tuple(counts[field] for field in fields) == expected
    roc_dir = output / "clinical/roc"
    roc = json.loads((roc_dir / "roc_summary.json").read_text())
    assert (roc["n"], roc["loss"], roc["comparator"], roc["n_adjusted"]) == (105, 57, 48, 105)
    np.testing.assert_allclose(roc["apparent_position_auc"], 0.7150950292397661, rtol=0, atol=1e-12)
    roster = pd.read_csv(roc_dir / "roc_input.csv").set_index("patient_id")
    assert roster.index.is_unique and roster.loc["STJ_P064", "chromosome_7_loss"] == 0
    assert roster.loc[["STJ_P093", "STJ_P094"], "chromosome_7_loss"].eq(1).all()
    s3 = sheets["S3 Clinical associations"]
    header, = [i for i, row in enumerate(s3) if row[:2] == ("Patient ID (S1)", "Primary protein label")]
    records = pd.DataFrame(s3[header + 1:header + 1 + len(roster)], columns=s3[header]).set_index("Patient ID (S1)")
    assert records.index.is_unique and set(records.index) == set(roster.index)
    reordered = roster.loc[records.index]
    np.testing.assert_array_equal(records["ROC outcome (loss = 1)"], reordered.chromosome_7_loss)
    np.testing.assert_array_equal(records["Position (aa)"], reordered.position_aa)
    np.testing.assert_array_equal(records["Age group"].str.lower(), reordered.age_group)
    np.testing.assert_allclose(records["Reported age for adjustment (years)"].astype(float), reordered.reported_age_years,
                               rtol=0, atol=1e-12, equal_nan=True)
    sensitivity = json.loads((roc_dir / "roc_sensitivities.json").read_text())
    checks = {
        "C30": summary["inflammation_broad"]["n"],
        "D30": summary["inflammation_broad"]["odds_ratio"],
        "E30": summary["inflammation_broad"]["p"],
        "C118": roc["loss"], "E118": roc["loss_median_aa"],
        "C119": roc["comparator"], "E119": roc["comparator_median_aa"],
        "C120": roc["n"], "E120": roc["exact_tied_rank_sum_two_sided_p"],
        "E121": roc["apparent_position_auc"],
        "F121": roc["apparent_position_auc_95_percentile_ci"][0],
        "G121": roc["apparent_position_auc_95_percentile_ci"][1],
        "E122": roc["cutoff"]["position_aa"], "C129": roc["n_adjusted"],
        "E129": roc["firth_adjusted_or_per_100_aa"],
        "F129": roc["firth_adjusted_or_95_percentile_ci"][0],
        "G129": roc["firth_adjusted_or_95_percentile_ci"][1],
        "C245": sensitivity["strict_normal"]["n"],
        "E245": sensitivity["strict_normal"]["p"],
        "E248": roc["held_out_auc_age_source"],
        "E249": roc["held_out_auc_position_age_source"],
    }
    all_carriers = output / "clinical/letter_roc_105"
    cv = pd.read_csv(all_carriers / "auc_105.csv", index_col=0)["AUC"]
    held = pd.read_csv(all_carriers / "held_out_105.csv")
    assert held.patient_id.is_unique and set(held.patient_id) == set(roster.index)
    for row, model in enumerate(["position", "position_age_source", "age_source"], 273):
        checks[f"E{row}"] = cv[model]
        checks[f"B{row}"] = len(held)
    omit = pd.read_csv(all_carriers / "omit_986.csv").iloc[0]
    checks.update({"B276": omit.N, "E276": omit.AUC, "F276": omit.optimal_residue})
    for row, key in [(246, "exclude_986"), (247, "exclude_partial")]:
        result = sensitivity[key]
        for column, value in zip("CEFG", (result["n"], result["odds_ratio"], *result["interval"])):
            checks[f"{column}{row}"] = value
    for row, key in enumerate(["affected_only", "without_P064", "marrow_documented", "single_diagnostic_sahoo"], 256):
        result = sensitivity[key]
        fields = ["n", "loss", "comparison", "auc", "p", "median_loss", "median_comparison"]
        checks.update({f"{column}{row}": result[field] for column, field in zip("BCDEFGH", fields)})
    for address, value in checks.items():
        r, c = coordinate_to_tuple(address)
        np.testing.assert_allclose(s3[r - 1][c - 1], value, rtol=1e-9, atol=1e-11,
                                   err_msg=f"S3!{address}")
    return {"affected_individuals": len(affected), "age_clinical_context": summary["age"]["clinical"]["n"],
            "pooled_roc_individuals": len(roster), "firth_adjusted_individuals": roc["n_adjusted"], "held_out_primary": len(held),
            "revised_result_cells": len(checks)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tables", type=Path, required=True)
    args = parser.parse_args()
    workbook = load_workbook(args.tables, read_only=True, data_only=True)
    try:
        sheets = {sheet.title: list(sheet.values) for sheet in workbook}
    finally:
        workbook.close()
    measured = json.loads((args.output / "clinical/table_statistics.json").read_text())
    clinical_cells = 0
    for sheet, values in measured.items():
        for address, value in values.items():
            r, c = coordinate_to_tuple(address)
            np.testing.assert_allclose(sheets[sheet][r-1][c-1], value, rtol=1e-9, atol=1e-10,
                                       err_msg=f"{sheet}!{address}")
            clinical_cells += 1
    stats = {"clinical_numeric_cells": clinical_cells,
             "sahoo_immune_positive_individuals": verify_sahoo_immune_phenotypes(sheets["S1 Cohort"]),
             "revised_clinical": verify_revised_clinical(args.output, sheets),
             "acquired_inventory": verify_acquired_inventory(args.output, sheets["S2 Variants"]),
             "molecular_numeric_cells": molecular_tables(args.output, sheets),
             "age_endpoint_observations": verify_membership(args.output, ROOT / "data/endpoint_membership.csv"),
             "Arg589Ter_profiles": arg589ter(args.output, sheets["S4 Regulatory predictions"])}
    (args.output / "table_verification.json").write_text(json.dumps(stats, indent=2) + "\n")
    print(json.dumps(stats))


if __name__ == "__main__":
    main()
