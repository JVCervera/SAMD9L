"""Reproduce the SAMD9L Letter from its supplementary tables and model exports."""
import argparse
import csv
from collections import Counter
import shutil
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def compare_json(new, expected):
    import numpy as np
    if isinstance(expected, dict):
        assert new.keys() == expected.keys(), new.keys() ^ expected.keys()
        for key, value in expected.items():
            compare_json(new[key], value)
    elif isinstance(expected, list):
        assert len(new) == len(expected)
        for a, b in zip(new, expected):
            compare_json(a, b)
    elif isinstance(expected, (int, float)):
        np.testing.assert_allclose(new, expected, rtol=1e-9, atol=1e-11)
    else:
        assert new == expected, (new, expected)


def compare_reference(output):
    import pandas as pd
    reference = ROOT / "reference"
    for directory in ("clinical", "regulatory"):
        for path in sorted((reference / directory).rglob("*")):
            new = output / path.relative_to(reference)
            if path.suffix == ".csv":
                pd.testing.assert_frame_equal(pd.read_csv(new), pd.read_csv(path), check_dtype=False,
                                              check_exact=False, rtol=1e-9, atol=1e-11)
            elif path.suffix == ".json":
                compare_json(json.loads(new.read_text()), json.loads(path.read_text()))


def copy_submission_figures(output, render_svg=False):
    reference = ROOT / "reference/figures"
    filenames = {f"{name}.{extension}" for name in
                 ("Figure_1", "Figure_2", "Figure_1B_Clinical_context_age_scatter", *[f"Figure_S{k}" for k in range(1, 9)])
                 for extension in ("pdf", "svg")}
    filenames.update({"Figure_1.tif", "Figure_2.tif"})
    if {path.name for path in reference.iterdir()} != filenames:
        raise ValueError("Unexpected inventory of submission figures")
    output.mkdir(parents=True, exist_ok=True)
    for filename in sorted(filenames):
        shutil.copyfile(reference / filename, output / filename)
        if sha256(output / filename) != sha256(reference / filename):
            raise ValueError(f"Figure checksum mismatch: {filename}")
    if render_svg:
        rendered = output.parent / "figures_rendered"
        rendered.mkdir()
        for source in sorted(output.glob("*.svg")):
            subprocess.run(["inkscape", str(source),
                            f"--export-filename={rendered / (source.stem + '.pdf')}"], check=True)
    return sorted(filenames)


def verify_figure_1b(analysis_dir):
    from scipy.stats import spearmanr

    with (analysis_dir / "age_rows.csv").open(newline="") as stream:
        clinical = [row for row in csv.DictReader(stream) if row["endpoint"] == "clinical"]
    expected = Counter((row["patient"], int(row["position"]), round(float(row["age"]), 12), row["group"])
                       for row in clinical)
    with (ROOT / "data/Figure_1B_points.csv").open(newline="") as stream:
        supplied = Counter((row["Patient_ID"], int(row["Primary_germline_residue"]),
                            round(float(row["Clinical_context_age_years"]), 12), row["Stratum"])
                           for row in csv.DictReader(stream))
    if supplied != expected:
        raise ValueError("Figure 1B plotted points disagree with the clinical-context records")
    result = spearmanr([float(row["age"]) for row in clinical],
                       [int(row["position"]) for row in clinical])
    return {"n": len(clinical), "rho": float(result.statistic), "p": float(result.pvalue)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "reproduced",
                        help="New or empty directory for this run.")
    parser.add_argument("--workbook", type=Path, default=ROOT / "data/SAMD9L_Supplementary_Tables.xlsx")
    parser.add_argument("--skip-track-checks", action="store_true",
                        help="Skip track-level verification and use the supplied per-allele and paired summaries.")
    parser.add_argument("--render-svg", action="store_true",
                        help="Render SVGs to a separate directory with local Inkscape.")
    parser.add_argument("--copy-figures-only", action="store_true",
                        help="Copy the exact submission graphics and skip all analyses.")
    parser.add_argument("--redraw-figures", action="store_true",
                        help="Redraw figures for data inspection with Matplotlib in a separate directory.")
    args = parser.parse_args()
    output, workbook = args.output.resolve(), args.workbook.resolve()
    if output.exists() and any(output.iterdir()):
        parser.error(f"Output directory contains files: {output}. Choose a new --output directory.")
    output.mkdir(parents=True, exist_ok=True)
    if args.copy_figures_only:
        copy_submission_figures(output / "figures", args.render_svg)
        return
    if not workbook.is_file():
        parser.error(f"Missing workbook: {workbook}")
    if not args.skip_track_checks and not (ROOT / "restore_tracks.py").is_file():
        parser.error("The compact archive uses frozen model summaries. Run with --skip-track-checks.")
    logs = output / "logs"
    logs.mkdir()
    inputs = ROOT / "archive"
    raw = inputs / "analysis"
    scripts = ROOT / "scripts"
    env = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
           "OMP_NUM_THREADS": "1", "MPLBACKEND": "Agg", "PYTHONOPTIMIZE": ""}
    steps = []

    def execute(label, script, *arguments):
        print(label, flush=True)
        command = [sys.executable, str(script), *map(str, arguments)]
        logfile = logs / (label + ".log")
        with logfile.open("w") as stream:
            status = subprocess.run(command, env=env, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
        if status.returncode:
            print(logfile.read_text(errors="replace")[-3000:], file=sys.stderr)
            raise SystemExit(f"Failed {label}. See {logfile}")
        steps.append(label)

    if not args.skip_track_checks:
        from restore_tracks import restore_tracks
        restore_tracks()
        for kind in ["isolated", "gene", "paired"]:
            execute("track_" + kind, raw / ("verify_" + kind + "_tracks.py"), "--output", output / "track_validation")
        execute("splice_junction_tracks", raw / "junctions/verify.py")
    execute("clinical", scripts / "clinical/clinical.py", "--workbook", workbook, "--output", output / "clinical")
    roc_output = output / "clinical/roc"
    roc_output.mkdir()
    roc_input = roc_output / "roc_input.csv"
    roc_input.write_bytes((ROOT / "data/pooled_roc_input.csv").read_bytes())
    execute("pooled_roc", scripts / "clinical/pooled_roc.py", "--input", roc_input,
            "--output-dir", roc_output)
    execute("roc_sensitivities", scripts / "clinical/roc_sensitivities.py",
            "--input", ROOT / "data/roc_audited_records.csv",
            "--baseline", output / "clinical/sahoo_baseline.csv",
            "--output", roc_output / "roc_sensitivities.json")
    execute("roc_posthoc", scripts / "clinical/roc_posthoc.py", "--input", roc_input,
            "--output", roc_output / "roc_posthoc.json")
    execute("letter_roc_105", scripts / "clinical/letter_cutoffs.py", "--output", output / "clinical/letter_roc_105")
    execute("molecular", scripts / "regulatory/recalculate.py", "--archive", inputs,
            "--tables", workbook, "--output", output / "regulatory")
    assets = copy_submission_figures(output / "figures", args.render_svg)
    clinical_context = verify_figure_1b(output / "clinical")
    steps.extend(["submission_figures", "figure_1b_points"])
    if args.redraw_figures:
        execute("figures_redrawn", scripts / "figures/build_letter_figures.py", "--workbook", workbook,
                "--source", raw / "source_data", "--clinical", output / "clinical",
                "--regulatory", output / "regulatory", "--events", output / "regulatory/events.csv",
                "--roc-summary", roc_output / "roc_summary.json", "--output", output / "figures_redrawn")
    execute("table_values", scripts / "verify_letter_results.py", "--output", output, "--tables", workbook)
    execute("figure_values", ROOT / "validation/verify_figure_data.py", output)
    compare_reference(output)
    result = {"status": "PASS", "steps": steps + ["reference_comparison"],
              "tables_sha256": sha256(workbook), "track_checks": not args.skip_track_checks,
              "submission_figure_assets": assets, "figure_1b_clinical_context": clinical_context,
              "table_checks": json.loads((output / "table_verification.json").read_text())}
    (output / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
