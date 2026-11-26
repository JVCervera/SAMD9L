"""Regenerate the publication artwork for Figure 1.

The template holds the fixed variant landscape and axis geometry. Analysis
outputs populate the clinical-context scatter and ROC curve.
"""

import argparse
from collections import Counter
from copy import deepcopy
import json
from pathlib import Path
import subprocess

from lxml import etree
import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import spearmanr


NS = {"svg": "http://www.w3.org/2000/svg"}
SUPERSCRIPT = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")

def node(root, identifier):
    return root.xpath(f'.//*[@id="{identifier}"]')[0]


def label(root, identifier):
    return node(root, identifier).find("svg:text", NS)


def write_path(path, coordinates):
    path.set("d", " ".join(f"{cmd} {x:.6f} {y:.6f}" for cmd, x, y in coordinates))


def scale(axis, direction):
    """Map data to in-panel SVG coordinates from the fixed artwork's ticks."""
    anchors = []
    for tick in axis.xpath(f'.//*[contains(@id,"{direction}tick_")]'):
        text = tick.find(".//svg:text", NS)
        marker = tick.find(".//svg:use", NS)
        if text is not None and marker is not None:
            anchors.append((float(text.text), float(marker.get(direction))))
    return np.polyfit(*np.array(anchors).T, deg=1)


def formatted_p(value):
    mantissa, power = f"{value:.2e}".split("e")
    return f"{mantissa} × 10{str(int(power)).translate(SUPERSCRIPT)}"


def age_scatter(svg, rows):
    ages = rows.loc[rows.endpoint.eq("clinical")].copy()
    assert ages.patient.is_unique and ages.age.notna().all() and ages.position.notna().all()
    axis = node(svg, "B_axes_1")
    x_map, y_map = (scale(axis, dimension) for dimension in "xy")
    for i, group in enumerate(("Pediatric", "Adult"), start=1):
        markers = node(svg, f"B_PathCollection_{i}").find("svg:g", NS)
        points = ages.loc[ages.group.eq(group)]
        template = deepcopy(markers[0])
        for child in list(markers):
            markers.remove(child)
        for point in points.itertuples():
            marker = deepcopy(template)
            marker.set("x", f"{np.polyval(x_map, point.position):.6f}")
            marker.set("y", f"{np.polyval(y_map, point.age):.6f}")
            markers.append(marker)
        label(svg, f"B_text_{12 + i}").text = f"{group} (n = {len(points)})"
        assert len(markers) == len(points)
    result = spearmanr(ages.position, ages.age)
    rho = f"{result.statistic:.3f}".replace("-", "−")
    label(svg, "B_text_184").text = f"Spearman ρ = {rho},  P = {formatted_p(result.pvalue)}"
    return ages.rename(columns={
        "patient": "Patient_ID", "position": "Primary_germline_residue",
        "age": "Clinical_context_age_years", "group": "Stratum",
    })[["Patient_ID", "Primary_germline_residue", "Clinical_context_age_years", "Stratum"]]


def roc_panel(svg, roc, roc_summary):
    axis = node(svg, "C_axes_2")
    x_map, y_map = (scale(axis, dimension) for dimension in "xy")
    path = node(svg, "C_line2d_34").find("svg:path", NS)
    x = np.polyval(x_map, roc.false_positive_rate)
    y = np.polyval(y_map, roc.true_positive_rate)
    write_path(path, [("M" if i == 0 else "L", xp, yp) for i, (xp, yp) in enumerate(zip(x, y))])
    lo, hi = roc_summary["apparent_position_auc_95_percentile_ci"]
    auc = roc_summary["apparent_position_auc"]
    assert np.isclose(np.trapezoid(roc.true_positive_rate, roc.false_positive_rate), auc)
    label(svg, "C_text_187").text = f"Apparent AUC {auc:.3f} (95% CI {lo:.3f}–{hi:.3f})"


def write_tree(svg, target):
    etree.ElementTree(svg).write(target, encoding="UTF-8", xml_declaration=True)


def export(output, with_tiff=False):
    pdf = output.with_suffix(".pdf")
    subprocess.run(["inkscape", str(output), f"--export-filename={pdf}"], check=True)
    if with_tiff:
        png = output.with_suffix(".working.png")
        subprocess.run(["inkscape", str(output), "--export-dpi=1000",
                        "--export-background=#ffffff", "--export-background-opacity=1",
                        f"--export-filename={png}"], check=True)
        with Image.open(png) as render:
            rgb = render.convert("RGB")
        rgb.save(output.with_suffix(".tif"), compression="tiff_lzw", dpi=(1000, 1000))
        png.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", type=Path,
                        default=Path(__file__).with_name("Figure_1_layout.svg"))
    parser.add_argument("--age-rows", required=True, type=Path)
    parser.add_argument("--roc-curve", required=True, type=Path)
    parser.add_argument("--roc-summary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--svg-only", action="store_true",
                        help="Write vector artwork and point coordinates without Inkscape conversion.")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    svg = etree.parse(args.template).getroot()
    points = age_scatter(svg, pd.read_csv(args.age_rows, float_precision="round_trip"))
    assert Counter(points.Stratum) == Counter({"Pediatric": 105, "Adult": 33})
    points.to_csv(args.output / "Figure_1B_points.csv", index=False)
    roc_panel(svg, pd.read_csv(args.roc_curve), json.loads(args.roc_summary.read_text()))

    main_svg = args.output / "Figure_1.svg"
    write_tree(svg, main_svg)
    panel = deepcopy(svg[1])
    panel.attrib.pop("x")
    panel.attrib.pop("y")
    for attr in ("width", "height"):
        panel.attrib[attr] += "pt"
    panel_svg = args.output / "Figure_1B_Clinical_context_age_scatter.svg"
    write_tree(panel, panel_svg)
    if not args.svg_only:
        export(main_svg, with_tiff=True)
        export(panel_svg)


if __name__ == "__main__":
    main()
