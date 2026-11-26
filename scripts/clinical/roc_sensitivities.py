"""Recompute the roster, specimen and adjusted sensitivities in Supplementary Table S3."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from firth_validation import contrast, fit_model, profile
from pooled_roc import auc, design_matrix, draw_indices, exact_rank_sum, firth, groups


PARTIAL_LOSS = {
    "STJ_P018", "STJ_P022", "STJ_P038", "STJ_P040",
    "STJ_P086", "STJ_P097", "STJ_P105", "VUELTA_N697Y",
}


def position_stats(data, outcome="chromosome_7_loss"):
    position = data.position_aa.to_numpy()
    positive = data[outcome].eq(1).to_numpy()
    return {
        "n": len(data),
        "loss": int(positive.sum()),
        "comparison": int((~positive).sum()),
        "auc": auc(positive, position),
        "p": exact_rank_sum(positive, position),
        "median_loss": float(np.median(position[positive])),
        "median_comparison": float(np.median(position[~positive])),
    }


def exclude_986(records):
    sample = records[records.position_aa.ne(986)]
    status = sample.chromosome_7_loss.to_numpy()
    design = design_matrix(
        sample.age_group.eq("adult").to_numpy(),
        sample.source_group.to_numpy(),
        sample.position_aa.to_numpy(),
    )
    blocks = groups(sample.get("resampling_group", sample.protein_variant).to_numpy())
    generator = np.random.default_rng(20261123)
    draws = []
    for _ in range(3000):
        drawn = draw_indices(generator, blocks)
        try:
            draws.append(np.exp(firth(status[drawn], design[drawn])[1]))
        except (RuntimeError, ValueError, np.linalg.LinAlgError):
            continue
    return {
        "n": len(sample),
        "loss": int(status.sum()),
        "odds_ratio": float(np.exp(firth(status, design)[1])),
        "interval": np.quantile(draws, [.025, .975]).tolist(),
        "valid_draws": len(draws),
    }


def exclude_partial(records):
    sample = records.loc[~records.index.isin(PARTIAL_LOSS)]
    rows = [
        dict(
            patient_id=patient,
            position=row.position_aa,
            age_reported_mixed=row.reported_age_years,
            age_group=row.age_group.capitalize(),
            source_group=row.source_group,
            cytogenetic_class="loss7" if row.chromosome_7_loss else "other",
        )
        for patient, row in sample.iterrows()
    ]
    fit = fit_model(rows, ["position", "age_group", "source"])
    limits = profile(fit, contrast(fit, "position_per100"))
    return {
        "n": len(sample),
        "loss": int(sample.chromosome_7_loss.sum()),
        "odds_ratio": limits["OR"],
        "interval": limits["profile95_OR"],
    }


def run(path, baseline):
    records = pd.read_csv(path).set_index("patient_id")
    assert len(records) == 105 and records.index.is_unique
    baseline = pd.read_csv(baseline).rename(
        columns={"position": "position_aa", "loss": "chromosome_7_loss"},
    )
    return {
        "strict_normal": position_stats(records[records.cytogenetic_group.ne("Other abnormality")]),
        "affected_only": position_stats(records[records.clinical_manifestations.eq("Yes")]),
        "without_P064": position_stats(records.drop(index="STJ_P064")),
        "marrow_documented": position_stats(records[records.marrow_documented], "marrow_loss"),
        "single_diagnostic_sahoo": position_stats(baseline),
        "exclude_986": exclude_986(records),
        "exclude_partial": exclude_partial(records),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summary = run(args.input, args.baseline)
    args.output.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
