#!/usr/bin/env python3
"""ROC, exact tied rank-sum, adjusted Firth OR, and held-out predictions."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize
from scipy.special import expit
from scipy.stats import rankdata


def auc(y, score):
    positive = y.astype(bool)
    n_loss = int(positive.sum())
    n_comparator = len(y) - n_loss
    if n_loss == 0 or n_comparator == 0:
        return np.nan
    return (rankdata(score)[positive].sum() - n_loss * (n_loss + 1) / 2) / (
        n_loss * n_comparator
    )


def youden(y, position):
    cutoffs = np.unique(position)
    positive = position[:, None] >= cutoffs
    upper = positive.sum(axis=0)
    valid = (upper >= 10) & (len(y) - upper >= 10)
    if not valid.any() or not y.any() or y.all():
        return None
    true_positive = y @ positive
    true_negative = (1 - y) @ (~positive)
    sensitivity = true_positive / y.sum()
    specificity = true_negative / (len(y) - y.sum())
    # Integer numerator preserves mathematically exact ties.
    index = true_positive * (len(y) - y.sum()) + true_negative * y.sum()
    best = np.argmax(np.where(valid, index, -1))
    return {
        "position_aa": int(cutoffs[best]),
        "youden_j": float(sensitivity[best] + specificity[best] - 1),
        "sensitivity": float(sensitivity[best]),
        "specificity": float(specificity[best]),
        "loss_at_or_above": int(true_positive[best]),
        "n_at_or_above": int(upper[best]),
        "loss_below": int(y.sum() - true_positive[best]),
        "n_below": int(len(y) - upper[best]),
    }


def exact_rank_sum(y, position):
    """Count every group assignment with the observed group sizes."""
    ranks = np.rint(2 * rankdata(position)).astype(int)
    n_group = min(int(y.sum()), len(y) - int(y.sum()))
    selected = y.astype(bool) if y.sum() == n_group else ~y.astype(bool)
    observed = ranks[selected].sum()
    null_center = n_group * (len(y) + 1)
    limit = n_group * ranks.max() + 1
    counts = np.zeros((n_group + 1, limit), dtype=object)
    counts[0, 0] = 1
    for rank in ranks:
        counts[1:, rank:] += counts[:-1, :-rank]
    distance = abs(observed - null_center)
    extreme = np.abs(np.arange(limit) - null_center) >= distance
    return float(sum(counts[n_group, extreme]) / sum(counts[n_group]))


def firth(y, design):
    """Fit Jeffreys-penalized logistic regression with its analytic score."""
    if np.unique(y).size < 2:
        raise ValueError("Firth fitting requires both endpoint classes.")
    coefficient_count = design.shape[1]
    if np.linalg.matrix_rank(design) != coefficient_count:
        # Bootstrap may omit a source entirely or leave age aliased with source.
        keep = []
        for i in range(design.shape[1]):
            candidate = keep + [i]
            if np.linalg.matrix_rank(design[:, candidate]) == len(candidate):
                keep = candidate
    else:
        keep = list(range(design.shape[1]))
    design = design[:, keep]

    def criterion(beta):
        predictor = design @ beta
        probability = expit(predictor)
        weight = probability * expit(-predictor)
        information = design.T @ (weight[:, None] * design)
        sign, log_det = np.linalg.slogdet(information)
        if sign <= 0:
            return np.inf, np.zeros_like(beta)
        solved = np.linalg.solve(information, design.T).T
        leverage = weight * (solved * design).sum(axis=1)
        log_likelihood = y @ predictor - np.logaddexp(0, predictor).sum()
        score = design.T @ (y - probability + leverage * (0.5 - probability))
        return -log_likelihood - 0.5 * log_det, -score

    result = minimize(
        criterion, np.zeros(len(keep)), jac=True, method="BFGS", options={"gtol": 1e-8}
    )
    if not result.success:
        gradient = np.max(np.abs(result.jac))
        if gradient > 1e-5:
            raise RuntimeError(f"Firth convergence failure, max absolute score {gradient:g}")
    full = np.zeros(coefficient_count)
    full[keep] = result.x
    return full


def design_matrix(adult, source, position=None):
    columns = [np.ones(len(adult))]
    # Include position before the nuisance variables if a draw induces aliasing.
    if position is not None:
        columns.append((position - 986) / 100)
    columns.append(adult.astype(float))
    for label in ("Sahoo A", "Sahoo B", "Nagata", "Other literature"):
        columns.append((source == label).astype(float))
    return np.column_stack(columns)


def groups(variant):
    return [np.flatnonzero(variant == name) for name in np.unique(variant)]


def draw_indices(rng, blocks):
    chosen = rng.integers(0, len(blocks), size=len(blocks))
    return np.concatenate([blocks[i] for i in chosen])


def held_out(y, design, blocks):
    predicted = np.empty(len(y))
    for tested in blocks:
        train = np.ones(len(y), dtype=bool)
        train[tested] = False
        beta = firth(y[train], design[train])
        train_prediction = design[train] @ beta
        prevalence = y[train].mean()
        intercept_adjustment = brentq(
            lambda shift: expit(train_prediction + shift).mean() - prevalence, -50, 50
        )
        predicted[tested] = expit(design[tested] @ beta + intercept_adjustment)
    return predicted


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path(__file__).with_name("roc_input.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--bootstrap", type=int, default=3000)
    parser.add_argument("--seed", type=int, default=20261123)
    args = parser.parse_args()
    with open(args.input, encoding="utf-8", newline="") as source:
        rows = list(csv.DictReader(source))
    y = np.array([int(row["chromosome_7_loss"]) for row in rows])
    position = np.array([int(row["position_aa"]) for row in rows])
    variant = np.array([row["protein_variant"] for row in rows])
    age_group = np.array([row["age_group"] for row in rows])
    if set(age_group) - {"adult", "pediatric"}:
        raise ValueError("A documented age group is required for adjusted models.")
    adult = age_group == "adult"
    source = np.array([row["source_group"] for row in rows])
    design = design_matrix(adult, source, position)
    reduced = design_matrix(adult, source)
    blocks = groups(np.array([row.get("resampling_group") or row["protein_variant"] for row in rows]))
    rng, adjusted_rng = [np.random.default_rng(seed) for seed in
                         np.random.SeedSequence(args.seed).spawn(2)]

    cutoff = youden(y, position)
    adjusted_or = np.exp(firth(y, design)[1])
    unadjusted_bootstrap = np.empty(args.bootstrap)
    adjusted_bootstrap = np.empty(args.bootstrap)
    cutoff_bootstrap = np.full(args.bootstrap, np.nan)
    for draw in range(args.bootstrap):
        idx = draw_indices(rng, blocks)
        unadjusted_bootstrap[draw] = auc(y[idx], position[idx])
        selected = youden(y[idx], position[idx])
        if selected is not None:
            cutoff_bootstrap[draw] = selected["position_aa"]
        idx = draw_indices(adjusted_rng, blocks)
        try:
            adjusted_bootstrap[draw] = np.exp(firth(y[idx], design[idx])[1])
        except (ValueError, RuntimeError, np.linalg.LinAlgError):
            adjusted_bootstrap[draw] = np.nan

    pred_with_position = held_out(y, design, blocks)
    pred_without_position = held_out(y, reduced, blocks)
    summary = {
        "n": len(y), "loss": int(y.sum()), "comparator": int(len(y) - y.sum()),
        "variant_groups": len(blocks),
        "loss_median_aa": float(np.median(position[y == 1])),
        "comparator_median_aa": float(np.median(position[y == 0])),
        "exact_tied_rank_sum_two_sided_p": exact_rank_sum(y, position),
        "apparent_position_auc": auc(y, position),
        "apparent_position_auc_95_percentile_ci": np.nanquantile(unadjusted_bootstrap, [.025, .975]).tolist(),
        "cutoff": cutoff,
        "bootstrap_selected_same_cutoff": int((cutoff_bootstrap == cutoff["position_aa"]).sum()),
        "bootstrap_valid_cutoff": int(np.isfinite(cutoff_bootstrap).sum()),
        "bootstrap_valid_apparent_auc": int(np.isfinite(unadjusted_bootstrap).sum()),
        "bootstrap_replicates_requested": args.bootstrap, "bootstrap_seed": args.seed,
        "n_adjusted": len(y), "loss_adjusted": int(y.sum()), "adjusted_variant_groups": len(blocks),
        "firth_adjusted_or_per_100_aa": float(adjusted_or),
        "firth_adjusted_or_95_percentile_ci": np.nanquantile(adjusted_bootstrap, [.025, .975]).tolist(),
        "bootstrap_valid_adjusted_or": int(np.isfinite(adjusted_bootstrap).sum()),
        "held_out_auc_position_age_source": auc(y, pred_with_position),
        "held_out_auc_age_source": auc(y, pred_without_position),
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with open(args.output_dir / "roc_summary.json", "w", encoding="utf-8") as output:
        json.dump(summary, output, indent=2)
        output.write("\n")
    with open(args.output_dir / "held_out_predictions.csv", "w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output)
        writer.writerow(["patient_id", "loss", "with_position", "without_position"])
        for row, observed, complete, reduced_score in zip(
            rows, y, pred_with_position, pred_without_position
        ):
            writer.writerow([row["patient_id"], observed, complete, reduced_score])
    with open(args.output_dir / "roc_curve.csv", "w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output)
        writer.writerow(["threshold_aa", "false_positive_rate", "true_positive_rate"])
        writer.writerow(["Inf", 0, 0])
        for threshold in np.unique(position)[::-1]:
            selected = position >= threshold
            writer.writerow([
                threshold, (selected & (y == 0)).sum() / (y == 0).sum(),
                (selected & (y == 1)).sum() / (y == 1).sum(),
            ])
    with open(args.output_dir / "bootstrap_draws.csv", "w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output)
        writer.writerow(["replicate", "apparent_auc", "youden_cutoff_aa", "adjusted_or_per_100_aa"])
        for index, draw in enumerate(
            zip(unadjusted_bootstrap, cutoff_bootstrap, adjusted_bootstrap), start=1
        ):
            writer.writerow([index, *draw])
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
