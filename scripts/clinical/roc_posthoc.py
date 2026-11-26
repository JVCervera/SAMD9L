import argparse
from collections import Counter
import csv
import json
from pathlib import Path

import numpy as np

from pooled_roc import auc, draw_indices, exact_rank_sum, groups, youden


def describe_cutoff(outcome, position, cutoff):
    positive = position >= cutoff
    tp = int((positive & outcome).sum())
    fp = int((positive & ~outcome).sum())
    tn = int((~positive & ~outcome).sum())
    fn = int((~positive & outcome).sum())
    sensitivity = tp / (tp + fn)
    specificity = tn / (tn + fp)
    return {
        'cutoff': cutoff, 'tp': tp, 'fp': fp, 'tn': tn, 'fn': fn,
        'sensitivity': sensitivity, 'specificity': specificity,
        'youden': sensitivity + specificity - 1,
    }


def run(source):
    with source.open(newline='', encoding='utf-8') as handle:
        records = list(csv.DictReader(handle))
    assert len({record['patient_id'] for record in records}) == len(records)
    outcome = np.array([record['chromosome_7_loss'] == '1' for record in records])
    position = np.array([int(record['position_aa']) for record in records])
    variant = np.array([record['protein_variant'] for record in records])
    omit = 'R1281K'
    keep = variant != omit
    sub_y, sub_x = outcome[keep], position[keep]
    block_labels = np.array([row.get("resampling_group") or row["protein_variant"] for row in records])
    bootstrap_groups = groups(block_labels[keep])
    rng = np.random.default_rng(np.random.SeedSequence(20261123).spawn(2)[0])
    boot = []
    for _ in range(3000):
        indices = draw_indices(rng, bootstrap_groups)
        boot.append(auc(sub_y[indices], sub_x[indices]))
    scan = []
    for label, count in Counter(variant).items():
        retained = variant != label
        scan.append({'variant': label, 'removed': count,
                     'auc': float(auc(outcome[retained], position[retained]))})
    scan.sort(key=lambda entry: (entry['auc'], entry['variant']))
    return {
        'main': {'n': len(records), 'auc': float(auc(outcome, position))},
        'omit_R1281K': {
            'n': int(keep.sum()), 'loss': int(sub_y.sum()), 'comparator': int((~sub_y).sum()),
            'omitted_loss': int(outcome[~keep].sum()), 'omitted_comparator': int((~outcome[~keep]).sum()),
            'auc': float(auc(sub_y, sub_x)), 'ci': np.nanquantile(boot, [.025, .975]).tolist(),
            'p': exact_rank_sum(sub_y.astype(int), sub_x),
            'median_loss': float(np.median(sub_x[sub_y])),
            'median_comparator': float(np.median(sub_x[~sub_y])),
            'optimal_cutoff': youden(sub_y.astype(int), sub_x)['position_aa'],
        },
        'cutoffs': [describe_cutoff(outcome, position, t) for t in (986, 876)],
        'influence_scan': scan,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Post hoc influence of protein variants and operating points.')
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(args.input), indent=2) + '\n', encoding='utf-8')
