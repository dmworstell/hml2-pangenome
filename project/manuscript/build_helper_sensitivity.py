"""Reproduce the two retained helper-dependence sensitivity tables."""
from pathlib import Path
import csv
import math

OUT = Path(__file__).resolve().parent / 'delta292_helper_abc_v1'


def helper_rows():
    rows = []
    for focal_fraction in (10 / 16, 17 / 23):
        for moi in (.5, 1, 2, 3, 5, 10):
            helper_fraction = 1 - focal_fraction
            complementation = 1 - math.exp(-moi * helper_fraction)
            rows.append({'target_focal_fraction': focal_fraction,
                'helper_fraction': helper_fraction, 'moi': moi,
                'complementation_probability': complementation,
                'minimum_replication_gain_for_equilibrium': 1 / complementation})
    return rows


def neutral_rows():
    rows = []
    for population in (10, 50, 100, 1000, 10000, 100000):
        for copies in (1, 5):
            for target in (10 / 16, 17 / 23):
                rows.append({'effective_population_size': population,
                    'starting_copies': copies, 'starting_fraction': copies / population,
                    'target_fraction': target,
                    'martingale_upper_bound_probability_ever_reach_target':
                    min(1, copies / population / target),
                    'note': 'Upper bound for an unbiased neutral-frequency martingale.'})
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for name, rows in [('helper_equilibrium_requirements.tsv', helper_rows()),
                       ('neutral_founder_hitting_bounds.tsv', neutral_rows())]:
        with (OUT / name).open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)


if __name__ == '__main__':
    main()
