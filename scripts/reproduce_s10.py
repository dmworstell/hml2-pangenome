"""Replot final Figure S10 from retained, record-level figure witnesses."""
from collections import Counter
from pathlib import Path
import argparse
import csv
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'manuscript_figures/python/figures'))
import atypical_fusion_figure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'outputs/Figure_S10')
    args = parser.parse_args()
    data = ROOT / 'Supplementary_Data/Figure_S9'
    summary = json.loads((data / 'figure_counts_and_provenance.json').read_text())
    counts = []
    for name, field, expected in [
        ('highest_annotation_tiers.tsv', 'highest_tier', 'annotation_tiers'),
        ('gag_region_longest_frame_witnesses.tsv', 'regional_reach', 'gag_region_reach'),
    ]:
        with (data / name).open(newline='') as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        assert len(rows) == len({row['ID_Full'] for row in rows}), name
        observed = Counter(row[field] for row in rows)
        assert observed == summary[expected], name
        counts.append(observed)
    args.output.mkdir(parents=True, exist_ok=True)
    atypical_fusion_figure.render_figure(*counts, args.output, 'Figure_S10')
    print('PASS: Figure S10 counts match both retained witness tables and provenance.')


if __name__ == '__main__':
    main()
