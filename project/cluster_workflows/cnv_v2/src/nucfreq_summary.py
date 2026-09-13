"""The site-counting and depth-reliability rules behind the July 2026 CNV table.

The frozen table used BED start coordinates directly against inclusive numeric
body_start/body_end bounds taken from the depth track. This module deliberately
reproduces that historical interval. It does not silently shift coordinates.
"""
from collections import defaultdict
import csv
import math
from pathlib import Path
from statistics import median


def count_region_sites(lines, body_start, body_end):
    """Count covered/mixed-base sites exactly as the original summary reducer.

    Standard input is rustybam '#chr start end A C G T' output. The original
    four-column chrom/position/top-count/second-count variant is also retained.
    Zero-count positions do not enter the denominator. There is no independent
    minimum total-depth cutoff. Mixed-base sites require d2>=3 and
    d2/(d1+d2)>=0.15, where d1/d2 are the two largest base counts.
    """
    total = 0
    mixed = 0
    second_sum = 0
    for line in lines:
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 4:
            continue
        try:
            position = int(fields[1])
        except ValueError:
            continue
        if not body_start <= position <= body_end:
            continue
        counts = []
        for field in fields[2:]:
            try:
                counts.append(int(field))
            except ValueError:
                pass
        if len(counts) >= 4:
            first, second = sorted(counts[-4:], reverse=True)[:2]
        elif len(counts) >= 2:
            first, second = counts[:2]
        else:
            continue
        if first + second <= 0:
            continue
        total += 1
        if second >= 3 and second / (first + second) >= 0.15:
            mixed += 1
            second_sum += second
    return {
        'nucfreq_total_sites': total,
        'nucfreq_het_sites': mixed,
        'nucfreq_het_frac': mixed / total if total else None,
        'nucfreq_mean_alt_depth': second_sum / mixed if mixed else None,
    }


def check_retained_summary(path):
    """Recalculate the retained table's site fractions and depth-reliability flags.

    This checks the archived summary, not a recount from raw BAMs or BED files.
    Low-baseline regions lose depth ratios only, not their NucFreq counts.
    """
    with Path(path).open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    flanks = defaultdict(list)
    for row in rows:
        flank = _number(row['flank_median_all'])
        if flank is not None and flank > 0:
            flanks[row['sample']].append(flank)
    numeric = []
    low_count = 0
    low_numeric = 0
    for row in rows:
        identity = (row['locus'], row['sample'], row['region'])
        reference = median(flanks[row['sample']]) if flanks[row['sample']] else None
        if not _same_four_decimals(reference, row['sample_ref_cov']):
            raise ValueError(f'Sample reference mismatch: {identity}')
        flank = _number(row['flank_median_all']) or 0
        low = flank <= 0 or (reference is not None and reference > 0 and flank < 0.4 * reference)
        if str(low) != row['low_baseline']:
            raise ValueError(f'Low-baseline flag mismatch: {identity}')
        low_count += int(low)
        body = _number(row['body_median_all'])
        expected_depth = body / (flank / 2) if not low and flank > 0 and body is not None else None
        if not _same_four_decimals(expected_depth, row['hapcopies_all']):
            raise ValueError(f'Haploid-copy depth mismatch: {identity}')
        total = int(row['nucfreq_total_sites'])
        mixed = int(row['nucfreq_het_sites'])
        if not 0 <= mixed <= total:
            raise ValueError(f'Invalid site counts: {identity}')
        fraction = mixed / total if total else None
        if not _same_four_decimals(fraction, row['nucfreq_het_frac']):
            raise ValueError(f'Regional site fraction mismatch: {identity}')
        if fraction is not None:
            numeric.append(float(row['nucfreq_het_frac']))
            low_numeric += int(low)
    return {
        'rows': len(rows),
        'low_baseline_rows': low_count,
        'nucfreq_available_rows': sum(row['nucfreq_available'] == 'True' for row in rows),
        'numeric_nucfreq_rows': len(numeric),
        'numeric_nucfreq_rows_with_low_baseline': low_numeric,
        'maximum_recorded_nucfreq_fraction': max(numeric) if numeric else None,
    }


def _number(text):
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _same_four_decimals(expected, text):
    actual = _number(text)
    if expected is None:
        return actual is None
    return actual is not None and format(expected, '.4f') == format(actual, '.4f')
