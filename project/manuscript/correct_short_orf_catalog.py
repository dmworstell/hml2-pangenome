#!/usr/bin/env python3
"""Apply the caller's existing short-product rule to previously exempt FS-end calls.

Raw sequence and mutation fields are retained. A conservative product-length
upper bound from the recorded residue deletions and insertions establishes that
each corrected product is below 60% without rerunning or changing its alignment.
"""
from pathlib import Path
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
import re

PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / 'results/resolved_manuscript_catalog_20260914/combined_hml2_orf_analysis.RESOLVED.tsv'
OUT = PROJECT / 'results/short_orf_rule_correction_20260915'
DEST = OUT / 'combined_hml2_orf_analysis.RESOLVED.SHORT_ORF_CORRECTED.tsv'
GENES = ('gag', 'pro', 'pol', 'env', 'np9', 'rec')
REFERENCE_LENGTHS = {
    'type1': dict(gag=667, pro=335, pol=874, env=582, np9=75),
    'type2': dict(gag=667, pro=335, pol=957, env=700, rec=105),
}


def product_length_upper_bound(row, gene):
    reference = REFERENCE_LENGTHS[row['provirus_type']][gene]
    deleted = set()
    inserted = 0
    signature = row['missense_' + gene]
    for first, last in re.findall(r'(?:^|,)(\d+)-(\d+)del(?:,|$)', signature):
        deleted.update(range(int(first), int(last) + 1))
    for residue in re.findall(r'(?:^|,)(\d+)del(?:,|$)', signature):
        deleted.add(int(residue))
    for sequence in re.findall(r'\d+ins([A-Z*]+)', signature):
        inserted += len(sequence)
    # Including the reference stop and one additional residue makes this a
    # conservative bound, not a claim about the exact translated length.
    return reference - len(deleted) + inserted + 1, reference


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with SOURCE.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = list(reader.fieldnames)
        rows = list(reader)
    changes = []
    delta = defaultdict(float)
    before_totals = Counter()
    after_totals = Counter()
    weighted_before = Counter()
    weighted_after = Counter()
    prov_before = []
    prov_after = []
    for row in rows:
        public = bool(re.fullmatch(r'(HG|NA)\d+', row['ID']))
        included = public and row['analysis_include'] == '1'
        if included:
            for field in ('expected_intact_gag_copies', 'expected_intact_accessory_copies'):
                before_totals[field] += float(row[field]) if row[field] not in ('', 'NA') else 0
                weighted_before[field] += float(row['observation_weight']) * float(row[field]) if row[field] not in ('', 'NA') else 0
        if included and row['Structure'] in {'Provirus', 'Provirus_from_Multi'}:
            prov_before.append(tuple(row[field] for field in fields))
        for gene in GENES:
            if row[gene] not in {'Intact', 'Intact_FS_End'}:
                continue
            coverage = float(row[gene + '_coverage'])
            if coverage >= 0.6:
                continue
            upper, reference = product_length_upper_bound(row, gene)
            if not upper < reference * 0.6:
                raise ValueError(f'Cannot establish short product: {row["ID_Full"]} {gene}')
            expected = 'expected_intact_accessory_copies' if gene in {'np9', 'rec'} else f'expected_intact_{gene}_copies'
            old_count = row[expected]
            changes.append(dict(Locus=row['Locus'], ID_Full=row['ID_Full'], ID=row['ID'], gene=gene,
                previous_verdict=row[gene], corrected_verdict='Fragment_Intact', reference_coverage=coverage,
                product_length_upper_bound=upper, reference_protein_length=reference,
                previous_expected_copies=old_count, corrected_expected_copies='0', public_included=int(included)))
            row[gene] = 'Fragment_Intact'
            row[expected] = '0'
            if included and old_count not in ('', 'NA'):
                delta[gene] += float(old_count) * float(row['observation_weight'])
        if included:
            for field in ('expected_intact_gag_copies', 'expected_intact_accessory_copies'):
                after_totals[field] += float(row[field]) if row[field] not in ('', 'NA') else 0
                weighted_after[field] += float(row['observation_weight']) * float(row[field]) if row[field] not in ('', 'NA') else 0
        if included and row['Structure'] in {'Provirus', 'Provirus_from_Multi'}:
            prov_after.append(tuple(row[field] for field in fields))
    assert prov_before == prov_after, 'Short-rule correction altered a proviral screening row'
    assert all(x['gene'] in {'gag', 'rec'} for x in changes)
    changed_ids = {x['ID_Full'] for x in changes}
    assert all(row['provirus_type'] == 'type2' for row in rows if row['ID_Full'] in changed_ids)
    with DEST.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    audit = OUT / 'short_orf_classification_changes.tsv'
    with audit.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(changes[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(changes)
    summary = dict(source=str(SOURCE), source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        output=str(DEST), output_sha256=hashlib.sha256(DEST.read_bytes()).hexdigest(),
        affected_records=len(changes), public_changes=dict(Counter(x['gene'] for x in changes if x['public_included'])),
        nonpublic_changes=dict(Counter(x['gene'] for x in changes if not x['public_included'])),
        public_expected_copies_before=dict(before_totals), public_expected_copies_after=dict(after_totals),
        public_weighted_expected_copies_before=dict(weighted_before), public_weighted_expected_copies_after=dict(weighted_after),
        mean_copy_decrement_per_292_donors={gene:n/292 for gene,n in delta.items()},
        proviral_screen_rows_unchanged=len(prov_after), type1_calls_unchanged=True,
        rule='The <60% translated-product length rule applies to Intact and Intact_FS_End. Short products become Fragment_Intact. Raw alignment and mutation fields are unchanged.')
    (OUT / 'verification.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
