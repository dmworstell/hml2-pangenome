#!/usr/bin/env python3
"""Re-call Rec from retained original alignments after fixing its stop boundary."""
from collections import Counter
import csv
from functools import lru_cache
import gzip
import hashlib
import json
import logging
from pathlib import Path
import sys

PROJECT = Path(__file__).resolve().parents[1]
SHARED = PROJECT.parent / 'HML2_ProjectResources/cluster_pipeline_source/shared'
sys.path.insert(0, str(SHARED))
import orf_analysis as caller
from Bio import SeqIO

SOURCE = PROJECT / 'results/short_orf_rule_correction_20260915/combined_hml2_orf_analysis.RESOLVED.SHORT_ORF_CORRECTED.tsv'
OUT = PROJECT / 'results/rec_exon_boundary_correction_20260915'
DEST = OUT / 'combined_hml2_orf_analysis.RESOLVED.REC_CORRECTED.tsv'
REC_FIELDS = ['rec', 'missense_rec', 'rec_dna_identity', 'rec_protein_identity', 'rec_coverage']
COMPATIBLE = {'Intact', 'Intact_FS_End', 'Frameshift_at_end'}
REFERENCE_FILE = PROJECT.parent / 'HML2_ProjectResources/data/ref/type2_KCON.fa'
REFERENCE = str(SeqIO.read(REFERENCE_FILE, 'fasta').seq).upper()
JARGON = {'Deletion': 'Protein_Missing', 'too_short': 'Sequence_Too_Short',
          'invalid_chars': 'Invalid_Characters', 'ref_map_error': 'Mapping_Error',
          'translation_error': 'Translation_Error', 'ComparisonError': 'Reference_Missing'}


@lru_cache(maxsize=None)
def call_rec(ref1, ref2, sample1, sample2, end):
    ref_gapped = ref1 + ref2
    sample_gapped = sample1 + sample2
    ref_dna = REFERENCE[6450:6711] + REFERENCE[8410:end]
    assert ref_gapped.replace('-', '') == ref_dna
    ref_protein, _ = caller.translate_and_check(ref_dna)
    clean_ref, clean_sample = caller.trim_local_gaps(ref_gapped, sample_gapped)
    status, protein, corrected, frameshifts, stop, insertions = caller.analyze_orf_structural_integrity(
        clean_sample, ref_protein, ref_dna, feature_name='rec', locus_type='TypeII',
        gapped_ref_dna=clean_ref)
    coverage = caller.calculate_reference_coverage(ref_gapped, sample_gapped)
    if coverage < caller.COVERAGE_FLOOR_FOR_IDENTITY:
        dna_identity = protein_identity = 'NA'
    else:
        dna_identity = f'{caller.calculate_alignment_identity(sample_gapped.replace("-", ""), ref_dna):.4f}'
        protein_identity = f'{caller.calculate_alignment_identity(protein, ref_protein):.4f}' if protein else '0.0000'
    mutations = caller.align_and_call_mutations(ref_protein, protein, corrected_protein=corrected,
        frameshift_locations=frameshifts, stop_at_aa=stop, large_insertions=insertions)
    mutation_string = ','.join(mutations)
    return dict(zip(REC_FIELDS, [JARGON.get(status, status), JARGON.get(mutation_string, mutation_string),
                                dna_identity, protein_identity, f'{coverage:.4f}']))


def equal_field(field, left, right):
    if field in REC_FIELDS[2:] and left not in ('NA','') and right not in ('NA',''):
        return float(left) == float(right)
    return left == right


def write_tsv(path, rows, fields):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    logging.disable(logging.CRITICAL)
    slices = {}
    for kind in ('cluster', 'local'):
        with gzip.open(OUT / f'{kind}_rec_alignment_slices.jsonl.gz', 'rt') as handle:
            for line in handle:
                record = json.loads(line)
                key = (record['source'], record['record_id'].split('|')[0])
                assert key not in slices, key
                slices[key] = record
    with SOURCE.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = list(reader.fieldnames)
        rows = list(reader)
    # These three retained alignments were regenerated on August 2, after the
    # June completed-call snapshot used by the catalog. Validate their paired
    # August call files, rather than pretend those annotations came from June.
    retained_calls = {}
    for line in (OUT / 'rec_replay_sources.jsonl').read_text().splitlines():
        receipt = json.loads(line)
        path = Path(receipt['path'])
        if path.name.endswith('_orf_integrity_results_type2_KCON.csv'):
            local = OUT / path.name
            assert hashlib.sha256(local.read_bytes()).hexdigest() == receipt['sha256']
            with local.open(newline='') as handle:
                for prior in csv.DictReader(handle):
                    retained_calls[path.parent.name, prior['ID_Full']] = prior
    mismatches, version_differences, changes, evidence = [], [], [], []
    transitions = Counter()
    corrected = []
    for index, row in enumerate(rows):
        if row['rec'] in ('', 'NA', 'nan'):
            continue
        key = row['row_source_path'], row['ID_Full']
        if key not in slices:
            mismatches.append(dict(Locus=row['Locus'], ID_Full=row['ID_Full'], field='alignment', catalog='present', replay='missing'))
            continue
        record = slices[key]
        assert caller.extract_source_id(record['record_id']) == row['Source_Identifier'], key
        ref, sample = record['ref_slices'], record['sample_slices']
        previous = call_rec(ref[0], ref[1], sample[0], sample[1], 8466)
        for field in REC_FIELDS:
            if not equal_field(field, row[field], previous[field]):
                discrepancy = dict(Locus=row['Locus'], ID_Full=row['ID_Full'], field=field, catalog=row[field], replay=previous[field])
                paired = retained_calls.get((Path(row['row_source_path']).parent.name, row['ID_Full']))
                paired_value = JARGON.get(paired[field], paired[field]) if paired else None
                if paired and equal_field(field, paired_value, previous[field]):
                    version_differences.append(discrepancy)
                else:
                    mismatches.append(discrepancy)
        current = call_rec(ref[0], ref[2], sample[0], sample[2], 8467)
        transitions[row['rec'] + ' -> ' + current['rec']] += 1
        corrected.append((index, current))
        evidence.append(dict(Locus=row['Locus'], ID_Full=row['ID_Full'],
            source_identifier=row['Source_Identifier'],
            historical_catalog_source_sha256=row['row_source_sha256'],
            alignment_sha256=record['alignment_sha256'],
            spliced_alignment_sha256=hashlib.sha256((ref[0]+ref[2]+'\n'+sample[0]+sample[2]).encode()).hexdigest()))
    write_tsv(OUT / 'old_rec_replay_discrepancies.tsv', mismatches,
              ['Locus', 'ID_Full', 'field', 'catalog', 'replay'])
    write_tsv(OUT / 'retained_alignment_version_differences.tsv', version_differences,
              ['Locus', 'ID_Full', 'field', 'catalog', 'replay'])
    print('REPLAY', len(corrected), 'rows, discrepancies', len(mismatches), 'unique calls', call_rec.cache_info(), flush=True)
    print('TRANSITIONS', dict(transitions), flush=True)
    if mismatches:
        raise RuntimeError('Original calls did not reproduce. Catalog not written.')
    for index, current in corrected:
        row = rows[index]
        updates = dict(current)
        old_compatible = row['rec'] in COMPATIBLE
        new_compatible = current['rec'] in COMPATIBLE
        if old_compatible != new_compatible:
            updates['expected_intact_accessory_copies'] = row['expected_biological_copy_count'] if new_compatible else '0'
        if (row['rec'] in caller.TRANSLATION_COMPETENT_STATUSES) != (current['rec'] in caller.TRANSLATION_COMPETENT_STATUSES):
            intact = all(row[g] in caller.TRANSLATION_COMPETENT_STATUSES for g in ('gag','pro','env')) and current['rec'] in caller.TRANSLATION_COMPETENT_STATUSES and row['pol_translation_compatible'] != 'N'
            updates['all_orfs_intact'] = 'Y' if intact else 'N'
        if old_compatible != new_compatible:
            intact_count = all(row[g] in COMPATIBLE for g in ('gag','pro','env')) and new_compatible and row['pol_translation_compatible'] != 'N'
            updates['expected_all_orfs_intact_copies'] = row['expected_biological_copy_count'] if intact_count else '0'
        for field, value in updates.items():
            if not equal_field(field, row[field], value):
                changes.append(dict(Locus=row['Locus'], ID_Full=row['ID_Full'], ID=row['ID'], field=field,
                    previous=row[field], corrected=value, analysis_include=row['analysis_include']))
                row[field] = value
    write_tsv(DEST, rows, fields)
    write_tsv(OUT / 'rec_boundary_cell_changes.tsv', changes,
              ['Locus','ID_Full','ID','field','previous','corrected','analysis_include'])
    write_tsv(OUT / 'rec_boundary_alignment_evidence.tsv', evidence,
              ['Locus','ID_Full','source_identifier','historical_catalog_source_sha256','alignment_sha256','spliced_alignment_sha256'])
    summary = dict(source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        output_sha256=hashlib.sha256(DEST.read_bytes()).hexdigest(), rows=len(rows),
        recalled_rows=len(corrected), baseline_discrepancies=len(mismatches),
        older_catalog_annotation_differences=len(version_differences),
        changed_cells=len(changes), changes_by_field=dict(Counter(x['field'] for x in changes)),
        status_transitions=dict(transitions), new_reference_nt=318, new_reference_amino_acids=105,
        reference_translation_characters_including_stop=106)
    (OUT / 'rec_boundary_verification.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
