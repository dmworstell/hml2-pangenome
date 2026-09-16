"""Reconstruct the corrected public testing families and catalog summaries."""
from collections import Counter
from pathlib import Path
import csv
import gzip
import math
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

def rows(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', newline='') as stream:
        yield from csv.DictReader(stream, delimiter='\t')

def number(value):
    return float(value) if value not in ('', 'NA', 'nan', 'NaN') else math.nan

def check_bh(path, p_field, q_field, expected_rows, expected_finite):
    p, q = [], []
    total = 0
    for row in rows(path):
        total += 1
        pv, qv = number(row[p_field]), number(row[q_field])
        if path.name.startswith('Table_S10g') and row['included_in_bh_family'] != 'True':
            assert math.isnan(qv)
            continue
        if math.isfinite(pv):
            assert 0 <= pv <= 1 and math.isfinite(qv)
            p.append(pv)
            q.append(qv)
        else:
            assert math.isnan(qv)
    assert total == expected_rows and len(p) == expected_finite, (path, total, len(p))
    p, q = np.array(p), np.array(q)
    order = np.argsort(p, kind='stable')
    corrected = np.minimum(1, np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1])
    error = np.max(np.abs(corrected-q[order]))
    assert error < 1e-11, (path, error)
    print(path.name, total, 'rows,', len(p), 'finite P values, max q difference', error)

def main():
    supp = ROOT/'Supplementary_Data'
    check_bh(supp/'Table_S10a_functional_237_model_results.tsv', 'analysis_p_pedigree_cluster', 'q_bh_finite_model_suite', 237, 174)
    check_bh(supp/'Table_S10e_MAGE_complete_discovery_family.tsv.gz', 'p_value', 'bh_global_targeted_family', 1289856, 1289856)
    check_bh(supp/'Table_S10g_GEUVADIS_complete_SLC44A5_followup_family.tsv', 'p_two_sided_hc3', 'q_bh_secondary_family', 98, 30)
    s7 = [r for r in rows(supp/'Table_S7_array_copy_number_distributions.tsv') if r['locus']=='7p22.1']
    called = [r for r in s7 if r['call_status']=='called']
    missing = [r for r in s7 if r['call_status']=='unknown']
    assert sum(int(r['haplotypes']) for r in called)==583
    assert sum(int(r['haplotypes']) for r in missing)==1 and all(r['copy_number']=='' for r in missing)
    assert sum(int(r['haplotypes']) for r in called if int(r['copy_number'])==0)==2
    assert sum(int(r['haplotypes']) for r in called if int(r['copy_number'])>=2)==247
    corrected = Counter()
    weighted = Counter()
    for row in rows(ROOT/'data/combined_hml2_orf_analysis.RESOLVED.SHORT_ORF_CORRECTED.tsv.gz'):
        if row['ID'].startswith(('HG', 'NA')) and row['analysis_include']=='1':
            for gene in ('gag', 'rec'):
                if row[gene]=='Fragment_Intact' and float(row[gene+'_coverage']) < .6:
                    corrected[gene] += 1
            for key in ('expected_intact_gag_copies', 'expected_intact_accessory_copies'):
                weighted[key] += number(row['observation_weight'])*number(row[key]) if row[key] not in ('', 'NA') else 0
            if row['Locus']=='HML-2_8q24.3a' and row['gag']=='Fragment_Intact':
                assert row['expected_intact_gag_copies']=='0'
    assert corrected['gag']>=482 and corrected['rec']>=584
    assert weighted=={'expected_intact_gag_copies':5002.0, 'expected_intact_accessory_copies':17675.0}, weighted
    print('PASS: catalog weighted copy totals, short-product reclassification and 7p22.1 missing calls')

if __name__=='__main__':
    main()
