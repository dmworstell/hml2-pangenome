"""Read-only checks for the compact publication release."""
from pathlib import Path
import ast
from collections import Counter
import csv
import hashlib
import importlib.util
import json
import math
import sys

ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def load_module(name, path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


def main():
    sources=list(ROOT.rglob('*.py'))
    for path in sources:
        if '.git' not in path.parts:
            ast.parse(path.read_text(),filename=str(path))
    tables=read_tsv(ROOT/'Supplementary_Data/source_manifest.tsv')
    assert len(tables)==18
    for row in tables:
        assert hashlib.sha256((ROOT/'Supplementary_Data'/row['file']).read_bytes()).hexdigest()==row['sha256'],row['file']
    for row in read_tsv(ROOT/'data/additional_manifest.tsv'):
        assert hashlib.sha256((ROOT/row['file']).read_bytes()).hexdigest()==row['sha256'],row['file']
    seven=read_tsv(ROOT/'Supplementary_Data/Figure_S7_donor_categories.tsv')
    assert len(seven)==len({r['donor'] for r in seven})==292
    assert Counter(r['category'] for r in seven)=={
        'Four intact ORFs, no Y195C':7,'Four intact ORFs, with Y195C':259,
        'No observed copy meeting either criterion':26}
    eight=read_tsv(ROOT/'Supplementary_Data/Figure_7_solo_LTR_haplotype_counts.tsv')
    assert len(eight)==13 and sum(int(r['observations']) for r in eight)==461
    assert int(eight[0]['observations'])==432
    provirus=json.loads((ROOT/'data/8q11_23_provirus_annotation.json').read_text())
    assert provirus['5_prime_TSD']==provirus['3_prime_TSD']=='CACAC'
    type1=[r for r in read_tsv(ROOT/'Supplementary_Data/Table_S11_direct_TypeI_locus_calls.tsv') if r['locus']!='TOTAL']
    assert len(type1)==20
    assert sum(int(r['callable_typeI']) for r in type1)==9697
    assert sum(int(r['callable_typeII']) for r in type1)==0
    assert sum(int(r['internal_bearing_unresolved']) for r in type1)==36
    cnv=read_tsv(ROOT/'data/cnv_depth_summary.tsv')
    numeric=[r for r in cnv if r['nucfreq_het_frac'] not in ('','NA','nan')]
    assert len(cnv)==416 and len(numeric)==351
    assert max(float(r['nucfreq_het_frac']) for r in numeric)==.0135
    for row in numeric:
        actual=int(row['nucfreq_het_sites'])/int(row['nucfreq_total_sites'])
        assert abs(actual-float(row['nucfreq_het_frac']))<=.00005000001
    helper=load_module('helper',ROOT/'project/manuscript/build_helper_sensitivity.py')
    for name,computed,key in [
        ('helper_equilibrium_requirements.tsv',helper.helper_rows(),'minimum_replication_gain_for_equilibrium'),
        ('neutral_founder_hitting_bounds.tsv',helper.neutral_rows(),'martingale_upper_bound_probability_ever_reach_target')]:
        saved=read_tsv(ROOT/'project/manuscript/delta292_helper_abc_v1'/name)
        assert len(saved)==len(computed)
        for a,b in zip(saved,computed):
            assert math.isclose(float(a[key]),b[key],rel_tol=1e-12,abs_tol=1e-14)
    model=load_module('source_effect',ROOT/'project/manuscript/build_delta292_source_vs_effect_simulation.py')
    assert model.ALPHAS==(.1,.25,.5,1.,2.,5.,20.)
    assert model.N_CLASSES==7 and model.OBSERVATIONS==(('orthology_aware_floor',10),('inclusive_ceiling',17))
    comparisons=read_tsv(ROOT/'Supplementary_Data/Figure_5_source_effect_model_comparison.tsv')
    for row in comparisons:
        alpha=float(row['dirichlet_alpha_per_lesion'])
        assert math.isclose(float(row['source_opportunity_cv']),1/math.sqrt(alpha))
        log_evidence=model.quadrature_log_evidence(alpha,
            float(row['focal_ascertainment_odds_low']),float(row['focal_ascertainment_odds_high']),
            int(row['focal_count']),row['model'])
        assert math.isclose(log_evidence,float(row['log_marginal_likelihood']),rel_tol=1e-11,abs_tol=1e-10),row
    print(f'PASS: Python syntax, {len(tables)} supplementary table hashes, compact denominators, helper formulas, and {len(comparisons)} marginal-evidence calculations')


if __name__=='__main__':
    main()
