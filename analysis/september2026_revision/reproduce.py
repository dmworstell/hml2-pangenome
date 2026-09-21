"""Reproduce September 21 manuscript tables from retained biological inputs."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse, csv, json, math, os, shutil, subprocess, sys

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DATA = REPO/'Supplementary_Data'
GROUPS = {
 'arrays': ('Table_S15', ['arrays_analyze.py','arrays_sites.py','arrays_clock.py'], [
  'Table_S15_pairwise_nucleotide_differences.tsv','Table_S15_array_summary.tsv',
  'Table_S15_variable_sites.tsv','Table_S15_indel_runs.tsv','Table_S15_duplication_order.tsv',
  'Table_S15_locus_summary.tsv','sequence_boundary_validation.tsv','Table_S15_site_frequencies.tsv',
  'Table_S15_conditional_clock_sensitivity.tsv']),
 'solo': ('Table_S16', ['solo_compare.py','solo_figure.py'], ['solo_LTR_diversity_comparison.tsv','solo_LTR_diversity_exclusions.tsv','solo_comparison_LTR_boundaries.tsv','Table_S16B_expected_differences.tsv']),
 'type1': ('Table_S17', ['type1_analyze.py'], ['window_pairwise_differences.tsv','window_summary.tsv','window_callability.tsv','candidate_similarity_to_TypeI.tsv','nearest_type2_by_window.tsv','human_cassette_site_frequencies.tsv','primate_species_counts.tsv','human_full_cassette_6000_7293.aligned.fa','primate_full_cassette_6000_7293.aligned.fa','human_1001bp_cassette_flanks.aligned.fa','human_cassette_consensus.fa']),
 'phylogeny': ('Table_S14/phylogeny_review', ['phylogeny_analyze.py --apd','phylogeny_analyze.py','phylogeny_tables.py'], ['bootstrap_splits.tsv','bootstrap_nearest_loci.tsv','bootstrap_nearest_loci_replicates.tsv','pairwise_differences_by_region.tsv','pairwise_differences_network.tsv','tree_pairwise_distances.tsv','focal_tree_distance_comparisons.tsv','network_sequence_denominators.tsv']),
 'tsd': ('Figure_S7_boundary_evidence', ['tsd_rebuild.py'], ['Figure_S7_TSD_per_locus.tsv','Figure_S7_TSD_population_denominators.tsv','Figure_S7_TSD_sequence_pairs.tsv','Figure_S7_record_level_boundary_calls.tsv','manuscript_terminal_boundary_measurements.tsv','manuscript_TSD_unmatched_records.tsv','comment90_current_boundary_disposition.tsv','comment90_retained_terminal_evidence.tsv','Fiberseq_cell_type_roster.tsv']),
}

def equivalent(actual, expected):
    if actual == expected: return True
    try:
        a, b = float(actual), float(expected)
        return math.isfinite(a) and math.isfinite(b) and math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-14)
    except ValueError:
        return False

def verify(actual, expected):
    if actual.suffix != '.tsv':
        assert actual.read_bytes() == expected.read_bytes(), actual
        return {'file':expected.name, 'comparison':'identical bytes'}
    with actual.open() as f: ar=list(csv.DictReader(f,delimiter='\t'))
    with expected.open() as f: er=list(csv.DictReader(f,delimiter='\t'))
    assert len(ar)==len(er), (expected.name,len(ar),len(er))
    assert not ar or set(ar[0])==set(er[0]), expected.name
    # Row order in original filesystem traversal is not a biological property.
    fields=sorted(ar[0]) if ar else []
    key=lambda r:tuple(r[k] for k in fields)
    ar.sort(key=key);er.sort(key=key)
    for i,(a,e) in enumerate(zip(ar,er)):
        for k in fields:
            assert equivalent(a[k],e[k]),(expected.name,i,k,a[k],e[k])
    return {'file':expected.name,'rows':len(ar),'comparison':'all fields equal; float tolerance 1e-12 relative / 1e-14 absolute'}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,required=True,help='New output directory, outside Supplementary_Data')
    ap.add_argument('--only',choices=list(GROUPS),help='Rerun one analysis group')
    ap.add_argument('--figures',action='store_true',help='Also rebuild the revised figures and 24-page cassette alignment')
    args=ap.parse_args();out=args.out.resolve()
    assert not out.exists(),f'Choose a new output directory: {out}'
    assert not out.is_relative_to(DATA.resolve()),'Never write into the published data'
    out.mkdir(parents=True)
    def run(name):
        rel,commands,targets=GROUPS[name];work=out/name
        shutil.copytree(DATA/rel,work)
        for target in targets: (work/target).unlink()
        env={**os.environ,'HML2_REVISION_OUTPUT':str(work),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MPLCONFIGDIR':str(out/'matplotlib')}
        commands=list(commands)
        if args.figures:
            commands += {'arrays':['arrays_figure.py'], 'type1':['type1_figures.py'], 'phylogeny':['phylogeny_figure.py']}.get(name,[])
        for command in commands:
            words=command.split();log=work/(words[0]+('.apd' if '--apd' in words else '')+'.log')
            with log.open('w') as f:
                subprocess.run([sys.executable,str(HERE/words[0]),*words[1:]],cwd=HERE,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
        checks=[verify(work/t,DATA/rel/t) for t in targets]
        result={'group':name,'status':'PASS','checks':checks}
        (work/'reproduction_checks.json').write_text(json.dumps(result,indent=2)+'\n')
        print(name,'PASS',len(checks),'tables or sequence files',flush=True)
        return result
    names=[args.only] if args.only else list(GROUPS)
    with ThreadPoolExecutor(max_workers=3) as pool: results=list(pool.map(run,names))
    (out/'reproduction_checks.json').write_text(json.dumps({'status':'PASS','python':sys.version,'groups':results},indent=2)+'\n')
    print('PASS',sum(len(r['checks']) for r in results),'published outputs reproduced')
if __name__=='__main__': main()
