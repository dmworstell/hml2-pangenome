"""Replay every Rec-bearing record from the retained alignment slices."""
from pathlib import Path
import gzip
import hashlib
import importlib.util
import json
import shutil
import tempfile

ROOT=Path(__file__).resolve().parents[1]
source=ROOT/'project/manuscript/correct_rec_exon_boundary.py'
spec=importlib.util.spec_from_file_location('rec_boundary',source)
owner=importlib.util.module_from_spec(spec)
spec.loader.exec_module(owner)
retained=ROOT/'project/results/rec_exon_boundary_correction_20260915'
with tempfile.TemporaryDirectory(prefix='hml2-rec-replay-') as tmp:
    output=Path(tmp)
    for path in retained.iterdir():
        if path.name.endswith(('_rec_alignment_slices.jsonl.gz','_orf_integrity_results_type2_KCON.csv')) or path.name=='rec_replay_sources.jsonl':
            shutil.copy2(path,output/path.name)
    owner.OUT=output
    owner.DEST=output/'combined_hml2_orf_analysis.RESOLVED.REC_CORRECTED.tsv'
    owner.main()
    report=json.loads((output/'rec_boundary_verification.json').read_text())
    assert report['recalled_rows']==36073 and report['baseline_discrepancies']==0
    expected=hashlib.sha256()
    with gzip.open(ROOT/'data/combined_hml2_orf_analysis.RESOLVED.REC_CORRECTED.tsv.gz','rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):expected.update(chunk)
    assert report['output_sha256']==expected.hexdigest(), 'Rec replay differs from published catalog'
    print('PASS: all 36,073 Rec-bearing records replay and final catalog bytes match the public distribution')
