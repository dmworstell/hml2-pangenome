"""Extract the current corrected catalog bundled with this code release."""
from pathlib import Path
import gzip
import shutil

ROOT = Path(__file__).resolve().parents[1]
source = ROOT/'data/combined_hml2_orf_analysis.RESOLVED.SHORT_ORF_CORRECTED.tsv.gz'
destination = ROOT/'project/results/short_orf_rule_correction_20260915/combined_hml2_orf_analysis.RESOLVED.SHORT_ORF_CORRECTED.tsv'
destination.parent.mkdir(parents=True, exist_ok=True)
if destination.exists():
    with gzip.open(source, 'rb') as expected, destination.open('rb') as current:
        while True:
            a, b = expected.read(1024 * 1024), current.read(1024 * 1024)
            if a != b:
                raise RuntimeError(f'Refusing to overwrite a different file: {destination}')
            if not a:
                break
else:
    with gzip.open(source, 'rb') as stream, destination.open('xb') as target:
        shutil.copyfileobj(stream, target)
print(destination.relative_to(ROOT))
