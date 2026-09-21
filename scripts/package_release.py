"""Build matching code/data archives from a committed release and prior input archive."""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import argparse
import csv
import gzip
import hashlib
import io
import json
import subprocess


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def table(rows, fields):
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--previous-data', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    version = next(line.split(':', 1)[1].strip() for line in (root/'CITATION.cff').read_text().splitlines() if line.startswith('version:'))
    assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=root).strip(), 'Commit release files before packaging'
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip()
    tracked = sorted(filter(None, subprocess.check_output(['git', 'ls-files', '-z'], cwd=root, text=True).split('\0')))
    args.output.mkdir(parents=True, exist_ok=True)
    code = args.output/f'HML2_code_v{version}.zip'
    subprocess.run(['git', 'archive', '--format=zip', '--output='+str(code.resolve()), commit], cwd=root, check=True)
    catalog = gzip.decompress((root/'data/combined_hml2_orf_analysis.RESOLVED.REC_CORRECTED.tsv.gz').read_bytes())
    catalog_path = 'project/results/rec_exon_boundary_correction_20260915/combined_hml2_orf_analysis.RESOLVED.REC_CORRECTED.tsv'
    override = {n: (root/n).read_bytes() for n in tracked if n.startswith('Supplementary_Data/')}
    override[catalog_path] = catalog
    override['Supplementary_Data/Catalog/HML2_structural_and_ORF_catalog.tsv'] = catalog
    supp_rows = [{'file': n.removeprefix('Supplementary_Data/'), 'bytes': len(raw), 'sha256': digest(raw)}
                 for n, raw in sorted(override.items()) if n.startswith('Supplementary_Data/') and n != 'Supplementary_Data/source_manifest.tsv']
    override['Supplementary_Data/source_manifest.tsv'] = table(supp_rows, ['file', 'bytes', 'sha256'])
    inputs = list(csv.DictReader((root/'external_inputs.tsv').open(), delimiter='\t'))
    for row in inputs:
        if row['path'] in tracked:
            override[row['path']] = (root/row['path']).read_bytes()
    archive = args.output/f'HML2_derived_data_v{version}.zip'
    rows, seen = [], set()
    with ZipFile(args.previous_data) as old, ZipFile(archive, 'w', ZIP_DEFLATED, compresslevel=6) as new:
        assert old.testzip() is None
        def put(name, raw, source):
            assert name not in seen and not name.startswith('/') and '..' not in Path(name).parts
            seen.add(name)
            new.writestr(name, raw)
            rows.append({'path': name, 'bytes': len(raw), 'sha256': digest(raw), 'source': source})
        for info in old.infolist():
            name = info.filename
            if info.is_dir() or name in {'file_manifest.tsv', 'README.txt'} or name.startswith('Supplementary_Data/') or name in override:
                continue
            raw = (root/name).read_bytes() if name in tracked else old.read(name)
            put(name, raw, commit if name in tracked else 'Retained prior public archive input')
        for name, raw in sorted(override.items()):
            put(name, raw, commit)
        put('README.txt', (
            f'HML-2 pangenome code and derived data, version {version}\n\nCode commit: {commit}\n\n'
            f'Extract HML2_code_v{version}.zip and this derived-data ZIP into the same directory. '
            'The code archive contains the new September revision inputs and reproduction commands. '
            'This data archive retains the complete corrected catalog, extracted sequences, '
            'phylogenetic inputs, testing families and phenotype-analysis inputs from the previous public release, '
            'with the current supplementary tables. See analysis/september2026_revision/README.md '
            'for revised tables and figure reproduction.\n\n'
            'file_manifest.tsv records the relative path, size and SHA-256 of every other file. '
            'The manuscript catalog is Supplementary_Data/Catalog/HML2_structural_and_ORF_catalog.tsv. '
            'The matching analysis-path catalog is byte-identical. Original code is MIT licensed. '
            'Original derived data are CC BY 4.0. Third-party sources retain their existing terms.\n'
        ).encode(), 'Release metadata')
        new.writestr('file_manifest.tsv', table(rows, ['path', 'bytes', 'sha256', 'source']))
    with ZipFile(archive) as z:
        assert z.testzip() is None
        for row in rows:
            raw = z.read(row['path'])
            assert len(raw) == row['bytes'] and digest(raw) == row['sha256'], row['path']
        for row in inputs:
            raw = z.read(row['path']) if row['path'] in z.namelist() else (root/row['path']).read_bytes()
            assert digest(raw) == row['sha256'], row['path']
    with ZipFile(code) as z:
        assert z.testzip() is None
        assert sorted(info.filename for info in z.infolist() if not info.is_dir()) == tracked
        for name in tracked:
            assert z.read(name) == (root/name).read_bytes(), name
    summary = {'version': version, 'commit': commit, 'code_files': len(tracked), 'data_files': len(rows),
               'declared_inputs_checked': len(inputs), 'catalog_sha256': digest(catalog),
               'archives': [{'file': p.name, 'bytes': p.stat().st_size, 'sha256': digest(p.read_bytes()),
                             'md5': hashlib.md5(p.read_bytes()).hexdigest()} for p in (code, archive)]}
    (args.output/'verification.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
