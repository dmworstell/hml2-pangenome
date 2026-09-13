#!/usr/bin/env python3
"""Materialize the exact direct-call GEUVADIS subset and its audit files."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
WORK = ROOT / "project/working/direct_geuvadis_function_screen_v1"
DATA = WORK / "data"
DATA.mkdir(parents=True, exist_ok=True)

EXPRESSION = DATA / "GEUVADIS_gene_RPKM_50FN_resk10.txt.gz"
SDRF = DATA / "E-GEUV-1.sdrf.txt"
DIRECT_MATRIX = ROOT / "project/working/direct_ebv_fitness_screen_v1/results/person_level_direct_matrix.tsv"
FEATURE_MANIFEST = ROOT / "project/working/direct_ebv_fitness_screen_v1/results/feature_manifest.tsv"
CATALOG = ROOT / "HML2_ProjectResources/data/catalog/combined_hml2_orf_analysis.tsv"
GTF = ROOT / "HML2_ProjectResources/data/ref/gencode/v38/gencode.v38.primary_assembly.annotation.gtf.gz"
MAGE_ALIGNMENT = ROOT / "project/working/onep31b_slc44a5_followup_v1/robust_stats/results/sample_alignment.tsv"


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


for needed in (EXPRESSION, SDRF, DIRECT_MATRIX, FEATURE_MANIFEST, CATALOG, GTF, MAGE_ALIGNMENT):
    if not needed.exists():
        raise FileNotFoundError(needed)

direct_rows = read_tsv(DIRECT_MATRIX)
if len(direct_rows) != 292 or len({row["sample"] for row in direct_rows}) != 292:
    raise AssertionError("Direct matrix must contain 292 unique biological people")
direct_by_id = {row["sample"]: row for row in direct_rows}
manifest = read_tsv(FEATURE_MANIFEST)
if len(manifest) != 79 or len({row["exposure_id"] for row in manifest}) != 79:
    raise AssertionError("Frozen feature manifest must contain 79 unique exposures")

with gzip.open(EXPRESSION, "rt", newline="") as handle:
    reader = csv.reader(handle, delimiter="\t")
    expression_header = next(reader)
expression_ids = expression_header[4:]
if len(expression_ids) != 462 or len(set(expression_ids)) != 462:
    raise AssertionError("Official expression header must contain 462 unique donors")

# The SDRF contains one row per FASTQ mate; require every biological identity to
# have a single consistent sex and exact Source/individual identity.
sdrf_rows = read_tsv(SDRF)
sdrf_by_id = defaultdict(list)
for row in sdrf_rows:
    source = row["Source Name"]
    individual = row["Characteristics[individual]"]
    factor_individual = row["Factor Value[individual]"]
    if source != individual or source != factor_individual:
        raise AssertionError(f"Non-identical GEUVADIS identity labels: {source}, {individual}, {factor_individual}")
    sdrf_by_id[source].append(row)
if len(sdrf_by_id) != 462:
    raise AssertionError(f"Expected 462 unique SDRF people, observed {len(sdrf_by_id)}")

overlap_ids = sorted(set(expression_ids) & set(direct_by_id))
if len(overlap_ids) != 33:
    raise AssertionError(f"Expected exact 33-person direct overlap, observed {len(overlap_ids)}")
if any(sample not in sdrf_by_id for sample in overlap_ids):
    raise AssertionError("Every direct overlap must be present in the official SDRF")

identity_rows = []
for sample in overlap_ids:
    official_sexes = {row["Characteristics[sex]"].lower() for row in sdrf_by_id[sample]}
    official_ancestries = {row["Factor Value[ancestry category]"] for row in sdrf_by_id[sample]}
    if len(official_sexes) != 1 or len(official_ancestries) != 1:
        raise AssertionError(f"Inconsistent SDRF metadata for {sample}")
    local = direct_by_id[sample]
    official_sex = next(iter(official_sexes))
    if official_sex != local["sex"].lower():
        raise AssertionError(f"Sex mismatch for {sample}: {official_sex} versus {local['sex']}")
    identity_rows.append(
        {
            "sample": sample,
            "expression_header_exact": "TRUE",
            "sdrf_source_exact": "TRUE",
            "sdrf_individual_exact": "TRUE",
            "official_sex": official_sex,
            "direct_metadata_sex": local["sex"],
            "official_ancestry_label": next(iter(official_ancestries)),
            "population": local["population"],
            "superpopulation": local["superpopulation"],
            "pedigree_component": local["pedigree_component"],
            "sdrf_rows": len(sdrf_by_id[sample]),
        }
    )
write_tsv(DATA / "identity_audit.tsv", list(identity_rows[0]), identity_rows)

mage_ids = {row["sample"] for row in read_tsv(MAGE_ALIGNMENT)}
mage_overlap_rows = [
    {"sample": sample, "in_prior_MAGE_direct39": "TRUE" if sample in mage_ids else "FALSE"}
    for sample in overlap_ids
]
write_tsv(DATA / "mage_overlap_audit.tsv", ["sample", "in_prior_MAGE_direct39"], mage_overlap_rows)

# Parse current GENCODE gene symbols and GRCh38 TSS coordinates. GEUVADIS uses
# an older Ensembl version, so the stable base ID is the identity key.
gene_annotation = {}
attribute_pattern = re.compile(r'(\w+) "([^"]+)"')
with gzip.open(GTF, "rt") as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9 or fields[2] != "gene":
            continue
        attributes = dict(attribute_pattern.findall(fields[8]))
        gene_id = attributes.get("gene_id", "")
        if not gene_id:
            continue
        base = gene_id.split(".", 1)[0]
        start, end = int(fields[3]), int(fields[4])
        gene_annotation[base] = {
            "gene_symbol": attributes.get("gene_name", ""),
            "chromosome_grch38": fields[0],
            "tss_grch38": start if fields[6] == "+" else end,
        }

selected_index = [expression_header.index(sample) for sample in overlap_ids]
expression_rows = []
with gzip.open(EXPRESSION, "rt", newline="") as handle:
    reader = csv.reader(handle, delimiter="\t")
    header = next(reader)
    if header != expression_header:
        raise AssertionError("Expression header changed between audit and extraction")
    for row in reader:
        if len(row) != len(header):
            raise AssertionError(f"Malformed expression row {row[0] if row else 'EMPTY'}")
        gene_base = row[0].split(".", 1)[0]
        annotation = gene_annotation.get(gene_base, {})
        expression_rows.append(
            {
                "gene_id_source": row[0],
                "gene_id_base": gene_base,
                "gene_symbol": annotation.get("gene_symbol", ""),
                "chromosome_grch38": annotation.get("chromosome_grch38", ""),
                "tss_grch38": annotation.get("tss_grch38", ""),
                "source_chromosome": row[2],
                "source_coordinate": row[3],
                **{sample: row[index] for sample, index in zip(overlap_ids, selected_index)},
            }
        )
if len(expression_rows) != 23722:
    raise AssertionError(f"Expected 23,722 GEUVADIS genes, observed {len(expression_rows)}")
if sum(row["gene_id_base"] == "ENSG00000137968" for row in expression_rows) != 1:
    raise AssertionError("SLC44A5 stable Ensembl ID must occur exactly once")
write_tsv(
    DATA / "direct33_expression.tsv",
    ["gene_id_source", "gene_id_base", "gene_symbol", "chromosome_grch38", "tss_grch38",
     "source_chromosome", "source_coordinate", *overlap_ids],
    expression_rows,
)

# Freeze the complete 79-exposure direct matrix plus every declared adjustment.
exposure_ids = [row["exposure_id"] for row in manifest]
adjustments = sorted(
    {row["adjustment"] for row in manifest if row["adjustment"] not in {"", "NA"}}
    - set(exposure_ids)
)
metadata = ["sample", "sex", "population", "superpopulation", "pedigree_component"]
direct33_rows = [{key: direct_by_id[sample].get(key, "NA") for key in [*metadata, *exposure_ids, *adjustments]}
                 for sample in overlap_ids]
write_tsv(DATA / "direct33_exposures.tsv", [*metadata, *exposure_ids, *adjustments], direct33_rows)
write_tsv(DATA / "feature_manifest.tsv", list(manifest[0]), manifest)

# Locus-resolved Type-I physical units support the Type-I-shared versus 1p31.1b-
# specific comparison. No catalog row is a biological zero, never an unknown.
catalog_rows = read_tsv(CATALOG)
provirus = [row for row in catalog_rows if row["Structure"] in {"Provirus", "Provirus_from_Multi"}
            and row["provirus_type"] in {"type1", "type2"}]
locus_types = defaultdict(set)
for row in provirus:
    locus_types[row["Locus"]].add(row["provirus_type"])
if any(len(types) != 1 for types in locus_types.values()):
    raise AssertionError("Every direct-catalog provirus locus must have one fixed type")
type1_loci = sorted(locus for locus, types in locus_types.items() if types == {"type1"})
if len(type1_loci) != 17:
    raise AssertionError(f"Expected 17 fixed Type-I loci, observed {len(type1_loci)}")
type1_counts = defaultdict(int)
for row in provirus:
    if row["ID"] in overlap_ids and row["Locus"] in type1_loci:
        type1_counts[(row["ID"], row["Locus"])] += 1
type1_rows = []
for sample in overlap_ids:
    type1_rows.append({"sample": sample, **{locus: type1_counts[(sample, locus)] for locus in type1_loci}})
write_tsv(DATA / "type1_locus_units.tsv", ["sample", *type1_loci], type1_rows)

# GRCh38 locus intervals for cis routing. Reference-assembly catalog rows are
# authoritative when present; exact focal-array coordinates are preserved from
# the prior direct cis-expression plan.
interval_pattern = re.compile(r"GRCh38#0#(chr[^:]+):(\d+)-(\d+)$")
intervals = defaultdict(list)
for row in catalog_rows:
    match = interval_pattern.match(row["Source_Identifier"])
    if row["ID"] == "GCA" and match:
        intervals[row["Locus"]].append((match.group(1), int(match.group(2)), int(match.group(3))))
coordinate_rows = []
for locus, values in sorted(intervals.items()):
    chromosomes = {value[0] for value in values}
    if len(chromosomes) != 1:
        raise AssertionError(f"Multiple GRCh38 chromosomes for {locus}")
    coordinate_rows.append({"locus": locus, "chromosome": next(iter(chromosomes)),
                            "start": min(value[1] for value in values),
                            "end": max(value[2] for value in values),
                            "coordinate_authority": "catalog_GRCh38_reference"})
by_locus = {row["locus"]: row for row in coordinate_rows}
for locus, chromosome, start, end in (
    ("HML-2_1p31.1b", "chr1", 75376492, 75384051),
    ("HML-2_7p22.1", "chr7", 4582432, 4600393),
):
    by_locus[locus] = {"locus": locus, "chromosome": chromosome, "start": start, "end": end,
                       "coordinate_authority": "preserved_frozen_GRCh38_direct_interval"}
write_tsv(DATA / "locus_coordinates_grch38.tsv",
          ["locus", "chromosome", "start", "end", "coordinate_authority"],
          [by_locus[locus] for locus in sorted(by_locus)])

summary = {
    "expression_people": len(expression_ids),
    "sdrf_people": len(sdrf_by_id),
    "direct_people": len(direct_rows),
    "exact_direct_overlap": len(overlap_ids),
    "expression_genes": len(expression_rows),
    "gencode_mapped_genes": sum(bool(row["chromosome_grch38"]) for row in expression_rows),
    "populations": sorted({direct_by_id[sample]["population"] for sample in overlap_ids}),
    "superpopulations": sorted({direct_by_id[sample]["superpopulation"] for sample in overlap_ids}),
    "females": sum(direct_by_id[sample]["sex"] == "female" for sample in overlap_ids),
    "males": sum(direct_by_id[sample]["sex"] == "male" for sample in overlap_ids),
    "repeated_pedigree_components": len(overlap_ids) - len({direct_by_id[sample]["pedigree_component"] for sample in overlap_ids}),
    "prior_MAGE_overlap_people": sum(sample in mage_ids for sample in overlap_ids),
    "GEUVADIS_only_people": sum(sample not in mage_ids for sample in overlap_ids),
    "frozen_exposures": len(manifest),
    "fixed_type1_loci": len(type1_loci),
}
(DATA / "materialization_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

source_rows = []
for role, path in (
    ("official_peer50_expression", EXPRESSION),
    ("official_sample_manifest", SDRF),
    ("direct_person_exposure_authority", DIRECT_MATRIX),
    ("frozen_exposure_manifest", FEATURE_MANIFEST),
    ("authoritative_orf_catalog", CATALOG),
    ("gencode_v38_annotation", GTF),
    ("prior_MAGE_direct39_alignment", MAGE_ALIGNMENT),
):
    source_rows.append({"role": role, "path": str(path.relative_to(ROOT)), "sha256": sha256(path)})
write_tsv(DATA / "source_manifest.tsv", ["role", "path", "sha256"], source_rows)

print(json.dumps(summary, indent=2))
