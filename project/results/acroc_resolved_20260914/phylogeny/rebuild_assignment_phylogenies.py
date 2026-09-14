#!/usr/bin/env python3
"""Build assignment-consistent representative trees from retained KCON panels.

All outputs remain in this script's directory. The manuscript owner supplies
the LTR selection, distance method, subfamily authority, and page renderer.
The source-disconnected historical ORF trees are retained only for comparison.
"""

from __future__ import annotations

import argparse
import csv
import inspect
import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from hashlib import sha256
from pathlib import Path
from functools import lru_cache
from Bio.Seq import Seq

OUT = Path(__file__).resolve().parent
PROJECT = OUT.parents[2]
CATALOG = PROJECT / "results/resolved_manuscript_catalog_20260914/combined_hml2_orf_analysis.RESOLVED.tsv"
BASELINE = PROJECT / "results/biological_orf_annotation_20260802/combined_hml2_orf_analysis.CNV_WEIGHTED.BIOLOGICALLY_ANNOTATED.v3.tsv"
ORIGINAL_PHY = PROJECT / "manuscript/figures/narrative/phylogeny"
ACROCENTRIC = {"13p13", "15p13a", "15p13b", "21p13", "22p13"}
# Positive published family identity, authenticated by the retained official
# sequence panel in fragment_subfamily_verification.json; not missing labels.
NON_HML2_LOCI = {"8p22", "17p13.1_hg38"}
# Zero-based, half-open KCON gene spans. The historical producer's end+1
# appended one nucleotide after each reference stop codon; do not repeat it.
ORF_COORDS = {"gag": (1111, 3112), "pro": (2913, 3918), "pol": (3878, 6749), "env": (6450, 8550)}
MISSING_GENE = {"", "NA", "N", "Protein_Missing", "ComparisonError_ProteinMissing"}
sys.dont_write_bytecode = True
sys.path.insert(0, str(PROJECT / "manuscript"))
os.environ["MPLCONFIGDIR"] = str(OUT / "mplcache")
import build_narrative_main_figures as owner


def read_catalog(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def retained_catalog_rows(rows):
    return [row for row in rows if row["analysis_include"] == "1" and owner.PUBLIC_ID.fullmatch(row["ID"])]


def canonical_distance_matrix(names, sequences):
    array = owner.np.array([list(sequence) for sequence in sequences])
    canonical = owner.np.isin(array, list("ACGT"))
    matrix = owner.np.zeros((len(names), len(names)), dtype=float)
    for i in range(len(names)):
        for j in range(i):
            valid = canonical[i] & canonical[j]
            matrix[i, j] = matrix[j, i] = float(owner.np.mean(array[i, valid] != array[j, valid])) if owner.np.any(valid) else 1.0
    return matrix


def write_region_tree(gene, per_locus):
    names, sequences, selected = [], [], []
    for locus in sorted(per_locus):
        for cluster, (sequence, count) in enumerate(per_locus[locus].most_common(2), start=1):
            name = f"{locus}__hap{cluster}__n{count}"
            names.append(name)
            sequences.append(sequence)
            selected.append({"gene": gene, "locus": locus, "tip": name, "cluster_observations": count,
                             "sequence_sha256": sha256(sequence.encode()).hexdigest(), "sequence": sequence})
    matrix = canonical_distance_matrix(names, sequences)
    # Report the historical >0.30 nearest p-distance flag, but do not discard
    # biological observations under an unvalidated artifact cutoff.
    nearest = matrix + owner.np.eye(len(names)) * 9.0
    diagnostics = [dict(selected[i], nearest_p_distance=float(nearest[i].min()),
                        exceeds_historical_0_30_flag=int(nearest[i].min() > 0.30))
                   for i in range(len(names))]
    lower = [list(matrix[i, :i + 1]) for i in range(len(names))]
    tree = owner.DistanceTreeConstructor().nj(owner.DistanceMatrix(names, lower))
    for clade in tree.find_clades():
        if clade.branch_length is not None and clade.branch_length < 0:
            clade.branch_length = 0.0
    filename = "hml2_pan_ltr_expanded_tree.nwk" if gene == "LTR" else f"hml2_pan_orf_{gene}_tree.nwk"
    owner.Phylo.write(tree, OUT / filename, "newick", format_branch_length="%1.15g")
    counts = Counter({locus: sum(values.values()) for locus, values in per_locus.items()})
    return filename, counts, selected, diagnostics


def alignment_metrics(line, query, source):
    """Replay a base-resolved alignment, including every mismatch and indel."""
    fields = line.split("\t")
    qs, qe, ts, te = map(int, (fields[2], fields[3], fields[7], fields[8]))
    assert int(fields[1]) == len(query) and int(fields[6]) == len(source)
    tags = {f.split(":", 2)[0]: f.split(":", 2)[2] for f in fields[12:]}
    cigar = tags["cg"]
    operations = [(int(n), op) for n, op in re.findall(r"(\d+)([=XID])", cigar)]
    assert "".join(f"{n}{op}" for n, op in operations) == cigar, cigar
    aligned_query = query[qs:qe]
    if fields[4] == "-":
        aligned_query = str(Seq(aligned_query).reverse_complement())
    elif fields[4] != "+":
        raise ValueError("Unknown source-alignment strand")
    q, t = 0, ts
    counts = Counter()
    for n, op in operations:
        if op in "=X":
            actual = sum(a != b for a, b in zip(aligned_query[q:q+n], source[t:t+n]))
            assert len(aligned_query[q:q+n]) == len(source[t:t+n]) == n
            assert actual == (0 if op == "=" else n), (op, actual, n)
            q += n
            t += n
        elif op == "I":
            q += n
        else:
            t += n
        counts[op] += n
    assert q == qe-qs and t == te
    return dict(query_start=qs, query_end=qe, query_bases=len(query),
        source_start=ts, source_end=te, strand=fields[4], cigar=cigar,
        matching_bases=counts["="], mismatching_bases=counts["X"],
        query_inserted_bases=counts["I"], omitted_source_bases=counts["D"],
        exact=int(qs == 0 and qe == len(query) and counts["X"] == 0 and counts["I"] == 0))


@lru_cache(maxsize=8192)
def source_alignment(query, source_path, source_file_sha256):
    """Require a full-query zero-mismatch alignment, not an arbitrary subsequence."""
    source_path = Path(source_path)
    assert sha256(source_path.read_bytes()).hexdigest() == source_file_sha256
    source = str(owner.SeqIO.read(source_path, "fasta").seq).upper()
    command = ["minimap2", "-x", "asm20", "-c", "--eqx", "-N", "20",
               "--secondary=yes", str(source_path), "-"]
    run = subprocess.run(command, input=">retained_projection\n"+query+"\n",
                         text=True, capture_output=True, check=True)
    hits = [alignment_metrics(line, query, source) for line in run.stdout.splitlines()]
    return hits


def bind_source(sequence, candidates, require_sequence, raw_cache, checks, alignment, record):
    if len(candidates) == 1 and not require_sequence:
        return candidates[0], "unique_original_catalog_slot"
    query = sequence.replace("-", "")
    matches = []
    for candidate in candidates:
        identifier = candidate["ID_Full"]
        if identifier not in raw_cache:
            local_record = Path(candidate["source_record_path_v3"])
            processed = owner.PROCESSED_LOCI / candidate["orig_Locus"] / f"{identifier}.fa"
            paths = [path for path in (local_record, processed)
                     if path.suffix.lower() in {".fa", ".fasta", ".fna"} and path.is_file()]
            source = None
            if paths:
                path = paths[0]
                records = list(owner.SeqIO.parse(path, "fasta"))
                if len(records) == 1:
                    source = (str(records[0].seq).upper(), str(records[0].seq.reverse_complement()).upper(),
                              str(path), sha256(path.read_bytes()).hexdigest())
            raw_cache[identifier] = source
        source = raw_cache[identifier]
        hits = source_alignment(query, source[2], source[3]) if source else []
        exact_positions = {(h["source_start"], h["source_end"], h["strand"]) for h in hits if h["exact"]}
        match = len(exact_positions) == 1
        checks.append({"alignment": alignment, "record": record, "candidate_ID_Full": identifier,
            "Source_Identifier": candidate["Source_Identifier"], "physical_locus": candidate["Locus"],
            "source_path": source[2] if source else "", "source_fasta_sha256": source[3] if source else "",
            "exact_alignment_projection": int(match), "source_available": int(source is not None),
            "source_alignment_metrics": json.dumps(hits)})
        if match:
            matches.append(candidate)
    return (matches[0], "unique_full_query_zero_mismatch_source_alignment") if len(matches) == 1 else (None, f"source_identity_ambiguous_{len(matches)}_exact_matches_of_{len(candidates)}_candidates")


def write_source_admission_partition(rows, admissions, exclusions, source_checks):
    admitted = {(row["sample"], row["haplotype"], row["Source_Identifier"]): row for row in admissions}
    slots = Counter()
    for row in [*admissions, *exclusions]:
        match = re.fullmatch(r"(.+)_(hap1|hap2|h1|h2|pat|mat)(#c\d+)?", row["record"])
        if match:
            sample, haplotype, _ = match.groups()
            haplotype = {"hap1": "h1", "hap2": "h2"}.get(haplotype, haplotype)
            slots[row["alignment"].removesuffix("_carrier_kcon.aln.fa"), sample, haplotype] += 1
    checks = defaultdict(list)
    for row in source_checks:
        checks[row["candidate_ID_Full"]].append(row)
    partition = []
    for row in rows:
        if row["observation_state"] != "PRESENT":
            continue
        original = row["orig_Locus"].removeprefix("HML-2_")
        locus = row["Locus"].removeprefix("HML-2_")
        key = row["ID"], row["Haplotype"], row["Source_Identifier"]
        candidate_count = slots[original, row["ID"], row["Haplotype"]]
        source_attempts = checks[row["ID_Full"]]
        if key in admitted:
            status, reason = "admitted", admitted[key]["source_binding_basis"]
        elif locus in NON_HML2_LOCI:
            status, reason = "excluded_biological_family_HML11", "Sequence-authenticated published HML-11 element, not HML-2"
        elif "unresolved" in locus.lower() or "pooled" in locus.lower() or not re.fullmatch(r"(?:[1-9]|1[0-9]|2[0-2]|X|Y)[pq]\d.*", locus):
            status, reason = "excluded_physical_assignment_unresolved", "No named physical locus; not an LTR-subfamily ambiguity"
        elif not candidate_count:
            status, reason = "no_retained_original_slot_alignment", "No retained alignment with this original extraction locus, public sample, and haplotype"
        elif not any(int(attempt["source_available"]) for attempt in source_attempts):
            status, reason = "source_identity_not_established", "Multiple/copy-suffixed candidate identity requires unavailable current source FASTA"
        elif any(int(attempt["exact_alignment_projection"]) for attempt in source_attempts):
            status, reason = "source_identity_not_established", "Retained sequence matches more than one current source; no unique copy identity"
        else:
            status, reason = "source_identity_not_established", "Available retained sequence does not exactly project to this current candidate source"
        partition.append({"ID_Full": row["ID_Full"], "ID": row["ID"], "Haplotype": row["Haplotype"],
            "Source_Identifier": row["Source_Identifier"], "original_locus": original, "physical_locus": locus,
            "retained_alignment_candidates_at_original_slot": candidate_count,
            "partition": status, "reason": reason,
            "admitted_alignment": admitted[key]["alignment"] if key in admitted else "",
            "admitted_record": admitted[key]["record"] if key in admitted else ""})
    write_tsv(OUT / "catalog_source_admission_partition.tsv", partition, list(partition[0]))
    cross = [row for row in partition if row["original_locus"] != row["physical_locus"] and row["physical_locus"] in ACROCENTRIC]
    write_tsv(OUT / "cross_locus_catalog_source_coverage.tsv", cross, list(partition[0]))
    unique_sources = {}
    for row in partition:
        key = row["ID"], row["Haplotype"], row["Source_Identifier"]
        if key in unique_sources:
            assert unique_sources[key] == row["partition"], key
        unique_sources[key] = row["partition"]
    assert sum(value == "admitted" for value in unique_sources.values()) == len(admissions)
    return {"included_public_PRESENT_catalog_rows": len(partition), "unique_current_source_copies": len(unique_sources),
        "catalog_row_partition": dict(Counter(row["partition"] for row in partition)),
        "unique_source_copy_partition": dict(Counter(unique_sources.values())),
        "named_cross_locus_public_catalog_rows": len(cross),
        "named_cross_locus_source_coverage": dict(Counter(row["partition"] for row in cross))}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", type=Path, default=CATALOG, help="Final combined catalog, including the authenticated HML-11 exclusions required by the primary reader.")
    parser.add_argument("--check-inputs", action="store_true", help="Exercise the real catalog-to-alignment binding without publishing outputs.")
    parser.add_argument("--additional-projections", type=Path, required=True,
                        help="Current source-bound KCON projections for recovered graph records")
    args = parser.parse_args()
    catalog = args.catalog.resolve()
    catalog_digest = sha256(catalog.read_bytes()).hexdigest()
    baseline_all, repaired_all = read_catalog(BASELINE), read_catalog(catalog)
    assert len(baseline_all) == len(repaired_all)
    assert set(baseline_all[0]) == set(repaired_all[0]), "Unreviewed catalog schema change"
    family_path = OUT.parent / "HML11_comparator_catalog_rows.tsv"
    family_rows = read_catalog(family_path)
    family_by_id = {row["ID_Full"]: row for row in family_rows}
    assert len(family_by_id) == len(family_rows) == 1176
    assert {row["orig_Locus"].removeprefix("HML-2_") for row in family_rows} == NON_HML2_LOCI
    assert sum(row["public_cohort"] == "1" for row in family_rows) == 1168
    family_reason = "non_HML2_HML11_sequence_identity"
    recovery_path = catalog.parent / "source_recovery_changes.tsv"
    recoveries = read_catalog(recovery_path) if recovery_path.is_file() else []
    recovery_by_id = {row["original_id"]: row for row in recoveries}
    assert len(recovery_by_id) == len(recoveries)
    renamed = []
    family_changes = set()
    changed_columns = Counter()
    changed_ids = set()
    for slot, (before, after) in enumerate(zip(baseline_all, repaired_all)):
        assert all(before[field] == after[field] for field in ("ID", "Haplotype", "orig_Locus")), slot
        allowed_changes = {"Locus", "physical_locus_assignment", "physical_assignment_evidence"}
        recovery = recovery_by_id.get(before["ID_Full"])
        if recovery is not None:
            assert recovery["current_id"] == after["ID_Full"]
            assert recovery["locus"] == after["orig_Locus"]
            assert recovery["observation_state"] == after["observation_state"]
            assert recovery["after"] == after["Structure"]
            assert recovery["sequence_sha256"] == after["observation_sequence_sha256"]
            allowed_changes.update(recovery["changed_fields"].split(","))
        if before["ID_Full"] in family_by_id:
            comparator = family_by_id[before["ID_Full"]]
            assert all(before[field] == after[field] == comparator[field]
                       for field in ("ID_Full", "ID", "Haplotype", "Source_Identifier", "orig_Locus", "Locus"))
            assert comparator["verified_biological_family"] == "HML-11"
            assert comparator["publication_figure_panel"] == "1E"
            assert comparator["published_locus"] == after["orig_Locus"].removeprefix("HML-2_").removesuffix("_hg38")
            assert before["analysis_include"] == "1" and after["analysis_include"] == "0"
            assert after["analysis_exclusion_reason"] == family_reason
            allowed_changes.update({"analysis_include", "analysis_exclusion_reason"})
            family_changes.add(before["ID_Full"])
        if before["ID_Full"] != after["ID_Full"]:
            assert recovery is not None, before["ID_Full"]
            renamed.append({"original_zero_based_slot": slot, "original_ID_Full": before["ID_Full"],
                            "current_ID_Full": after["ID_Full"], "ID": after["ID"],
                            "Haplotype": after["Haplotype"], "orig_Locus": after["orig_Locus"]})
        for column in before:
            if before[column] != after[column]:
                assert column in allowed_changes, (before["ID_Full"], column)
                changed_columns[column] += 1
                changed_ids.add(after["ID_Full"])
    assert family_changes == set(family_by_id)
    baseline = retained_catalog_rows(baseline_all)
    owner.CATALOG = catalog
    rows, roster = owner.load_catalog()
    assert {row["ID_Full"] for row in rows} <= {row["ID_Full"] for row in retained_catalog_rows(repaired_all)}
    assert len({(row["ID"], row["Haplotype"]) for row in rows}) == 584
    kcon_path = PROJECT / "inputs/references/type2_KCON.fa"
    reference = str(owner.SeqIO.read(kcon_path, "fasta").seq).upper()
    coordinate_evidence = {}
    for gene, (start, end) in ORF_COORDS.items():
        coding = reference[start:end]
        assert len(coding) % 3 == 0 and coding[-3:] in {"TAA", "TAG", "TGA"}, gene
        coordinate_evidence[gene] = {"zero_based_start_inclusive": start, "zero_based_end_exclusive": end,
            "reference_nucleotides": len(coding), "terminal_stop": coding[-3:],
            "legacy_extra_nucleotide": reference[end], "corrected_sequence_sha256": sha256(coding.encode()).hexdigest()}
    old_allowed = {(row["Locus"].removeprefix("HML-2_"), row["ID"], row["Haplotype"])
                   for row in baseline if row["observation_state"] == "PRESENT"}
    source_rows = defaultdict(list)
    for row in rows:
        if row["observation_state"] == "PRESENT":
            source_rows[row["orig_Locus"].removeprefix("HML-2_"), row["ID"], row["Haplotype"]].append(row)
    per_locus = {region: defaultdict(Counter) for region in ("LTR", *ORF_COORDS)}
    old_ltr = defaultdict(Counter)
    observation_rows, inputs, admissions, exclusions = [], [], [], []
    bound_records, source_checks, raw_cache = [], [], {}
    total_alignment_records = 0
    for path in sorted((PROJECT / "data").glob("*_carrier_kcon.aln.fa")):
        original_locus = path.name.removesuffix("_carrier_kcon.aln.fa")
        inputs.append({"path": str(path), "sha256": sha256(path.read_bytes()).hexdigest()})
        for record in owner.SeqIO.parse(path, "fasta"):
            total_alignment_records += 1
            match = re.fullmatch(r"(.+)_(hap1|hap2|h1|h2|pat|mat)(#c\d+)?", record.id)
            if not match:
                exclusions.append({"alignment": path.name, "record": record.id, "reason": "unrecognized_or_reference_alignment_record"})
                continue
            sample, haplotype, copy_suffix = match.groups()
            haplotype = {"hap1": "h1", "hap2": "h2"}.get(haplotype, haplotype)
            key = original_locus, sample, haplotype
            sequence = str(record.seq).upper()
            ltr = owner.kcon_ltr_sequence(sequence)
            if not copy_suffix and key in old_allowed and ltr is not None:
                old_ltr[original_locus][ltr] += 1
            if original_locus in NON_HML2_LOCI:
                exclusions.append({"alignment": path.name, "record": record.id, "reason": "biological_family_HML11_not_HML2"})
                continue
            candidates = [row for row in source_rows.get(key, [])
                          if row.get("v3_row_origin") != "BROAD_GRAPH_SOURCE_RECOVERY_20260914"]
            if not candidates:
                exclusions.append({"alignment": path.name, "record": record.id, "reason": "no_included_public_present_original_source_slot"})
                continue
            require_sequence = bool(copy_suffix) or any(row["Locus"].removeprefix("HML-2_") in ACROCENTRIC for row in candidates)
            candidate, basis = bind_source(sequence, candidates, require_sequence, raw_cache,
                                            source_checks, path.name, record.id)
            if candidate is None:
                exclusions.append({"alignment": path.name, "record": record.id, "reason": basis})
                continue
            locus = candidate["Locus"].removeprefix("HML-2_")
            if locus in NON_HML2_LOCI:
                exclusions.append({"alignment": path.name, "record": record.id, "reason": "biological_family_HML11_not_HML2"})
                continue
            if "unresolved" in locus.lower() or "pooled" in locus.lower() or not re.fullmatch(r"(?:[1-9]|1[0-9]|2[0-2]|X|Y)[pq]\d.*", locus):
                exclusions.append({"alignment": path.name, "record": record.id, "reason": "physical_locus_remains_unresolved_or_pooled"})
                continue
            newly_admitted = bool(copy_suffix) or key not in old_allowed
            bound_records.append({"alignment": path.name, "record": record.id, "original_locus": original_locus, "locus": locus,
                "sample": sample, "haplotype": haplotype, "newly_admitted": int(newly_admitted),
                "candidate_source_count": len(candidates),
                "ID_Full": candidate["ID_Full"], "Source_Identifier": candidate["Source_Identifier"],
                "source_binding_basis": basis, "copy_suffixed_record": int(bool(copy_suffix)),
                "physical_locus_changed_from_extraction": int(locus != original_locus),
                "aligned_canonical_bases": sum(base in "ACGT" for base in sequence),
                "aligned_sequence_sha256": sha256(sequence.encode()).hexdigest(), "sequence": sequence})
    current_by_id = {row["ID_Full"]: row for row in rows}
    additional = read_catalog(args.additional_projections)
    expected_additional = {row["ID_Full"] for row in rows
        if row.get("v3_row_origin") == "BROAD_GRAPH_SOURCE_RECOVERY_20260914"}
    assert {row["ID_Full"] for row in additional} == expected_additional
    assert len(additional) == len(expected_additional)
    for projected in additional:
        candidate = current_by_id[projected["ID_Full"]]
        assert all(projected[key] == candidate[key] for key in
                   ("ID", "Haplotype", "Locus", "orig_Locus", "Source_Identifier"))
        assert sha256(Path(projected["source_fasta"]).read_bytes()).hexdigest() == projected["source_file_sha256"]
        assert sha256(Path(projected["alignment_path"]).read_bytes()).hexdigest() == projected["alignment_sha256"]
        sequence = projected["typeII_KCON_sequence"]
        assert len(sequence) == len(reference) == 9472
        assert sha256(sequence.encode()).hexdigest() == projected["typeII_KCON_sequence_sha256"]
        locus = candidate["Locus"].removeprefix("HML-2_")
        assert locus not in NON_HML2_LOCI and locus not in {"acro_type1", "acro_type2"}
        bound_records.append(dict(alignment=Path(projected["alignment_path"]).name,
            record=projected["alignment_record_id"], original_locus=candidate["orig_Locus"].removeprefix("HML-2_"),
            locus=locus, sample=candidate["ID"], haplotype=candidate["Haplotype"], newly_admitted=1,
            candidate_source_count=1, ID_Full=candidate["ID_Full"], Source_Identifier=candidate["Source_Identifier"],
            source_binding_basis="current_source_exact_owning_caller_alignment",
            copy_suffixed_record=0, physical_locus_changed_from_extraction=0,
            aligned_canonical_bases=sum(base in "ACGT" for base in sequence),
            aligned_sequence_sha256=projected["typeII_KCON_sequence_sha256"], sequence=sequence))
        total_alignment_records += 1
    inputs.append({"path": str(args.additional_projections),
                   "sha256": sha256(args.additional_projections.read_bytes()).hexdigest()})
    seen_sources = {}
    # A source/copy contributes once. If its retained projection is repeated,
    # retain greatest canonical coverage, then deterministic filename/ID order.
    for entry in sorted(bound_records, key=lambda value: (-value["aligned_canonical_bases"], value["alignment"], value["record"])):
        source_key = entry["sample"], entry["haplotype"], entry["Source_Identifier"]
        assert entry["Source_Identifier"] not in {"", "NA"}, entry["ID_Full"]
        if source_key in seen_sources:
            assert seen_sources[source_key]["locus"] == entry["locus"], source_key
            exclusions.append({"alignment": entry["alignment"], "record": entry["record"], "reason": "duplicate_alignment_same_current_source_copy"})
            continue
        seen_sources[source_key] = entry
        sequence = entry.pop("sequence")
        admissions.append(entry)
        regions = {"LTR": owner.kcon_ltr_sequence(sequence)}
        for gene, (start, end) in ORF_COORDS.items():
            subsequence = sequence[start:end]
            regions[gene] = subsequence if sum(base in "ACGT" for base in subsequence) >= 0.5 * (end - start) else None
        for region, subsequence in regions.items():
            if subsequence is None:
                continue
            per_locus[region][entry["locus"]][subsequence] += 1
            observation_rows.append({"region": region, "locus": entry["locus"], "original_locus": entry["original_locus"],
                "alignment_record": entry["record"], "ID_Full": entry["ID_Full"], "Source_Identifier": entry["Source_Identifier"],
                "sample": entry["sample"], "haplotype": entry["haplotype"], "newly_admitted": entry["newly_admitted"],
                "sequence_sha256": sha256(subsequence.encode()).hexdigest()})
    admissions.sort(key=lambda value: (value["alignment"], value["record"]))
    assert len(admissions) + len(exclusions) == total_alignment_records
    assert len(admissions) == len(seen_sources)
    assert ACROCENTRIC <= set(per_locus["LTR"])
    assert sha256(catalog.read_bytes()).hexdigest() == catalog_digest, "Catalog changed during source binding"
    summary = {"catalog": str(catalog), "catalog_sha256": catalog_digest,
        "baseline_catalog": str(BASELINE), "baseline_catalog_sha256": sha256(BASELINE.read_bytes()).hexdigest(),
        "changed_catalog_columns": dict(changed_columns), "changed_catalog_records": len(changed_ids),
        "authenticated_recovery_ID_renames": len(renamed),
        "authenticated_HML11_family_exclusion_rows": len(family_changes),
        "HML11_comparator_roster": str(family_path),
        "HML11_comparator_roster_sha256": sha256(family_path.read_bytes()).hexdigest(),
        "source_recovery_mapping": str(recovery_path) if recoveries else None,
        "source_recovery_mapping_sha256": sha256(recovery_path.read_bytes()).hexdigest() if recoveries else None,
        "KCON_gene_coordinate_evidence": coordinate_evidence,
        "KCON_reference_sha256": sha256(kcon_path.read_bytes()).hexdigest(),
        "included_public_catalog_records": len(rows), "public_sample_haplotypes": 584,
        "retained_alignment_records_admitted": len(admissions),
        "alignment_records_inspected": total_alignment_records,
        "source_bound_copy_suffixed_records": sum(row["copy_suffixed_record"] for row in admissions),
        "cross_locus_reassigned_alignment_records": sum(row["physical_locus_changed_from_extraction"] for row in admissions),
        "source_bindings_by_evidence": dict(Counter(row["source_binding_basis"] for row in admissions)),
        "source_sequence_disambiguation_checks": len(source_checks),
        "excluded_non_HML2_loci": sorted(NON_HML2_LOCI),
        "all_admitted_acrocentric_sequences_have_full_query_zero_mismatch_current_source_alignment": True,
        "known_nonidentical_or_source_ambiguous_sequences_admitted": 0,
        "retained_sequence_limit": "All acrocentric, copy-suffixed, and multiple-source-slot admissions require a unique full-query zero-mismatch alignment to a current source FASTA. Query insertions and clipping are forbidden. Source insertions omitted by the KCON projection are explicitly counted in the CIGAR. Other unique original source slots retain their existing KCON sequences and were not globally re-generated or proven byte-identical to raw FASTAs. The documented mixed-orientation publication identity segments at 1q21.3 and 7q22.2 are not silently reoriented.",
        "one_observation_per_current_source_copy": True,
        "admitted_excluded_partition_conserved": True,
        "newly_admitted_records": dict(Counter(row["locus"] for row in admissions if row["newly_admitted"])),
        "newly_admitted_LTR_records": dict(Counter(row["locus"] for row in observation_rows if row["newly_admitted"] and row["region"] == "LTR")),
        "original_LTR_loci": len(old_ltr), "resolved_LTR_loci": len(per_locus["LTR"]),
        "new_admissions_have_unique_source_identity": True,
        "every_admitted_alignment_binds_one_current_source_identity": True,
        "admitted_multiple_candidate_slots_disambiguated_by_sequence": sum(row["candidate_source_count"] > 1 for row in admissions),
        "alignment_record_exclusion_reasons": dict(Counter(row["reason"] for row in exclusions)),
        "tree_count_definition": "n = retained source-bound aligned element observations passing the region's sequence-coverage threshold, counting each current source/copy once; up to two modal exact sequence clusters per physical locus are drawn. Includes exact-source-disambiguated copy-suffixed records, unlike the historical single-haplotype-name parser. Counts are not full catalog copy counts or protein-functional counts.",
        "retained_historical_ORF_tree_limitation": "The original pooled aligned FASTA and its tip-to-source map are absent locally. Historical locus__hapN__nN tips cannot authenticate repaired source assignment; the supplied current ORF trees are newly inferred from the retained KCON alignments.",
        "historical_ORF_locus_counts": {}, "regions": {},
        "methods": {"binding": "Bind alignment file's original extraction locus, sample, and normalized haplotype to included public PRESENT catalog source rows. Multiple-source slots, all #c records, and all acrocentric sources require a unique full-query zero-mismatch minimap2 asm20 alignment into one current source FASTA, with no query insertions or clipping. Omitted source insertions are counted explicitly. Other single unsuffixed source slots retain the existing source-bound KCON sequence. Bin by that source's repaired physical locus. Exclude source identity failures, unresolved physical loci, and sequence-authenticated HML-11 loci 8p22/17p13.1_hg38. Deduplicate current sample/haplotype/Source_Identifier, retaining maximal canonical coverage then filename/record lexical order.",
                    "LTR": "Existing owner's sequence/statistical method with corrected source-first physical-locus binding: 3-prime KCON bases 8504:9472, else 5-prime bases 0:968; at least 600 A/C/G/T bases; two most frequent exact aligned clusters per locus; pairwise canonical-site p-distance; Biopython neighbor joining; negative branch lengths clipped to zero.",
                    "ORFs": "Same catalog-to-KCON admission as LTR. Corrected zero-based half-open gag [1111,3112), pro [2913,3918), pol [3878,6749), env [6450,8550); each reference span ends at its stop codon. At least 50% canonical bases; up to two most frequent exact aligned gene clusters per locus; pairwise canonical-site nucleotide p-distance; Biopython neighbor joining; negative branch lengths clipped to zero. The historical nearest p-distance >0.30 flag (more than 30% differing canonical aligned sites) is diagnostic only, with no tip exclusions."}}
    if args.check_inputs:
        print(json.dumps({"status": "READY_TO_EXECUTE", **summary}, indent=2))
        return
    comparison = OUT / "historical_tree_comparison"
    comparison.mkdir(exist_ok=True)
    for path in ORIGINAL_PHY.glob("*.nwk"):
        if not (comparison / path.name).exists():
            shutil.copy2(path, comparison / path.name)
    for name in ("locus_subfamily.tsv", "expanded_ltr_subfamily_assignments.tsv"):
        shutil.copy2(ORIGINAL_PHY / name, OUT / name)
    owner.PHY = OUT
    region_trees = {}
    count_authority = {}
    selected_rows, distance_rows = [], []
    for region in per_locus:
        filename, counts, selected, distances = write_region_tree(region, per_locus[region])
        region_trees[region] = filename
        count_authority[region] = counts
        selected_rows.extend(selected)
        distance_rows.extend(distances)
        if region != "LTR":
            historical = owner.Phylo.read(comparison / filename, "newick")
            summary["historical_ORF_locus_counts"][region] = len({tip.name.split("__")[0] for tip in historical.get_terminals()})
    for region, filename in region_trees.items():
        tree = owner.Phylo.read(OUT / filename, "newick")
        tips = tree.get_terminals()
        assert {tip.name.split("__")[0] for tip in tips} == set(count_authority[region])
        assert all(clade.branch_length is None or clade.branch_length >= 0 for clade in tree.find_clades())
        summary["regions"][region] = {"tree": filename, "loci": len(count_authority[region]), "tree_tips": len(tips),
            "counted_aligned_observations": sum(count_authority[region].values()),
            "selected_cluster_observations": sum(owner.observation_count_from_tip(tip.name) for tip in tips),
            "acrocentric_counts": {locus: count_authority[region].get(locus, 0) for locus in sorted(ACROCENTRIC)},
            "newick_sha256": sha256((OUT / filename).read_bytes()).hexdigest()}
    from resolve_fragment_subfamilies import resolve as resolve_fragment_labels
    resolved_subfamilies = resolve_fragment_labels(catalog)
    with (OUT / "locus_subfamily.tsv").open("a") as handle:
        for row in resolved_subfamilies:
            handle.write(f"{row['locus']}\t{row['subfamily']}\n")
    known_subfamilies = owner.load_subfamily_authority()
    displayed_loci = set().union(*(set(counts) for counts in count_authority.values()))
    unassigned = sorted(displayed_loci - set(known_subfamilies))
    assert not unassigned, unassigned
    assert {known_subfamilies[locus] for locus in displayed_loci} <= {"LTR5_Hs", "LTR5A", "LTR5B", "non-LTR5"}
    summary["loci_without_LTR_subfamily_authority"] = []
    summary["subfamily_authority_supplement"] = resolved_subfamilies
    summary["displayed_subfamily_locus_counts"] = dict(Counter(known_subfamilies[locus] for locus in displayed_loci))
    fields = ["gene", "locus", "tip", "cluster_observations", "sequence_sha256", "sequence"]
    write_tsv(OUT / "selected_sequence_clusters.tsv", selected_rows, fields)
    write_tsv(OUT / "nearest_sequence_distance_diagnostics.tsv", distance_rows,
              fields + ["nearest_p_distance", "exceeds_historical_0_30_flag"])
    write_tsv(OUT / "outlier_exclusions.tsv", [], fields + ["nearest_p_distance"])
    write_tsv(OUT / "authenticated_source_ID_recoveries.tsv", renamed,
              ["original_zero_based_slot", "original_ID_Full", "current_ID_Full", "ID", "Haplotype", "orig_Locus"])
    write_tsv(OUT / "alignment_admission.tsv", admissions, list(admissions[0]))
    write_tsv(OUT / "alignment_exclusions.tsv", exclusions, list(exclusions[0]))
    write_tsv(OUT / "source_sequence_disambiguation.tsv", source_checks, list(source_checks[0]))
    write_tsv(OUT / "cross_locus_reassigned_alignment_records.tsv",
              [row for row in admissions if row["physical_locus_changed_from_extraction"]], list(admissions[0]))
    summary["source_admission_partition"] = write_source_admission_partition(rows, admissions, exclusions, source_checks)
    write_tsv(OUT / "region_observations.tsv", observation_rows, list(observation_rows[0]))
    write_tsv(OUT / "alignment_inputs.tsv", inputs, ["path", "sha256"])
    with (OUT / "selected_sequence_clusters.fasta").open("w") as handle:
        for row in selected_rows:
            handle.write(f">{row['gene']}|{row['tip']}\n{row['sequence']}\n")
    (OUT / "count_authority.json").write_text(json.dumps({region: dict(counts) for region, counts in count_authority.items()}, indent=2) + "\n")
    catalog_counts = []
    for gene in ORF_COORDS:
        before_counts, after_counts = [], []
        for selection, destination in ((baseline, before_counts), (rows, after_counts)):
            destination.extend(row["Locus"].removeprefix("HML-2_") for row in owner.proviral_rows(selection) if row[gene] not in MISSING_GENE)
        old, new = Counter(before_counts), Counter(after_counts)
        for locus in sorted(set(old) | set(new)):
            catalog_counts.append({"gene": gene, "locus": locus, "baseline_catalog_provirus_gene_calls": old[locus],
                                   "resolved_catalog_provirus_gene_calls": new[locus], "delta": new[locus] - old[locus]})
    write_tsv(OUT / "catalog_gene_call_count_changes.tsv", catalog_counts, list(catalog_counts[0]))
    # Reuse the current scientific renderer at its existing page-scale typography.
    renderer = inspect.getsource(owner.build_representative_phylogeny)
    renderer = renderer.replace('        "LTR5": PURPLE,', '        "LTR5": PURPLE,\n        "non-LTR5": "#343A40",')
    # Two alternating label columns leave room for 9-point type on a 9-inch page.
    # Every original terminal and branch remains in place; only leader lines move.
    renderer = renderer.replace("    x_label = x_tree * 1.07", "    x_label = x_tree * 1.03")
    renderer = renderer.replace("    for (tip_y, tip_x, label, color), target_y in zip(labels, label_y):", "    for label_index, ((tip_y, tip_x, label, color), target_y) in enumerate(zip(labels, label_y)):\n        target_x = x_label + (label_index % 2) * x_tree * 0.33")
    renderer = renderer.replace("            [tip_x, x_label * 0.985],\n            [tip_y, target_y],", "            [tip_x, x_label * 0.985, target_x * 0.985],\n            [tip_y, target_y, target_y],")
    renderer = renderer.replace("            x_label,\n            target_y,", "            target_x,\n            target_y,")
    renderer = renderer.replace("    ax.set_xlim(-0.065 * x_tree, x_tree * 1.45)", "    ax.set_xlim(-0.065 * x_tree, x_tree * 1.73)")
    renderer = renderer.replace("    path = PHY / output_name", "    fig.subplots_adjust(top=0.956, bottom=0.075)\n    if output_name.startswith('Main_'):\n        fig.text(0.012, 0.989, 'A' if 'LTR' in output_name else 'B', va='top', fontsize=12, weight='bold')\n    fig.text(0.51, 0.024, 'n = retained source-bound sequences; up to two modal clusters per locus', ha='center', va='bottom', fontsize=7)\n    fig.text(0.51, 0.006, 'Blue: LTR5Hs   Purple: LTR5A   Orange: LTR5B; HML-11 comparators excluded', ha='center', va='bottom', fontsize=7)\n    path = PHY / output_name")
    renderer = renderer.replace('    plt.close(fig)\n    return path', '    fig.savefig(path.with_suffix(".svg"), facecolor="white")\n    fig.savefig(path.with_suffix(".pdf"), facecolor="white")\n    fig.canvas.draw()\n    label_boxes = [(text.get_text(), text.get_window_extent(fig.canvas.get_renderer())) for text in ax.texts if re.search(r"\\(n=\\d+\\)", text.get_text())]\n    for index, (first_name, first_box) in enumerate(label_boxes):\n        for second_name, second_box in label_boxes[index + 1:]:\n            assert not first_box.overlaps(second_box), (output_name, first_name, second_name)\n    RENDER_QA.append({"panel": output_name, "locus_labels": len(label_boxes), "label_font_points": label_fontsize, "overlapping_locus_labels": 0})\n    plt.close(fig)\n    return path')
    owner.RENDER_QA = []
    exec(compile(renderer, "page_sized_phylogeny_renderer", "exec"), owner.__dict__)
    panels = [("gag", "Gag", "Supplement_Gag"), ("pro", "Pro", "Supplement_Pro"),
              ("env", "Env", "Supplement_Env"), ("LTR", "LTR", "Supplement_LTR"),
              ("pol", "Pol", "Supplement_Pol")]
    for region, title, basename in panels:
        owner.build_representative_phylogeny(region_trees[region], title, basename + ".png",
            mark_ltr_clades=True, count_authority=count_authority[region], figsize=(7.1, 9.0),
            label_fontsize=9.0, title_fontsize=11, axis_fontsize=9,
            clade_label_fontsize=8, star_size=100)
    summary["status"] = "REBUILT_RETAINED_SEQUENCE_PHYLOGENIES"
    summary["rendered_panels"] = [basename for _, _, basename in panels] + ["Main_LTR", "Main_Pol"]
    summary["render_geometry"] = owner.RENDER_QA
    summary["outlier_exclusions"] = []
    summary["nearest_sequence_distance_diagnostic_flags"] = sum(row["exceeds_historical_0_30_flag"] for row in distance_rows)
    assert sha256(catalog.read_bytes()).hexdigest() == catalog_digest, "Catalog changed during tree rebuild"
    (OUT / "verification.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    from build_focused_main_panels import main as build_focused_views
    build_focused_views()


if __name__ == "__main__":
    main()
