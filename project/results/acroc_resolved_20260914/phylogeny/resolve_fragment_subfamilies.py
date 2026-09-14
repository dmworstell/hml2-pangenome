#!/usr/bin/env python3
"""Authenticate omitted figure labels against the retained publication panel.

Sequence alignment authenticates the named published element; the subfamily
comes from the publication's LTR category, never an ORF-tree neighborhood.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
from hashlib import sha256
from pathlib import Path

sys.dont_write_bytecode = True
OUT = Path(__file__).resolve().parent
PROJECT = OUT.parents[2]
SOURCE = PROJECT / "working/type1_ltr5hs_matched_control_agent/sources"
PANEL = SOURCE / "PMC3228705_supplementary/1742-4690-8-90-S1.FAS"
FIGURE = SOURCE / "1742-4690-8-90-1_4x.jpg"
CATALOG = PROJECT / "results/resolved_manuscript_catalog_20260914/combined_hml2_orf_analysis.RESOLVED.tsv"
sys.path.insert(0, str(PROJECT / "manuscript"))
import build_narrative_main_figures as owner

# Figure 1 B-E was visually inspected. In particular, the current 4p16.3b
# element is the published 4p16.3a, and both HML-11 members are non-LTR5.
PUBLISHED = {
    "11q12.1": ("11q12.1", "LTR5B", "D"),
    "17p13.1_hg38": ("17p13.1", "non-LTR5", "E"),
    "1q21.3": ("1q21.3", "LTR5_Hs", "B"),
    "3q24_hg38": ("3q24", "LTR5_Hs", "B"),
    "4p16.3b": ("4p16.3a", "LTR5B", "D"),
    "6q25.1": ("6q25.1", "LTR5B", "D"),
    "7q22.2": ("7q22.2", "LTR5_Hs", "B"),
    "7q34": ("7q34", "LTR5_Hs", "B"),
    "8p22": ("8p22", "non-LTR5", "E"),
    "Xq11.1": ("Xq11.1", "LTR5B", "D"),
    "Xq12": ("Xq12", "LTR5A", "C"),
}


def write_tsv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def resolve(catalog):
    owner.CATALOG = catalog.resolve()
    rows, _ = owner.load_catalog()
    allowed = {(row["Locus"].removeprefix("HML-2_"), row["ID"], row["Haplotype"])
               for row in rows if row["observation_state"] == "PRESENT"}
    # The two positively identified HML-11 elements are comparator evidence,
    # never primary HML-2 admissions. Bind their sequence queries to the exact
    # authorized source roster, preserving the primary reader for all others.
    comparator_path = OUT.parent / "HML11_comparator_catalog_rows.tsv"
    with comparator_path.open(newline="") as handle:
        comparators = {row["ID_Full"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert len(comparators) == 1176
    assert {row["orig_Locus"] for row in comparators.values()} == {"HML-2_8p22", "HML-2_17p13.1_hg38"}
    comparator_queries = set()
    with owner.CATALOG.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["ID_Full"] not in comparators:
                continue
            comparator = comparators[row["ID_Full"]]
            assert all(row[field] == comparator[field] for field in ("ID_Full", "ID", "Haplotype", "Source_Identifier", "orig_Locus", "Locus"))
            assert row["analysis_include"] == "0" and row["analysis_exclusion_reason"] == "non_HML2_HML11_sequence_identity"
            if owner.PUBLIC_ID.fullmatch(row["ID"]) and row["observation_state"] == "PRESENT":
                comparator_queries.add((row["Locus"].removeprefix("HML-2_"), row["ID"], row["Haplotype"]))
    panel_records = list(owner.SeqIO.parse(PANEL, "fasta"))
    panel_path = OUT / "published_subfamily_identity_panel.fasta"
    with panel_path.open("w") as handle:
        for record in panel_records:
            sequence = str(record.seq).upper().replace("-", "")
            handle.write(f">{record.id}\n{sequence}\n")
    queries, evidence = {}, []
    query_path = OUT / "fragment_subfamily_identity_queries.fasta"
    with query_path.open("w") as handle:
        for locus, (published, family, panel) in PUBLISHED.items():
            path = PROJECT / "data" / f"{locus}_carrier_kcon.aln.fa"
            candidates = []
            maxima = [0, 0]
            for record in owner.SeqIO.parse(path, "fasta"):
                match = re.fullmatch(r"(.+)_(hap1|hap2|h1|h2|pat|mat)", record.id)
                if match is None:
                    continue
                sample, haplotype = match.groups()
                haplotype = {"hap1": "h1", "hap2": "h2"}.get(haplotype, haplotype)
                if (locus, sample, haplotype) not in (comparator_queries if family == "non-LTR5" else allowed):
                    continue
                aligned = str(record.seq).upper()
                sequence = aligned.replace("-", "")
                candidates.append((sum(base in "ACGT" for base in sequence), record.id, sequence))
                for index, (start, end) in enumerate(((0, 968), (8504, 9472))):
                    maxima[index] = max(maxima[index], sum(base in "ACGT" for base in aligned[start:end]))
            assert candidates, locus
            bases, identifier, query = max(candidates, key=lambda value: (value[0], value[1]))
            queries[locus] = query
            handle.write(f">{locus}\n{query}\n")
            evidence.append({"locus": locus, "published_locus": published, "subfamily": family,
                "evidence_class": "published_LTR_category_sequence_authenticated_identity",
                "publication_figure_panel": f"1{panel}", "admitted_alignment_records_checked": len(candidates),
                "maximum_5prime_LTR_canonical_bases": maxima[0], "maximum_3prime_LTR_canonical_bases": maxima[1],
                "identity_query_record": identifier, "identity_query_canonical_bases": bases,
                "identity_query_sha256": sha256(query.encode()).hexdigest(),
                "alignment_input": str(path), "alignment_input_sha256": sha256(path.read_bytes()).hexdigest(),
                "publication_source": "https://doi.org/10.1186/1742-4690-8-90",
                "classification_basis": "Published LTR category of sequence-authenticated element; not inferred from ORF clade"})
    command = ["minimap2", "-x", "asm20", "-c", "--secondary=yes", "-N", "100", "-p", "0.1",
               "-t", "2", str(panel_path), str(query_path)]
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    (OUT / "fragment_subfamily_identity.paf").write_text(result.stdout)
    (OUT / "fragment_subfamily_identity.stderr.txt").write_text(result.stderr)
    alignments = []
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        tags = {tag.split(":", 2)[0]: tag.split(":", 2)[2] for tag in fields[12:]}
        alignments.append({"query": fields[0], "query_length": int(fields[1]), "query_start": int(fields[2]),
            "query_end": int(fields[3]), "strand": fields[4], "target": fields[5],
            "target_length": int(fields[6]), "target_start": int(fields[7]), "target_end": int(fields[8]),
            "matching_bases": int(fields[9]), "alignment_block_length": int(fields[10]),
            "mapping_quality": int(fields[11]), "alignment_score": int(tags["AS"]),
            "CIGAR": tags["cg"], "edit_distance": int(tags["NM"])})
    for row in evidence:
        matches = sorted((match for match in alignments if match["query"] == row["locus"]),
                         key=lambda match: match["alignment_score"], reverse=True)
        assert matches, row["locus"]
        best = matches[0]
        runner = next((match for match in matches[1:] if match["target"] != best["target"]), None)
        assert best["target"] == "HML-2_" + row["published_locus"], (row["locus"], best)
        # Element identity does not require an unproven collinear projection.
        # Two retained queries have distinct near-exact opposite-strand pieces
        # of the same published element; keep that limitation in the receipt.
        same_target = [match for match in matches if match["target"] == best["target"]]
        covered = set()
        for match in same_target:
            assert match["matching_bases"] / match["alignment_block_length"] >= 0.90, match
            covered.update(range(match["query_start"], match["query_end"]))
        coverage = len(covered) / best["query_length"]
        identity = sum(match["matching_bases"] for match in same_target) / sum(match["alignment_block_length"] for match in same_target)
        assert coverage >= 0.90 and identity >= 0.90, (row["locus"], coverage, identity)
        assert runner is None or best["alignment_score"] > runner["alignment_score"], row["locus"]
        row.update({"verified_published_target": best["target"], "best_score": best["alignment_score"],
            "query_coverage": coverage, "identity": identity, "edit_distance": best["edit_distance"],
            "identity_alignment_segments": len(same_target),
            "query_strands_in_published_element": ";".join(sorted({match["strand"] for match in same_target})),
            "projection_collinearity": "mixed_strand_segments_not_claimed_collinear" if len({match["strand"] for match in same_target}) > 1 else "same_strand_alignment",
            "runner_up_target": runner["target"] if runner else "none_reported",
            "runner_up_score": runner["alignment_score"] if runner else "",
            "identity_crosswalk_status": "PASS"})
    write_tsv(OUT / "fragment_subfamily_identity_alignments.tsv", alignments)
    write_tsv(OUT / "fragment_subfamily_resolution.tsv", evidence)
    receipt = {"status": "ALL_11_PUBLICATION_IDENTITIES_AUTHENTICATED", "catalog": str(owner.CATALOG),
        "catalog_sha256": sha256(owner.CATALOG.read_bytes()).hexdigest(), "command": command,
        "official_panel": str(PANEL), "official_panel_sha256": sha256(PANEL.read_bytes()).hexdigest(),
        "visually_inspected_publication_figure": str(FIGURE), "publication_figure_sha256": sha256(FIGURE.read_bytes()).hexdigest(),
        "reference_panel_sequences": len(panel_records), "resolved_subfamily_labels": evidence,
        "HML11_comparator_roster": str(comparator_path), "HML11_comparator_roster_sha256": sha256(comparator_path.read_bytes()).hexdigest(),
        "HML11_role": "The exact roster's 8p22 and 17p13.1_hg38 sources are retained only as sequence-identity comparators; both loci are excluded from every primary HML-2 tree. The legacy admitted_alignment_records_checked field counts eligible identity-panel records, not primary HML-2 admission for these comparators.",
        "method": "Map the longest canonical retained eligible element per locus against every ungapped official Subramanian Additional File 1 element using minimap2 asm20, retaining secondary hits; require a unique best-scoring named target with union query coverage >=90% across >=90%-identity segments. The nine HML-2 loci use the actual primary reader, and the two HML-11 loci use only the exact authorized comparator source roster. Opposite-strand segments are explicitly reported and do not establish collinearity of the retained KCON projection. These gates authenticate element identity, not a new LTR-classification probability. Read the matched element's published LTR5_Hs/A/B or HML-11 category from visually inspected Figure 1 B-E."}
    (OUT / "fragment_subfamily_verification.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": receipt["status"], "results": evidence}, indent=2))
    return evidence


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", type=Path, default=CATALOG)
    resolve(parser.parse_args().catalog)
