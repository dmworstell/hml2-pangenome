"""Reproduce Table 2 from the frozen comparison inputs and locus type map.

Run with --data pointing to the extracted Supplementary_Data directory.
The added results preserve the original donor-by-locus unit and inclusion
rules. Intact means the original caller's Intact category. The combined
screen also admits its additional stop-free frameshift category. Each gene
is scored in its own annotated frame; canonical gag/pro/pol frame transitions
are not additional sequence frameshifts.
"""
import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path


def locus(value):
    return re.sub(r"HML-2_7p22\.1[ab]$", "HML-2_7p22.1", value)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    type_path = args.data / "Catalog/HML2_structural_and_ORF_catalog.tsv"
    types = defaultdict(set)
    with type_path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["analysis_include"] == "1" and row["provirus_type"] in ("type1", "type2"):
                types[locus(row["Locus"])].add(row["provirus_type"])
    assert all(len(values) == 1 for values in types.values()), "Mixed locus types require row-level typing"
    types = {key: next(iter(values)) for key, values in types.items()}
    catalogs, universes, sources = {}, {}, []
    for platform in ("long_read", "short_read"):
        path = args.data / f"Table_S6/inputs/{platform}.tsv.gz"
        sources.append({"path": str(path.relative_to(args.data)), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
        cells, people, loci = defaultdict(set), set(), set()
        with gzip.open(path, "rt") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["analysis_include"] != "1" or not re.fullmatch(r"(?:HG|NA)\d+", row["ID"]):
                    continue
                key = row["ID"], locus(row["Locus"])
                people.add(key[0])
                loci.add(key[1])
                flags = cells[key]
                structure = row["Structure"].lower()
                if "provirus" in structure or "multi" in structure:
                    flags.add("provirus")
                    if "multi" in structure:
                        flags.add("multi_copy")
                    for screen, allowed in (("intact", {"Intact"}), ("combined", {"Intact", "Intact_FS_End", "Frameshift_at_End"})):
                        for feature, genes in (("gag", ("gag",)), ("gag_pro", ("gag", "pro")), ("gag_pro_pol", ("gag", "pro", "pol")), ("env", ("env",))):
                            if all(row[gene] in allowed for gene in genes):
                                flags.add(screen + "_" + feature)
                                if feature == "env":
                                    flags.add(screen + "_" + types[key[1]] + "_env")
                elif "solo-ltr" in structure or "solo_ltr" in structure:
                    flags.add("solo_ltr")
        catalogs[platform], universes[platform] = cells, (people, loci)
    people = universes["long_read"][0] & universes["short_read"][0]
    loci = universes["long_read"][1] & universes["short_read"][1]
    assert (len(people), len(loci)) == (282, 83)
    results, per_locus, observations = [], [], []
    features = [(screen, feature) for feature in ("gag", "gag_pro", "gag_pro_pol", "env", "type2_env", "type1_env") for screen in ("intact", "combined")]
    features += [("structural", feature) for feature in ("provirus", "multi_copy", "solo_ltr")]
    for screen, feature in features:
        flag = feature if screen == "structural" else screen + "_" + feature
        keys = sorted(key for key, values in catalogs["long_read"].items() if key[0] in people and key[1] in loci and flag in values)
        by_locus = defaultdict(lambda: [0, 0])
        for key in keys:
            recovered = int(flag in catalogs["short_read"].get(key, set()))
            by_locus[key[1]][0] += recovered
            by_locus[key[1]][1] += 1
            observations.append({"donor": key[0], "locus": key[1], "screen": screen, "feature": feature, "recovered": recovered})
        n, N = sum(v[0] for v in by_locus.values()), len(keys)
        results.append({"feature": feature, "screen": screen, "n": n, "N": N, "percent": round(100 * n / N, 1) if N else None})
        for site, (n_site, N_site) in sorted(by_locus.items()):
            per_locus.append({"locus": site, "feature": feature, "screen": screen, "n": n_site, "N": N_site, "percent": round(100 * n_site / N_site, 1)})
    expected = {"gag": (1855,2534), "gag_pro": (745,1556), "gag_pro_pol": (266,569), "env": (1623,1944), "provirus": (10173,10587), "multi_copy": (0,225), "solo_ltr": (455,2168)}
    for row in results:
        if row["feature"] in expected and row["screen"] in ("combined", "structural"):
            assert (row["n"], row["N"]) == expected[row["feature"]], row
    for screen in ("intact", "combined"):
        lookup = {r["feature"]: r for r in results if r["screen"] == screen}
        for count in ("n", "N"):
            assert lookup["env"][count] == lookup["type1_env"][count] + lookup["type2_env"][count]
    for filename, rows in (("Table_2_ORF_recovery.tsv", results), ("Table_2_ORF_recovery_by_locus.tsv", per_locus), ("Table_2_ORF_recovery_observations.tsv", observations)):
        with (args.out / filename).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    sources.append({"path": str(type_path.relative_to(args.data)), "sha256": hashlib.sha256(type_path.read_bytes()).hexdigest()})
    (args.out / "provenance.json").write_text(json.dumps({"sources": sources, "donors": len(people), "loci": len(loci), "type_stratification": "Unique retained catalog type for each locus; all eligible loci are type-invariant in the catalog", "original_totals_reproduced": True, "results": results}, indent=2))
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
