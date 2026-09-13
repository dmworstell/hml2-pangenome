#!/usr/bin/env python3
"""Build the deterministic exact-four structural CNV/array manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import pathlib
import sys


HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "frozen" / "cnv_array_exact4_manifest.v1.json"
STAGING = ROOT / "inputs" / "cluster_execution_plans" / "cnv_array_staging_plan.json"
TRUTH = ROOT / "working" / "sevenp22_proxy_resolution_agent" / "haplotype_copy_number_truth.tsv"
REAUDIT = ROOT / "working" / "array_structural_reaudit_v1" / "haplotype_copy_truth.tsv"
REAUDIT_SUMMARY = ROOT / "working" / "array_structural_reaudit_v1" / "summary.json"
PHYSICAL = ROOT / "handoff" / "HML2_PARALOG_SYNTENY_PHYSICAL_COPY_ADJUDICATION_20260718.json"
PHYSICAL_EVIDENCE = ROOT / "handoff" / "HML2_PARALOG_SYNTENY_PHYSICAL_COPY_EVIDENCE_20260718.tsv"
CLUSTER_PYTHON = (
    "/path/to/hml2_workspace/condaenv/"
    "repeatmaskerenv/bin/python3.11"
)

SAMPLES = ("HG02027", "HG02178", "HG03669", "HG03834")
TRUTH_HAPS = {
    "HG02027": ("mat", "pat"),
    "HG02178": ("h1", "h2"),
    "HG03669": ("mat", "pat"),
    "HG03834": ("mat", "pat"),
}
ASSEMBLY_HAPS = {
    "HG02027": ("mat", "pat"),
    "HG02178": ("hap1", "hap2"),
    "HG03669": ("mat", "pat"),
    "HG03834": ("mat", "pat"),
}

# Read-only live authentication captured on Rocky9 on 2026-07-18.  The BAM
# records are assembly contigs aligned to T2T, not read alignments.
BAM_PINS = {
    ("HG02027", "mat"): (1088427831, "a69f9b9b184c601c51ea644a9d9fd799499d0c2bff64a6683bbe2ee80f4eba78", 1542496, "20758128b2ed15d7f4810ce12c111e7b2928c6e40a3a361ed9a3ba1779662c85", 1),
    ("HG02027", "pat"): (1098101353, "c5733c704e5556d0a808201ecabd74bd3cc264730796809b14743eb35f9cef4c", 1544328, "7310002b88ce434b8c7f7e2a81cfbbf3cb2a985a09f31b7debf0aded9d81a1a7", 1),
    ("HG02178", "hap1"): (1113230121, "4fc7bcf77c3c052ac2b49b6673dc42c0c31e1779375239c764a57859c335f2d8", 1510656, "cc958757dc79211e1a7cebb16f8697a9edc29ed5de502f7cca65ed55383b1ecc", 2),
    ("HG02178", "hap2"): (1147808962, "ba46f65663274f516120718f851dc6db9e555f1fde9ddfb60647737426453df1", 1542672, "fc9262c2443bedd5f60ccc4bd2b5dba7d57a3a6c41db57bb2a3d9b5b1fb82b50", 1),
    ("HG03669", "mat"): (1129166095, "953479979b0ab6542619ef837caa24e446bea485157d4f054116400d76969afb", 1540072, "a9905a50ac2f63d0e47cd4e60ee21ad2217c6c344a42831989958d056d0419c1", 1),
    ("HG03669", "pat"): (1112054154, "6dfe9eb6d3beacb4f7d81a0e73855fdaf63b1401c0481f588f81b7b87de6659f", 1522120, "ae324cf036ff7eb5e9a28df262bbce68e46a80fad641b60990db6cb97c24024e", 1),
    ("HG03834", "mat"): (1142737600, "0d7a4cec2e9577632ad64b9e333f3675ec5e4d859d3304a5de62a8d80d1d9328", 1539560, "ddd7c63f3ce8de6a7a2ba4e5190804aa2f159889941cc37ea2ff59ed739f6d28", 1),
    ("HG03834", "pat"): (1127544601, "b35625a38cedf49d2a8d6029984d4a51b4260352ac5eb3129ba7eb54994d39f2", 1540632, "e5dd164833cd7c31bb572ec2166503ea39c3af94a1a9ab1596ed2b609a4e418b", 1),
}


def canonical_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n").encode("utf-8")


def sha(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def artifact(path: pathlib.Path, role: str) -> dict[str, object]:
    return {"role": role, "repository_path": str(path.relative_to(ROOT)), "size_bytes": path.stat().st_size, "sha256": sha(path)}


def load_truth() -> dict[tuple[str, str], dict[str, object]]:
    with TRUTH.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    out: dict[tuple[str, str], dict[str, object]] = {}
    for row in rows:
        key = (row["sample"], row["haplotype"])
        if row["sample"] in SAMPLES:
            out[key] = {
                "truth_haplotype": row["haplotype"],
                "array_state": row["array_state"],
                "physical_copy_number": int(row["array_copy_number"]),
                "multi_vs_single_binary": int(row["multi_vs_single_binary"]),
                "high_copy_ge3_binary": int(row["high_copy_ge3_binary"]),
            }
    expected = {(s, h) for s in SAMPLES for h in TRUTH_HAPS[s]}
    if set(out) != expected:
        raise SystemExit("exact-four truth rows are missing or duplicated")
    return out


def build() -> dict[str, object]:
    staging = json.loads(STAGING.read_text(encoding="utf-8"))
    assemblies = {(x["sample_id"], x["haplotype"]): x for x in staging["exact_four_inventory"]["assembly_authorities"]}
    truth = load_truth()
    bam_root = "/path/to/hml2_workspace/HML2_project/global_alignments_bam"
    run_root = "/path/to/hml2_workspace/HML2_project/cnv_array_v3/runs/run_20260718_003"
    bundle_root = f"{run_root}/bundle"
    units: list[dict[str, object]] = []
    for index, sample in enumerate(SAMPLES):
        haps: list[dict[str, object]] = []
        for assembly_hap, truth_hap in zip(ASSEMBLY_HAPS[sample], TRUTH_HAPS[sample]):
            assembly = assemblies[(sample, assembly_hap)]
            stem = pathlib.Path(str(assembly["fasta_path"])).name[:-3]
            bam_size, bam_sha, bai_size, bai_sha, records = BAM_PINS[(sample, assembly_hap)]
            record = dict(truth[(sample, truth_hap)])
            record.update({
                "assembly_haplotype": assembly_hap,
                "haplotype_alias_is_explicit": assembly_hap != truth_hap,
                "assembly_fasta": {"path": assembly["fasta_path"], "sha256": assembly["fasta_sha256"]},
                "assembly_fai": {"path": assembly["fai_path"], "sha256": assembly["fai_sha256"]},
                "assembly_alignment_bam": {
                    "path": f"{bam_root}/{stem}.global_t2t.sorted.bam",
                    "size_bytes": bam_size,
                    "sha256": bam_sha,
                    "semantic_type": "assembly_contig_to_t2t_alignment_not_read_depth",
                },
                "assembly_alignment_bai": {
                    "path": f"{bam_root}/{stem}.global_t2t.sorted.bam.bai",
                    "size_bytes": bai_size,
                    "sha256": bai_sha,
                },
                "expected_indexed_window_record_count": records,
            })
            haps.append(record)
        units.append({
            "index": index,
            "sample_id": sample,
            "selection_reason": "both_legacy_coverage_pdfs_absent_and_no_accepted_cancelled_run_completion",
            "haplotypes": haps,
            "outputs": {
                "evidence_tsv": f"{run_root}/units/{sample}/assembly_structural_evidence.tsv",
                "checkpoint_json": f"{run_root}/units/{sample}/CHECKPOINT.json",
            },
        })

    source_files = [
        HERE / "src" / "cnv_array_plan.py",
        HERE / "src" / "01_cnv_array_unit.sh",
        HERE / "src" / "02_cnv_array_closure.sh",
        HERE / "src" / "03_cnv_array_receipt_controller.sh",
    ]
    manifest: dict[str, object] = {
        "schema_version": "hml2.cnv-array-structural-closure-manifest.v1",
        "run_id": "run_20260718_003",
        "scope": {
            "locus": "HML-2_7p22.1",
            "reference_build": "T2T-CHM13v2.0",
            "core_region_1based": "chr7:4699540-4717514",
            "validation_region_1based": "chr7:4649540-4767514",
            "ordered_samples": list(SAMPLES),
            "selected_sample_count": 4,
            "selected_haplotype_count": 8,
        },
        "selection": {
            "cancelled_queue_is_empty": True,
            "cancelled_run_output_status": "missing_or_invalid_no_accepted_completion_receipt",
            "legacy_pdf_status_by_sample": {s: {"coverage": "absent", "dedup_coverage": "absent"} for s in SAMPLES},
            "accepted_historical_samples_excluded": ["HG02040", "HG02293", "HG02841", "HG03471", "HG04115", "NA20282", "NA21110"],
            "missing_source_samples_excluded": ["HG00512", "HG00514", "HG01457", "HG02011", "HG02018", "HG02059"],
        },
        "authority_bindings": {
            "exact_copy_truth": artifact(TRUTH, "accepted_7p22_exact_haplotype_copy_truth"),
            # Blind-spot review (Lane G, 2026-07-19; Codex fourth-rulings s3.3):
            # for HML-2_7p22.1 run_audit.py copies the copy number straight from
            # the exact-truth TSV (``cn = int(truth["array_copy_number"])`` at
            # run_audit.py:118-121), so the reaudit is a same-substrate
            # corroboration of that truth, NOT an independent measurement, and
            # must NOT be counted as a second authority for 7p22.1.  Only
            # HML-2_1p31.1b is derived independently (from the ORF catalog /
            # part structure).  Relabelled from the prior "independent_" role.
            "structural_reaudit_truth": artifact(
                REAUDIT, "same_substrate_structural_reaudit_corroboration_7p22_passthrough_1p31_independent"),
            "structural_reaudit_summary": artifact(
                REAUDIT_SUMMARY, "same_substrate_structural_reaudit_summary_not_independent_for_7p22"),
            "structural_reaudit_independence_note": {
                "HML-2_7p22.1": "same_substrate_passthrough_of_exact_copy_truth_not_independent_not_a_second_authority",
                "HML-2_1p31.1b": "independent_from_orf_catalog_part_structure",
                "provenance": "run_audit.py:118-121 copies array_copy_number from the exact-truth TSV for 7p22.1",
            },
            "physical_copy_adjudication": artifact(PHYSICAL, "paralog_synteny_physical_copy_adjudication"),
            "physical_copy_evidence": artifact(PHYSICAL_EVIDENCE, "paralog_synteny_physical_copy_evidence"),
            "cluster_assembly_authority": staging["exact_four_inventory"]["resident_authorities"][0],
            "policy": {
                "unit": "physical_copy",
                "seven_p22_tandem_units_count_once": True,
                "merge_with_paralog_or_segdup_families": False,
                "alignment_record_count_is_copy_number": False,
                "read_depth_inference_authorized": False,
            },
        },
        "runtime": {
            "host": "login.cluster.example",
            "observed_login_fqdn": "login-p03.pax.tufts.edu",
            "os": "Rocky Linux 9.6",
            "partition": "batch",
            "module_policy": "lowercase_non_deprecated_only",
            "modules": ["samtools/1.21"],
            "required_environment": {
                "APPTAINER_BIND": "/cluster",
                "HML2_CLUSTER_PYTHON": CLUSTER_PYTHON,
                "LC_ALL": "C",
                "PYTHONDONTWRITEBYTECODE": "1",
            },
            "tools": {
                "bash": {"path": "/usr/bin/bash", "sha256": "ec6d007d48ef11bc47ad3f372b4b20ff2f0d4e63867e7e4cc0f1b17b19fa88b2"},
                "python": {
                    "path": CLUSTER_PYTHON,
                    "realpath": CLUSTER_PYTHON,
                    "size_bytes": 32882744,
                    "sha256": "d96a013a2278ecb1172d21ca695024ff1e4467828da246c45b3eba9772f7e72a",
                    "version": "Python 3.11.0",
                },
                "samtools": {"path": "/path/to/software/samtools/1.21/bin/samtools", "sha256": "bb3a30b5203ba308f9309187b1d5558caf614196f885114736ae9d20f55cf23f", "version": "1.21"},
                "sbatch": {"path": "/usr/bin/sbatch", "sha256": "a9905ca17089454bed9aad58ddd040110fcb83c3483c688cba4b089375d33499", "version": "23.11.11"},
            },
        },
        "execution_paths": {
            "bundle_root": bundle_root,
            "manifest": f"{bundle_root}/frozen/cnv_array_exact4_manifest.v1.json",
            "planner": f"{bundle_root}/src/cnv_array_plan.py",
            "unit_script": f"{bundle_root}/src/01_cnv_array_unit.sh",
            "closure_script": f"{bundle_root}/src/02_cnv_array_closure.sh",
            "receipt_controller": f"{bundle_root}/src/03_cnv_array_receipt_controller.sh",
        },
        "execution_policy": {
            "network_allowed": False,
            "aws_allowed": False,
            "download_allowed": False,
            "minimap2_allowed": False,
            "new_reference_or_index_construction_allowed": False,
            "combined_reference_construction_allowed": False,
            "overwrite_allowed": False,
            "submission_authorized": False,
            "submission_tool_in_package": False,
        },
        "storage": {
            "run_root": run_root,
            "protected_roots": [
                "/path/to/hml2_workspace/HML2_project/scripts/hprc/check_CNVs/cnv_bams",
                "/path/to/hml2_workspace/HML2_project/scripts/hprc/check_CNVs/nucfreq_results",
                "/path/to/hml2_workspace/HML2_project/scripts/hprc/check_CNVs/cnv_summary_output",
                "/path/to/hml2_workspace/HML2_project/scripts/hprc/check_CNVs/slurm_logs",
            ],
            "preflight_receipt": f"{run_root}/preflight/PREFLIGHT.json",
            "closure_tsv": f"{run_root}/closure/array_truth_exact4.tsv",
            "closure_json": f"{run_root}/closure/RUN_CLOSURE.json",
        },
        "scheduling": {
            "unit_array": {"array": "0-3%4", "task_count": 4, "cpus_per_task": 1, "memory_gib": 4, "wall_minutes": 20, "dependency": None},
            "closure": {"task_count": 1, "cpus_per_task": 1, "memory_gib": 1, "wall_minutes": 5, "dependency": "receipt_census:unit_checkpoints"},
            "total_scheduled_tasks": 5,
            "afterany_allowed": False,
        },
        "bundle_files": [{"path": str(p.relative_to(HERE)), "size_bytes": p.stat().st_size, "sha256": sha(p)} for p in source_files],
        "units": units,
        "plan_sha256": "",
    }
    material = dict(manifest)
    material.pop("plan_sha256")
    manifest["plan_sha256"] = hashlib.sha256(canonical_bytes(material)).hexdigest()
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    payload = canonical_bytes(build())
    if args.check:
        if not OUT.exists() or OUT.read_bytes() != payload:
            print(f"stale or missing manifest: {OUT}", file=sys.stderr)
            return 1
        print(f"OK {OUT.relative_to(ROOT)} {hashlib.sha256(payload).hexdigest()}")
        return 0
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_bytes(payload)
    print(f"wrote {OUT.relative_to(ROOT)} {hashlib.sha256(payload).hexdigest()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
