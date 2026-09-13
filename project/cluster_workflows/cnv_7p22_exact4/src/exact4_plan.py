#!/usr/bin/env python3
"""Fail-closed exact-four planning and release-gated local publication.

This module is intentionally offline.  It contains no provider, network,
scheduler, SSH, deployment, submission, retry, or cancellation client.  The
only filesystem mutations it exposes are run-bound create-only claims and
same-filesystem no-replace publication after validation of a distinct
execution-release document.
"""

from __future__ import annotations

import argparse
from collections import Counter
import copy
import csv
import ctypes
import hashlib
import importlib.util
import json
import math
import os
import pathlib
import re
import stat
import subprocess
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

sys.dont_write_bytecode = True


class ExactFourPlanError(RuntimeError):
    """An exact-four schema, provenance, identity, or release gate failed."""


SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
SAFE_TOKEN_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.:-]{0,127}$")
SOURCE_RE = re.compile(r"^([^:\s]+):([0-9]+)-([0-9]+)$")

DEPLOYMENT_SCHEMA = "hml2_7p22_exact4_deployment_manifest_1"
RESIDENT_SCHEMA = "hml2_7p22_exact4_resident_receipt_1"
PROVIDER_SCHEMA = "hml2_7p22_exact4_provider_receipt_2"
TOOL_FIXTURE_SCHEMA = "hml2_7p22_exact4_resident_tool_fixture_receipt_1"
TOOL_FIXTURE_AUTH_SCHEMA = "hml2_7p22_exact4_resident_tool_fixture_authorization_1"
TOOL_FIXTURE_AUTH_ACCEPTANCE_SCHEMA = "hml2_7p22_exact4_resident_tool_fixture_authorization_acceptance_1"
PLAN_SCHEMA = "hml2_7p22_exact4_target_only_plan_1"
RELEASE_SCHEMA = "hml2_7p22_exact4_execution_release_1"
CLAIM_SCHEMA = "hml2_7p22_exact4_claim_1"
PREPARED_SCHEMA = "hml2_7p22_exact4_prepared_1"
PART_SCHEMA = "hml2_7p22_exact4_part_1"
ANALYSIS_SCHEMA = "hml2_7p22_exact4_analysis_1"
RUN_CLOSURE_SCHEMA = "hml2_7p22_exact4_run_closure_1"
SUMMARY_SCHEMA = "hml2_7p22_exact4_coverage_ratio_summary_1"
OUTPUT_SCHEMAS = {
    "prepared": PREPARED_SCHEMA,
    "part": PART_SCHEMA,
    "analysis": ANALYSIS_SCHEMA,
    "run_closure": RUN_CLOSURE_SCHEMA,
    "depth_summary": SUMMARY_SCHEMA,
}

COVERAGE_RATIO_SCHEMA_IDENTITY: dict[str, object] = {
    "callability_pairs": [["callable", "callable"], ["uncallable", "empty_own_flank_set"], ["uncallable", "own_flank_mean_below_5x"]],
    "element_record_fields": ["schema_id", "schema_version", "run_id", "sample_id", "locus_id", "haplotype_id", "contig", "source_start_0based", "source_end_0based", "product", "element_record_key", "flank_record_key", "source_depth_sum", "source_base_count", "S_exact", "R_exact", "callability", "callability_reason", "read_inferred_copy_number"],
    "element_record_key_fields": ["sample_id", "locus_id", "haplotype_id", "contig", "source_start_0based", "source_end_0based", "product"],
    "element_record_schema_id": "hml2_7p22_exact4_element_ratio_record_1",
    "element_record_schema_version": 1,
    "element_tsv_columns": ["schema_id", "schema_version", "run_id", "sample_id", "locus_id", "haplotype_id", "contig", "source_start_0based", "source_end_0based", "product", "element_record_key", "flank_record_key", "source_depth_sum", "source_base_count", "S_numerator", "S_denominator", "R_numerator", "R_denominator", "callability", "callability_reason", "read_inferred_copy_number"],
    "element_tsv_schema_id": "hml2_7p22_exact4_element_ratio_tsv_1",
    "flank_record_fields": ["schema_id", "schema_version", "run_id", "sample_id", "locus_id", "haplotype_id", "product", "flank_record_key", "flank_depth_sum", "flank_base_count", "F_exact", "callability", "callability_reason", "read_inferred_copy_number"],
    "flank_record_key_fields": ["sample_id", "locus_id", "haplotype_id", "product"],
    "flank_record_schema_id": "hml2_7p22_exact4_haplotype_flank_record_1",
    "flank_record_schema_version": 1,
    "flank_tsv_columns": ["schema_id", "schema_version", "run_id", "sample_id", "locus_id", "haplotype_id", "product", "flank_record_key", "flank_depth_sum", "flank_base_count", "F_numerator", "F_denominator", "callability", "callability_reason", "read_inferred_copy_number"],
    "flank_tsv_schema_id": "hml2_7p22_exact4_haplotype_flank_tsv_1",
    "json_missing_rational": None,
    "product_order": ["raw", "mapq10"],
    "rational_fields": ["schema_id", "schema_version", "numerator", "denominator"],
    "rational_schema_id": "hml2_7p22_exact4_reduced_rational_1",
    "rational_schema_version": 1,
    "record_order": ["accepted_target_table_sample_order", "declared_assembly_authority_haplotype_order", "utf8_bytewise_contig_order", "numeric_source_start_order", "numeric_source_end_order", "raw_then_mapq10"],
    "summary_fields": ["schema_id", "schema_version", "run_id", "sample_id", "locus_id", "haplotype_flank_records", "element_ratio_records", "read_inferred_copy_number"],
    "summary_schema_id": "hml2_7p22_exact4_coverage_ratio_summary_1",
    "summary_schema_version": 1,
    "tsv_missing_rational": "NA",
}
COVERAGE_RATIO_SCHEMA_SHA256 = "0938c88fafb56729f6062066357bc1e2af12d687fd2ad5f9a4ba4c093bbbe4b5"
EXECUTION_CONTRACT_SHA256 = "2c94ca36f31eb140f3feabfe6db535deb433372314dff153ddc9571b5e65c285"
EXECUTION_CONTRACT_ACCEPTANCE_SHA256 = "227fd968a03c4a2fefa180c8d0de01ea1f593d6620c06addca22055c190b0ef3"
COVERAGE_RATIO_DECISION_SHA256 = "ae81fd2e8825039f995cdc083cd76b4e5bf306b35655395b565d44cc86cc2851"
COVERAGE_RATIO_DECISION_ACCEPTANCE_SHA256 = "647b46014ebc9af5df38a3f2518ce815f4dfb34508e83ddcfb9f680a908f0902"
IMPLEMENTATION_CONTRACT_SHA256 = "936fb03189d8abcbaccb610eb107aa6dd60c7b95df5149c8362c78714a603040"
IMPLEMENTATION_CONTRACT_ACCEPTANCE_SHA256 = "3cabbb52b19d6aba6e178ae577a6b5e27158b49b74c0be10c1c26db9b5ce7773"
CONTRACT_BINDINGS = {
    "execution_contract_sha256": EXECUTION_CONTRACT_SHA256,
    "execution_contract_acceptance_sha256": EXECUTION_CONTRACT_ACCEPTANCE_SHA256,
    "coverage_ratio_decision_sha256": COVERAGE_RATIO_DECISION_SHA256,
    "coverage_ratio_decision_acceptance_sha256": COVERAGE_RATIO_DECISION_ACCEPTANCE_SHA256,
    "implementation_contract_sha256": IMPLEMENTATION_CONTRACT_SHA256,
    "implementation_contract_acceptance_sha256": IMPLEMENTATION_CONTRACT_ACCEPTANCE_SHA256,
}

LOCUS = "HML-2_7p22.1"
TARGETS = ("HG02027", "HG02178", "HG03669", "HG03834")
PART_COUNTS = (3, 3, 4, 3)
RUN_ROOT = "/path/to/hml2_workspace/HML2_project/cnv_7p22_exact4/runs"
TARGET_TABLE_BYTES = (
    b"sample_id\tlocus\n"
    b"HG02027\tHML-2_7p22.1\n"
    b"HG02178\tHML-2_7p22.1\n"
    b"HG03669\tHML-2_7p22.1\n"
    b"HG03834\tHML-2_7p22.1\n"
)
TARGET_TABLE_SHA256 = "c90e51fa2b8186dd9fa2ec501afa476c4d3dd2526c194d2a93acfb04ea3b3b0d"

PREREQUISITE_PINS = (
    ("handoff/7P22_EXACT_FOUR_TARGET_CLUSTER_EXECUTION_CONTRACT.md", EXECUTION_CONTRACT_SHA256),
    ("handoff/7P22_EXACT_FOUR_TARGET_CLUSTER_EXECUTION_CONTRACT_ACCEPTANCE.md", "6536cf0218f40634a7e6b20d0551414029476b43697961364a8fae4e712990c1"),
    ("handoff/7P22_EXACT_FOUR_TARGET_CLUSTER_EXECUTION_CONTRACT_AMENDED_ACCEPTANCE.md", EXECUTION_CONTRACT_ACCEPTANCE_SHA256),
    ("handoff/7P22_EXACT_FOUR_COVERAGE_RATIO_DECISION.md", COVERAGE_RATIO_DECISION_SHA256),
    ("handoff/7P22_EXACT_FOUR_COVERAGE_RATIO_DECISION_ACCEPTANCE.md", COVERAGE_RATIO_DECISION_ACCEPTANCE_SHA256),
    ("handoff/7P22_EXACT_FOUR_LIVE_INPUT_PREFLIGHT_RECEIPT.md", "662138df41bb5212a3008641e38d19e028cb48de83ea89775b6703ce590ad35a"),
    ("handoff/7P22_FOUR_TARGET_DRY_RUN_CONTRACT.md", "42abd8eb287ac8a7c40a783a600aabad77628c3248e2cbe094be679e82254d4b"),
    ("handoff/7P22_FOUR_TARGET_DRY_RUN_CONTRACT_ACCEPTANCE.md", "70ab1a506ea65f866c041f46d197b5ec1a77f39b0eb78e8024354ad9aaf1d5b8"),
    ("handoff/7P22_FOUR_TARGET_DRY_RUN_IMPLEMENTATION_ACCEPTANCE.md", "fe16e7d860b30859a4901053be93c9abe7bbdc1ce98632ba9a0feee53db296c6"),
    ("handoff/GENERIC_CNV_V2_WORKFLOW_CORRECTION_IMPLEMENTATION_ACCEPTANCE.md", "5460652d0173acc7b4e4445aff1326ba407b267d99c3c3bd55f2bfaaaecdf761"),
    ("handoff/7P22_EXACT_FOUR_EXECUTION_IMPLEMENTATION_CONTRACT.md", IMPLEMENTATION_CONTRACT_SHA256),
    ("handoff/7P22_EXACT_FOUR_EXECUTION_IMPLEMENTATION_CONTRACT_ACCEPTANCE.md", IMPLEMENTATION_CONTRACT_ACCEPTANCE_SHA256),
    ("generative/cluster_cnv_retry/render_7p22_exact_four_target_dry_run.py", "34c93e4c0ab36e42e886b4a4b3b724d67a04cdcb131224f43017f7ef2700dcbd"),
    ("generative/cluster_cnv_retry/test_render_7p22_exact_four_target_dry_run.py", "c8338301ce32db8d636361484a0dd6af9c975dfb9970dfb37451fa02053e413d"),
    ("handoff/7P22_EXACT_FOUR_EXPECTED_PROVIDER_METADATA_AUTHORITY_ACCEPTANCE.md", "cd18055697f2433933d2ea154718237d41e70265300c024d5d27d35a24e4117b"),
    ("handoff/7P22_EXACT_FOUR_AWS_TOOL_AUTHORITY_ACCEPTANCE.md", "59710dc0a77ca92f52b4a531a6e974b3706a1d45649da48fba193896ab994177"),
    ("handoff/7P22_EXACT_FOUR_LIVE_RECEIPT_CAPTURE_AUTHORIZATION.md", "b7d2a267e44291590868e2f679beb4c7f62790f1b89778e804f105266ab04d58"),
    ("handoff/7P22_EXACT_FOUR_LIVE_RECEIPT_CAPTURE_AUTHORIZATION_ACCEPTANCE.md", "78210dfb79668a0bb01a059d493b7796765b821b03955ea448227eab43ed3902"),
    ("handoff/7P22_EXACT_FOUR_CURRENT_RESIDENT_REHASH_RECEIPT.md", "228677eca4e5730bfc89f154404138481d62a94160d0e824e17354d1cbc68b75"),
    ("handoff/7P22_EXACT_FOUR_CURRENT_RESIDENT_REHASH_RECEIPT_ACCEPTANCE.md", "47d45db21fd22f7c4c6ecdcb310079113a6cff7c0740107c5c2e95f0ff7e8b93"),
    ("handoff/7P22_EXACT_FOUR_RECEIPT_SCHEMA_AND_RESIDENT_TOOL_REPIN_AMENDMENT.md", "c15c10fba7bf9cfc8aa6d247e719b85c19d9b44a760ee0f30d20682bbc8deeb3"),
    ("handoff/7P22_EXACT_FOUR_RECEIPT_SCHEMA_AND_RESIDENT_TOOL_REPIN_AMENDMENT_ACCEPTANCE.md", "5df78fbefa5977780cd25cc14fca2f2077da85c39d5c3f63fe971aeca3ba89b1"),
    ("handoff/7P22_EXACT_FOUR_V2_IMPLEMENTATION_REVIEW_CORRECTION_ADDENDUM.md", "58c55f29ee4eb6bda035f15a60e0c2e12529ca05d1483c12ba27046e44e7a89a"),
    ("handoff/7P22_EXACT_FOUR_V2_IMPLEMENTATION_REVIEW_CORRECTION_ADDENDUM_ACCEPTANCE.md", "1efa58eb0fe121ceeee0d298325cd19f5c2e0d9aac231f442b81764230fe3482"),
)
CLASSIFICATION_PROOF_SHA256 = "e35bbd6595b05f6a54f6fb1b9613dc7e2f4489d5e755b174529996960631bb7c"

RESIDENT_AUTHORITIES = (
    ("assembly_authority", "/path/to/hml2_workspace/HML2_project/cnv_retry_preflight_v1/run_20260713_001/artifacts/assembly_authority_11.v1.tsv", 7266, "467438608bbffa5a416d4d467fcb4a7d62ce26628cda0365ed716d71ea0bc31c"),
    ("combined_table", "/path/to/hml2_workspace/HML2_project/orf_analysis_combined/combined_hml2_orf_analysis.tsv", 85204182, "dcd78954ae6d4b09e41433367b6892c09d704a09b7031025401b9798ec5ebc29"),
    ("type2_kcon", "/path/to/hml2_workspace/HML2_1q22_project/ref/type2_KCON.fa", 9485, "a1734229fd22e804267f3763c1d38b837bd12ce329d7b823252d9a3c6abdb7ec"),
    ("ont_index", "/path/to/hml2_workspace/HML2_project/scripts/hprc/check_CNVs/data_ont_pre_release.index.csv", 615215, "de44adeb50f2b58af4c607ebf9fab4150979b204af6427cf89633a5ea76a9868"),
)

ASSEMBLY_AUTHORITIES = (
    ("HG02027", "mat", "/path/to/hml2_workspace/MOVED/unzipped_assemblies/HG02027_mat_hprc_r2_v1.0.1.fa", "26458dd7163e6ccfbe5f138a30261887492e9f0353809f72f9cf161df40255dc", "353dba6c0f9e493d0648b6ef81858c46633f30a288ad43166d295c87e467b35b"),
    ("HG02027", "pat", "/path/to/hml2_workspace/MOVED/unzipped_assemblies/HG02027_pat_hprc_r2_v1.0.1.fa", "667b2a52afdc72326f3d8bd3d55dc9ad38b7a3fe46d43e002ad7c4c7b493ea15", "52e7524a59a11f5a753aecdea534a48f333d6b4e3f88f2744a37e964f57fd826"),
    ("HG02178", "hap1", "/path/to/hml2_workspace/MOVED/unzipped_assemblies/HG02178_hap1_hprc_r2_v1.0.1.fa", "28a3cbc2a212e9123728c9003cd6d8a01755b87e1860382992013abdcfaff127", "82e50e5e6f0ece3f57ecee88354a8ffe9fe12c627ccd79590fe107680982ab04"),
    ("HG02178", "hap2", "/path/to/hml2_workspace/MOVED/unzipped_assemblies/HG02178_hap2_hprc_r2_v1.0.1.fa", "65914f792c82bb8dff43d80bd950419793142369ae48cf08b8ed90f0a13daaff", "28745f5134df2014f897e528dbb41f1ce8bf682883d8361fb00d8be7f158574d"),
    ("HG03669", "mat", "/path/to/hml2_workspace/MOVED/unzipped_assemblies/HG03669_mat_hprc_r2_v1.0.1.fa", "3a81300dd49b14ca8867bd17f80285e831f037f557498cc3f265634d2baa9df3", "ea52d66fc3df56482835d8cf6550c429b87b4c3c023d092c07a4024f3690d644"),
    ("HG03669", "pat", "/path/to/hml2_workspace/MOVED/unzipped_assemblies/HG03669_pat_hprc_r2_v1.0.1.fa", "63a573ec91cd7f9e584a22d6aed9a345d4daa204108ebd5c54e981a9b0d6cf8b", "311ceae9acfebef2d148433b98dd1bdf64ec6f99fb4caf4d2b73618344bbf283"),
    ("HG03834", "mat", "/path/to/hml2_workspace/MOVED/unzipped_assemblies/HG03834_mat_hprc_r2_v1.0.1.fa", "6b1ca9533db953db99e4d7d485883ca97ab0d64e581cc30abad54737e40d57ad", "42300d9c7baf17b1fc6c6774c319e8df8077e1e2b784579ee3e7a24b6b9d0f25"),
    ("HG03834", "pat", "/path/to/hml2_workspace/MOVED/unzipped_assemblies/HG03834_pat_hprc_r2_v1.0.1.fa", "8b9e5386d23061c3e28c27b90acd8f8959333927ec8edd217e3acf4c313ca56c", "c40a721abc3b96a990994e070339712ead9de840503a34f3b14eb677a0ecdeb6"),
)

PROVIDER_OBJECTS = (
    ("HG02027",0,47590004589,'"7bbd769c0a86a70dc06952494a9ec069-710"',"Mon, 09 Oct 2023 19:37:34 GMT","918gdxbLxtTUKnA16qerH.pqikwAfv3o","s3://human-pangenomics/working/HPRC/HG02027/raw_data/nanopore/guppy_6/12_08_21_R941_HG02027_1_Guppy_6.4.6_450bps_modbases_5mc_cg_sup_prom_pass.bam"),
    ("HG02027",1,60282765028,'"f37d2b282b6f335e82a257107b29adbd-899"',"Mon, 09 Oct 2023 19:38:58 GMT","y20FJriQ9aahhKkkQvchw3vs4DHvU8Oz","s3://human-pangenomics/working/HPRC/HG02027/raw_data/nanopore/guppy_6/12_08_21_R941_HG02027_2_Guppy_6.4.6_450bps_modbases_5mc_cg_sup_prom_pass.bam"),
    ("HG02027",2,59286498259,'"3dea717f177a93e8b5f2dd9f9ffd1755-884"',"Mon, 09 Oct 2023 19:40:32 GMT","Rvmnlx8tyIVq8cG7zZuf7UoHofNo3IZS","s3://human-pangenomics/working/HPRC/HG02027/raw_data/nanopore/guppy_6/12_08_21_R941_HG02027_3_Guppy_6.4.6_450bps_modbases_5mc_cg_sup_prom_pass.bam"),
    ("HG02178",0,53595063895,'"9319ed58c8ba87f63f69dd71f64dcb41-799"',"Thu, 20 Jun 2024 22:27:28 GMT","S9dWMoizJ4OFXlBJiInE.lKFki4tc0JU","s3://human-pangenomics/working/HPRC/HG02178/raw_data/nanopore/dorado0.6.0_sup4.3.0_5mCG_5hmCG/12_12_23_R1041_HPRC_HG02178_2_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam"),
    ("HG02178",1,61982761460,'"9ba12208cd05654573bed02fe2ec5a38-924"',"Thu, 20 Jun 2024 22:28:09 GMT","nzbAGO7CUHR9Cy28iIE32oznb1GktDJj","s3://human-pangenomics/working/HPRC/HG02178/raw_data/nanopore/dorado0.6.0_sup4.3.0_5mCG_5hmCG/12_12_23_R1041_HPRC_HG02178_3_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam"),
    ("HG02178",2,59762798806,'"6b6809e8153e4e430c98a0788eeb078a-891"',"Thu, 20 Jun 2024 22:26:42 GMT","PyjLp7nG12aIK1xXyTiSrXJSz4vce1Je","s3://human-pangenomics/working/HPRC/HG02178/raw_data/nanopore/dorado0.6.0_sup4.3.0_5mCG_5hmCG/12_12_23_R1041_HPRC_HG02178_1_dorado0.6.0_sup4.3.0_5mCG_5hmCG.bam"),
    ("HG03669",0,68854722019,'"438e4963eb79c40a1a0ccf37fced5706-1027"',"Mon, 09 Oct 2023 21:08:55 GMT","GHQWM2Dch2SN7BzOCl9Z8oeu9kKGYfLV","s3://human-pangenomics/working/HPRC/HG03669/raw_data/nanopore/guppy_6/07_27_21_R941_HG03669_1_Guppy_6.5.7_450bps_modbases_5mc_cg_sup_prom_pass.bam"),
    ("HG03669",1,9268899967,'"b321099515e00235bd049639bfb09854-139"',"Mon, 09 Oct 2023 21:10:15 GMT","Y9OueuEniRX9_bVsJiODwhnd7oi2rfjF","s3://human-pangenomics/working/HPRC/HG03669/raw_data/nanopore/guppy_6/07_27_21_R941_HG03669_2_Guppy_6.5.7_450bps_modbases_5mc_cg_sup_prom_pass.bam"),
    ("HG03669",2,69401065162,'"d7aac9b9483f2da69877183b5cc10263-1035"',"Mon, 09 Oct 2023 21:10:41 GMT","IfBlG9Pu27k2uBeL.8sZ6X2Zkd_lHcUD","s3://human-pangenomics/working/HPRC/HG03669/raw_data/nanopore/guppy_6/07_27_21_R941_HG03669_3_Guppy_6.5.7_450bps_modbases_5mc_cg_sup_prom_pass.bam"),
    ("HG03669",3,20285218453,'"6fd90ce729c472033c0f943acbdbbce1-303"',"Mon, 09 Oct 2023 21:12:06 GMT","LoPGmNJ7TXhhVLqLofXolxTBgHJWZNBX","s3://human-pangenomics/working/HPRC/HG03669/raw_data/nanopore/guppy_6/07_27_21_R941_HG03669_4_Guppy_6.5.7_450bps_modbases_5mc_cg_sup_prom_pass.bam"),
    ("HG03834",0,44723785655,'"bfde28623d099c04afd091415cbb52c0-667"',"Tue, 17 Oct 2023 22:31:31 GMT","EQ5NYdMbDRIlAX3.QVyjX0.ub.gW2tyW","s3://human-pangenomics/working/HPRC/HG03834/raw_data/nanopore/guppy_6/02_08_22_R941_HG03834_1_Guppy_6.5.7_450bps_modbases_5mc_cg_sup_prom_pass.bam"),
    ("HG03834",1,45452516869,'"88b3427e15f0199b50063f9ccfc36c91-678"',"Tue, 17 Oct 2023 22:32:46 GMT","UhjJz9gJiqlDM8YYg9UsvfnuPSbAsHHv","s3://human-pangenomics/working/HPRC/HG03834/raw_data/nanopore/guppy_6/02_08_22_R941_HG03834_2_Guppy_6.5.7_450bps_modbases_5mc_cg_sup_prom_pass.bam"),
    ("HG03834",2,42936745540,'"c7d6e921e296975ffc7a9caa97e9d7e0-640"',"Tue, 17 Oct 2023 22:33:57 GMT","oavKXKgP5U9TgO8.pD6mC9iFbYAEUHfM","s3://human-pangenomics/working/HPRC/HG03834/raw_data/nanopore/guppy_6/02_08_22_R941_HG03834_3_Guppy_6.5.7_450bps_modbases_5mc_cg_sup_prom_pass.bam"),
)

RESOURCES = {
    "preparation": {"array": "0-3%2", "cpus": 5, "memory_gib": 64, "wall_hours": 12, "max_concurrent": 2},
    "alignment": {"array": "0-12%2", "cpus": 5, "memory_gib": 128, "wall_hours": 48, "max_concurrent": 2, "dependency": "afterok:preparation"},
    "analysis": {"array": "0-3%2", "cpus": 4, "memory_gib": 32, "wall_hours": 4, "max_concurrent": 2, "dependency": "afterok:alignment"},
}
ACTION_GATES = {
    "cluster_contact_authorized": False,
    "network_authorized": False,
    "deployment_authorized": False,
    "download_authorized": False,
    "submission_authorized": False,
    "execution_authorized": False,
    "cancellation_authorized": False,
    "retry_authorized": False,
    "resubmission_authorized": False,
    "data_access_authorized": False,
    "production_authorized": False,
    "copy_number_inference_authorized": False,
    "result_use_authorized": False,
    "model_use_authorized": False,
    "manuscript_use_authorized": False,
}
PLAN_ACTION_GATES = {
    "cluster_contact_performed": False,
    "network_contact_performed": False,
    "deployment_authorized": False,
    "deployment_performed": False,
    "download_authorized": False,
    "download_performed": False,
    "submission_authorized": False,
    "submission_performed": False,
    "execution_authorized": False,
    "execution_performed": False,
    "cancellation_authorized": False,
    "cancellation_performed": False,
    "retry_authorized": False,
    "retry_performed": False,
    "resubmission_authorized": False,
    "resubmission_performed": False,
    "data_access_authorized": False,
    "production_authorized": False,
    "copy_number_inference_authorized": False,
    "result_use_authorized": False,
    "model_use_authorized": False,
    "manuscript_use_authorized": False,
}

BUNDLE_FILES = (
    ("README.md", "0644", "documentation"),
    ("src/exact4_plan.py", "0755", "planner_and_release_gate"),
    ("src/01_exact4_prepare.sh", "0755", "preparation_entrypoint"),
    ("src/02_exact4_alignment_worker.sh", "0755", "alignment_entrypoint"),
    ("src/03_exact4_analysis.sh", "0755", "analysis_entrypoint"),
    ("src/plot_exact4_depth.py", "0755", "plotter"),
    ("tests/test_exact4_workflow.py", "0644", "fixture_test_suite"),
    ("src/capture_exact4_provider_receipt.py", "0755", "provider_receipt_v2_emitter"),
    ("tests/test_capture_exact4_provider_receipt.py", "0644", "provider_receipt_v2_test_suite"),
    ("frozen/exact4_expected_provider_metadata.v1.json", "0644", "expected_provider_metadata_authority"),
    ("frozen/exact4_aws_tool_authority.v1.json", "0644", "aws_tool_authority"),
)

AMENDED_AUTHORIZATION_PATH = "handoff/7P22_EXACT_FOUR_AMENDED_LIVE_RECEIPT_CAPTURE_AUTHORIZATION.md"
AMENDED_AUTHORIZATION_ACCEPTANCE_PATH = "handoff/7P22_EXACT_FOUR_AMENDED_LIVE_RECEIPT_CAPTURE_AUTHORIZATION_ACCEPTANCE.md"
DEPLOYMENT_MANIFEST_PATH = "cluster_workflows/cnv_7p22_exact4/frozen/exact4_deployment_manifest.v1.json"
AMENDED_READ_ONLY_SCOPE = "read_only_one_aws_version_plus_13_ordered_versioned_head_object_no_sign_request"
PROVIDER_CAPTURE_LABEL = "read_only_exact13_metadata_head:b7d2a267e44291590868e2f679beb4c7f62790f1b89778e804f105266ab04d58"
EXPECTED_PROVIDER_AUTHORITY_SHA256 = "5406e35fed4d44233591b978914b8befbb90d11f6c647584fa6464edca9697ff"
EXPECTED_PROVIDER_ACCEPTANCE_SHA256 = "cd18055697f2433933d2ea154718237d41e70265300c024d5d27d35a24e4117b"
AWS_TOOL_AUTHORITY_SHA256 = "23728802f3249480adb3d45041a75fb38e0c6c52cadb24e4e03580e89bfeba63"
AWS_TOOL_ACCEPTANCE_SHA256 = "59710dc0a77ca92f52b4a531a6e974b3706a1d45649da48fba193896ab994177"
LEGACY_CAPTURE_AUTHORIZATION_SHA256 = "b7d2a267e44291590868e2f679beb4c7f62790f1b89778e804f105266ab04d58"
RESIDENT_REHASH_RECEIPT_SHA256 = "228677eca4e5730bfc89f154404138481d62a94160d0e824e17354d1cbc68b75"
MINIMAP2_RESIDENT_SHA256 = "e8e2e9ccad47e31120a807dea825609df00c1e376066b815690ef347d856f99e"
SAMTOOLS_RESIDENT_SHA256 = "d886be95e9983a2f392d648d35238a16f0a2d06311590c015817f6463d7425ca"
FIXTURE_AUTHORIZATION_PATH = "handoff/7P22_EXACT_FOUR_RESIDENT_TOOL_FIXTURE_EXECUTION_AUTHORIZATION.md"
FIXTURE_AUTHORIZATION_ACCEPTANCE_PATH = "handoff/7P22_EXACT_FOUR_RESIDENT_TOOL_FIXTURE_EXECUTION_AUTHORIZATION_ACCEPTANCE.md"
V2_CORRECTION_ADDENDUM_PATH = "handoff/7P22_EXACT_FOUR_V2_IMPLEMENTATION_REVIEW_CORRECTION_ADDENDUM.md"
V2_CORRECTION_ADDENDUM_SHA256 = "58c55f29ee4eb6bda035f15a60e0c2e12529ca05d1483c12ba27046e44e7a89a"
V2_CORRECTION_ADDENDUM_SIZE = 25825
V2_CORRECTION_ACCEPTANCE_PATH = "handoff/7P22_EXACT_FOUR_V2_IMPLEMENTATION_REVIEW_CORRECTION_ADDENDUM_ACCEPTANCE.md"
V2_CORRECTION_ACCEPTANCE_SHA256 = "1efa58eb0fe121ceeee0d298325cd19f5c2e0d9aac231f442b81764230fe3482"
V2_CORRECTION_ACCEPTANCE_SIZE = 1986
RESIDENT_REHASH_RECEIPT_PATH = "handoff/7P22_EXACT_FOUR_CURRENT_RESIDENT_REHASH_RECEIPT.md"
RESIDENT_REHASH_ACCEPTANCE_PATH = "handoff/7P22_EXACT_FOUR_CURRENT_RESIDENT_REHASH_RECEIPT_ACCEPTANCE.md"
RESIDENT_REHASH_ACCEPTANCE_SHA256 = "47d45db21fd22f7c4c6ecdcb310079113a6cff7c0740107c5c2e95f0ff7e8b93"
AWS_TOOL_AUTHORITY_PATH = "cluster_workflows/cnv_7p22_exact4/frozen/exact4_aws_tool_authority.v1.json"
AWS_TOOL_ACCEPTANCE_PATH = "handoff/7P22_EXACT_FOUR_AWS_TOOL_AUTHORITY_ACCEPTANCE.md"
FIXTURE_READ_ONLY_SCOPE = "offline_exact_two_resident_binaries_seven_toy_fixtures_once"


def sha256_file(path: os.PathLike[str] | str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _duplicates_rejected(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ExactFourPlanError(f"duplicate JSON key: {key!r}")
        result[key] = value
    return result


def _reject_constant(token: str) -> object:
    raise ExactFourPlanError(f"nonfinite JSON token is forbidden: {token}")


def _strict_builtin(value: object, label: str = "JSON") -> None:
    if type(value) in {str, int, bool} or value is None:
        return
    if type(value) is float:
        if not math.isfinite(value):
            raise ExactFourPlanError(f"{label} contains a nonfinite float")
        raise ExactFourPlanError(f"{label} floats are forbidden; use exact scaled integers")
    if type(value) is list:
        for index, item in enumerate(value):
            _strict_builtin(item, f"{label}[{index}]")
        return
    if type(value) is dict:
        for key, item in value.items():
            if type(key) is not str:
                raise ExactFourPlanError(f"{label} object key is not an exact string")
            _strict_builtin(item, f"{label}.{key}")
        return
    raise ExactFourPlanError(f"{label} contains forbidden type {type(value).__name__}")


def strict_json_bytes(payload: bytes, label: str) -> object:
    if type(payload) is not bytes:
        raise ExactFourPlanError(f"{label} must be exact bytes")
    if payload.startswith(b"\xef\xbb\xbf"):
        raise ExactFourPlanError(f"{label} must not contain a UTF-8 BOM")
    try:
        text = payload.decode("utf-8", errors="strict")
        value = json.loads(
            text,
            object_pairs_hook=_duplicates_rejected,
            parse_constant=_reject_constant,
        )
    except ExactFourPlanError:
        raise
    except (UnicodeError, json.JSONDecodeError) as error:
        raise ExactFourPlanError(f"invalid strict JSON in {label}: {error}") from error
    _strict_builtin(value, label)
    return value


def canonical_json_bytes(value: object) -> bytes:
    _strict_builtin(value, "canonical JSON")
    try:
        return (json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode("utf-8")
    except (TypeError, ValueError) as error:
        raise ExactFourPlanError(f"cannot encode canonical JSON: {error}") from error


def digest_value(value: object) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def _require_dict(value: object, label: str) -> dict[str, object]:
    if type(value) is not dict:
        raise ExactFourPlanError(f"{label} must be an exact JSON object")
    return value


def _exact_keys(value: Mapping[str, object], keys: Sequence[str], label: str) -> None:
    observed = set(value)
    expected = set(keys)
    if observed != expected:
        raise ExactFourPlanError(f"{label} keys differ; missing={sorted(expected-observed)}, unknown={sorted(observed-expected)}")


def _require_str(value: object, label: str, *, nonempty: bool = True) -> str:
    if type(value) is not str or (nonempty and not value):
        raise ExactFourPlanError(f"{label} must be an exact{' nonempty' if nonempty else ''} string")
    return value


def _require_int(value: object, label: str, *, minimum: int = 0) -> int:
    if type(value) is not int or value < minimum:
        raise ExactFourPlanError(f"{label} must be an exact integer >= {minimum}")
    return value


def _require_bool(value: object, expected: bool, label: str) -> None:
    if type(value) is not bool or value is not expected:
        raise ExactFourPlanError(f"{label} must be exact built-in {expected}")


def _require_sha(value: object, label: str) -> str:
    text = _require_str(value, label)
    if not SHA256_RE.fullmatch(text):
        raise ExactFourPlanError(f"{label} must be lowercase 64-hex SHA-256")
    return text


def _load_json_file(path: os.PathLike[str] | str, label: str, *, canonical: bool = True) -> dict[str, object]:
    item = Path(path)
    if not item.is_file() or item.is_symlink():
        raise ExactFourPlanError(f"{label} is not a regular non-symlink file: {item}")
    payload = item.read_bytes()
    value = _require_dict(strict_json_bytes(payload, label), label)
    if canonical and canonical_json_bytes(value) != payload:
        raise ExactFourPlanError(f"{label} is not canonical JSON bytes")
    return value


def _fingerprint(path: Path, *, allow_empty: bool = False) -> dict[str, object]:
    if not path.is_file() or path.is_symlink():
        raise ExactFourPlanError(f"required regular file is missing or symlinked: {path}")
    size = path.stat().st_size
    if not allow_empty and size <= 0:
        raise ExactFourPlanError(f"required file is empty: {path}")
    return {"path": str(path), "size_bytes": size, "sha256": sha256_file(path)}


def _project_path_from_relative(relative: str) -> Path:
    """Resolve one fixed project-relative authority path without caller fallback."""
    if type(relative) is not str or not relative or Path(relative).is_absolute() or ".." in Path(relative).parts:
        raise ExactFourPlanError("authority path must be one fixed normalized project-relative path")
    project = Path(__file__).resolve(strict=True).parents[3]
    return project / relative


def _require_lexical_normalized_absolute(path_value: object, label: str) -> Path:
    text = _require_str(path_value, f"{label} path")
    if "\x00" in text or not os.path.isabs(text) or text != os.path.normpath(text) or text.endswith(os.sep) or os.sep * 2 in text:
        raise ExactFourPlanError(f"{label} path must use one lexical normalized absolute spelling")
    if any(part in ("", ".", "..") for part in Path(text).parts[1:]):
        raise ExactFourPlanError(f"{label} path contains a lexical alias")
    return Path(text)


def _load_fixed_deployment_manifest_argument(path_value: object, label: str) -> dict[str, object]:
    """Gate a process-capable public CLI's raw manifest spelling."""
    path = _require_lexical_normalized_absolute(path_value, label)
    expected = _project_path_from_relative(DEPLOYMENT_MANIFEST_PATH)
    if path != expected:
        raise ExactFourPlanError(f"{label} differs from the fixed deployment-manifest role")
    return _validate_deployment_manifest(_load_json_file(str(path), label))


def _stable_reopen_absolute_pinned_file(path_value: object, size_bytes: int, sha256: str, label: str) -> bytes:
    """Stable-read a uniquely-spelled absolute external authorization file."""
    path = _require_lexical_normalized_absolute(path_value, label)
    _require_int(size_bytes, f"{label} size", minimum=1)
    _require_sha(sha256, f"{label} SHA-256")

    leaf_fields = (
        "st_dev", "st_ino", "st_mode", "st_nlink", "st_uid", "st_gid",
        "st_size", "st_mtime_ns", "st_ctime_ns",
    )

    def walk() -> tuple[tuple[object, ...], ...]:
        identities = []
        current = Path(path.anchor)
        for component in (None, *path.parts[1:]):
            if component is not None:
                current /= component
            info = os.lstat(current)
            leaf = current == path
            if stat.S_ISLNK(info.st_mode) or (leaf and not stat.S_ISREG(info.st_mode)) or (not leaf and not stat.S_ISDIR(info.st_mode)):
                raise ExactFourPlanError(f"{label} path walk found unsafe component")
            if leaf:
                identities.append((str(current), *(getattr(info, field) for field in leaf_fields)))
            else:
                # Directory mtime/ctime is child-sensitive.  Bind the stable
                # path/object/type/permission identity so unrelated sibling
                # churn cannot invalidate an unchanged pinned leaf.
                identities.append((str(current), info.st_dev, info.st_ino, info.st_mode))
        return tuple(identities)

    before_walk = walk()
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags)
    except OSError as error:
        raise ExactFourPlanError(f"cannot safely open {label}") from error
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode) or before.st_size != size_bytes:
            raise ExactFourPlanError(f"{label} type/size drift")
        chunks: list[bytes] = []
        remaining = size_bytes + 1
        while remaining:
            chunk = os.read(descriptor, min(1024 * 1024, remaining))
            if not chunk:
                break
            chunks.append(chunk)
            remaining -= len(chunk)
        payload = b"".join(chunks)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if any(getattr(before, field) != getattr(after, field) for field in leaf_fields):
        raise ExactFourPlanError(f"{label} changed during stable read")
    if len(payload) != size_bytes or hashlib.sha256(payload).hexdigest() != sha256:
        raise ExactFourPlanError(f"{label} content drift")
    if before_walk != walk():
        raise ExactFourPlanError(f"{label} root/ancestor/leaf identity drift")
    reopened = os.open(path, flags)
    try:
        reopened_info = os.fstat(reopened)
    finally:
        os.close(reopened)
    if any(getattr(reopened_info, field) != getattr(before, field) for field in leaf_fields):
        raise ExactFourPlanError(f"{label} changed before post-read reopen")
    return payload


def _stable_reopen_pinned_file(relative: str, size_bytes: int, sha256: str, label: str) -> bytes:
    """Stable-read then reopen/re-walk root, ancestors and leaf.

    The second walk deliberately occurs after the descriptor read so a renamed
    ancestor or replaced leaf cannot be hidden by a still-valid descriptor.
    """
    _require_int(size_bytes, f"{label} size", minimum=1)
    _require_sha(sha256, f"{label} SHA-256")
    path = _project_path_from_relative(relative)

    leaf_fields = (
        "st_dev", "st_ino", "st_mode", "st_nlink", "st_uid", "st_gid",
        "st_size", "st_mtime_ns", "st_ctime_ns",
    )

    def walk() -> tuple[tuple[object, ...], ...]:
        identities = []
        current = Path(path.anchor)
        root_info = os.lstat(current)
        if stat.S_ISLNK(root_info.st_mode) or not stat.S_ISDIR(root_info.st_mode):
            raise ExactFourPlanError(f"{label} filesystem root is unsafe")
        identities.append((str(current), root_info.st_dev, root_info.st_ino, root_info.st_mode))
        for component in path.parts[1:]:
            current /= component
            info = os.lstat(current)
            is_leaf = current == path
            if stat.S_ISLNK(info.st_mode) or (is_leaf and not stat.S_ISREG(info.st_mode)) or (not is_leaf and not stat.S_ISDIR(info.st_mode)):
                raise ExactFourPlanError(f"{label} path component is unsafe: {current}")
            if is_leaf:
                identities.append((str(current), *(getattr(info, field) for field in leaf_fields)))
            else:
                identities.append((str(current), info.st_dev, info.st_ino, info.st_mode))
        return tuple(identities)

    before_walk = walk()
    flags = os.O_RDONLY | (getattr(os, "O_NOFOLLOW", 0))
    try:
        descriptor = os.open(path, flags)
    except OSError as error:
        raise ExactFourPlanError(f"cannot safely open {label}") from error
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode) or before.st_size != size_bytes:
            raise ExactFourPlanError(f"{label} type/size drift")
        chunks = []
        remaining = size_bytes + 1
        while remaining:
            chunk = os.read(descriptor, min(1024 * 1024, remaining))
            if not chunk:
                break
            chunks.append(chunk)
            remaining -= len(chunk)
        payload = b"".join(chunks)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if any(getattr(before, field) != getattr(after, field) for field in leaf_fields):
        raise ExactFourPlanError(f"{label} changed during stable read")
    if len(payload) != size_bytes or hashlib.sha256(payload).hexdigest() != sha256:
        raise ExactFourPlanError(f"{label} content drift")
    after_walk = walk()
    if before_walk != after_walk:
        raise ExactFourPlanError(f"{label} root/ancestor/leaf identity drift")
    try:
        reopened = os.open(path, flags)
    except OSError as error:
        raise ExactFourPlanError(f"cannot reopen {label}") from error
    try:
        reopened_info = os.fstat(reopened)
    finally:
        os.close(reopened)
    if any(getattr(reopened_info, field) != getattr(before, field) for field in leaf_fields):
        raise ExactFourPlanError(f"{label} changed before post-read reopen")
    return payload


def _project_root(bundle_root: Path) -> Path:
    root = bundle_root.resolve(strict=True)
    if root.name != "cnv_7p22_exact4" or root.parent.name != "cluster_workflows":
        raise ExactFourPlanError("bundle_root must be the exact cnv_7p22_exact4 directory")
    return root.parents[1]


def _provider_manifest_rows() -> list[dict[str, object]]:
    return [
        {"sample_id": sample, "part_index": part, "content_length_bytes": length, "etag": etag,
         "last_modified": modified, "version_id": version, "uri": uri}
        for sample, part, length, etag, modified, version, uri in PROVIDER_OBJECTS
    ]


def _assembly_manifest_rows() -> list[dict[str, object]]:
    return [
        {"sample_id": sample, "haplotype": hap, "fasta_path": fasta,
         "fasta_sha256": fsha, "fai_path": fasta + ".fai", "fai_sha256": isha}
        for sample, hap, fasta, fsha, isha in ASSEMBLY_AUTHORITIES
    ]


def _manifest_identity(manifest: Mapping[str, object]) -> dict[str, object]:
    return {key: value for key, value in manifest.items() if key != "bundle_id"}


def _validate_deployment_manifest(value: object) -> dict[str, object]:
    manifest = _require_dict(value, "deployment manifest")
    keys = (
        "schema_version", "bundle_id", "bundle_files", "prerequisite_pins",
        "classification_proof_sha256", "target_table", "resident_authorities",
        "assembly_authorities", "expected_provider_objects", "targets", "dependencies",
        "output_schemas", "resources", "scientific_laws", "coverage_ratio_schema_identity",
        "coverage_ratio_schema_sha256", *CONTRACT_BINDINGS.keys(),
        "action_gates",
    )
    _exact_keys(manifest, keys, "deployment manifest")
    if manifest["schema_version"] != DEPLOYMENT_SCHEMA:
        raise ExactFourPlanError("deployment manifest schema drift")
    bundle_files = manifest["bundle_files"]
    if type(bundle_files) is not list or len(bundle_files) != 11:
        raise ExactFourPlanError("amended bundle must contain exactly 11 listed files")
    for row_value, (path, mode, role) in zip(bundle_files, BUNDLE_FILES):
        row = _require_dict(row_value, "bundle file row")
        _exact_keys(row, ("path", "size_bytes", "mode", "sha256", "role"), "bundle file row")
        if (row["path"], row["mode"], row["role"]) != (path, mode, role):
            raise ExactFourPlanError("amended bundle path/mode/role order drift")
        _require_int(row["size_bytes"], "bundle file size", minimum=1)
        _require_sha(row["sha256"], "bundle file SHA-256")
    pins = manifest["prerequisite_pins"]
    if type(pins) is not list or pins != [{"path": path, "sha256": sha} for path, sha in PREREQUISITE_PINS]:
        raise ExactFourPlanError("amended prerequisite pin set/order drift")
    if manifest["target_table"] != {"bytes": 100, "sha256": TARGET_TABLE_SHA256, "text": TARGET_TABLE_BYTES.decode()}:
        raise ExactFourPlanError("target authority drift")
    if manifest["targets"] != [
        {"sample_id": s, "locus": LOCUS, "part_count": c} for s, c in zip(TARGETS, PART_COUNTS)
    ]:
        raise ExactFourPlanError("exact-four target/order/count drift")
    if manifest["expected_provider_objects"] != _provider_manifest_rows():
        raise ExactFourPlanError("expected provider object identity drift")
    if manifest["assembly_authorities"] != _assembly_manifest_rows():
        raise ExactFourPlanError("assembly authority drift")
    if manifest["output_schemas"] != OUTPUT_SCHEMAS or manifest["resources"] != RESOURCES:
        raise ExactFourPlanError("output schema or resource topology drift")
    laws = _require_dict(manifest["scientific_laws"], "scientific laws")
    if laws.get("minimap2_version") != "2.26-r1175" or laws.get("samtools_version") != "1.21":
        raise ExactFourPlanError("resident tool scientific-law repin drift")
    gates = _require_dict(manifest["action_gates"], "deployment action_gates")
    if gates != ACTION_GATES or any(type(gates[key]) is not bool or gates[key] is not False for key in ACTION_GATES):
        raise ExactFourPlanError("deployment action gates must remain exactly false")
    if manifest["coverage_ratio_schema_identity"] != COVERAGE_RATIO_SCHEMA_IDENTITY or manifest["coverage_ratio_schema_sha256"] != COVERAGE_RATIO_SCHEMA_SHA256:
        raise ExactFourPlanError("coverage-ratio schema identity drift")
    for key, expected in CONTRACT_BINDINGS.items():
        if manifest[key] != expected:
            raise ExactFourPlanError(f"deployment contract binding drift: {key}")
    schema_bytes = canonical_json_bytes(manifest["coverage_ratio_schema_identity"])
    if len(schema_bytes) != 2549 or hashlib.sha256(schema_bytes).hexdigest() != COVERAGE_RATIO_SCHEMA_SHA256:
        raise ExactFourPlanError("coverage-ratio schema byte identity drift")
    forbidden_manifest_text = canonical_json_bytes(manifest).lower()
    for phrase in (b"coverage-compatible", b"coverage-discordant", b"ratio_min", b"ratio_max", b"compatibility_tolerance"):
        if phrase in forbidden_manifest_text:
            raise ExactFourPlanError("forbidden ratio-derived category/tolerance entered manifest")
    expected_id = "exact4-bundle-" + digest_value(_manifest_identity(manifest))[:24]
    if manifest["bundle_id"] != expected_id:
        raise ExactFourPlanError("bundle_id does not match canonical bundle identity")
    return manifest


def build_deployment_manifest(bundle_root: pathlib.Path) -> dict[str, object]:
    """Build deterministic bundle identity from the eleven local files.

    The frozen manifest is deliberately not part of its own identity.
    """
    root = Path(bundle_root)
    project = _project_root(root)
    files = []
    for relative, mode, role in BUNDLE_FILES:
        path = root / relative
        fingerprint = _fingerprint(path)
        observed_mode = stat.S_IMODE(path.stat().st_mode)
        expected_mode = int(mode, 8)
        if observed_mode != expected_mode:
            raise ExactFourPlanError(f"bundle file mode drift for {relative}: expected {mode}, observed {observed_mode:04o}")
        files.append({"path": relative, "size_bytes": fingerprint["size_bytes"], "mode": mode, "sha256": fingerprint["sha256"], "role": role})
    pins = []
    for relative, expected in PREREQUISITE_PINS:
        path = project / relative
        observed = sha256_file(path)
        if observed != expected:
            raise ExactFourPlanError(f"prerequisite pin drift for {relative}: {observed}")
        pins.append({"path": relative, "sha256": expected})
    generic = project / "cluster_workflows/cnv_v2/src/cnv_plan.py"
    generic_sha = sha256_file(generic)
    if generic_sha != "392a9a8bc63e323c06f0abdb0d013e190082a7b0fed0fbfaba7f42480d75363c":
        raise ExactFourPlanError("frozen generic parity oracle drift")
    manifest: dict[str, object] = {
        "schema_version": DEPLOYMENT_SCHEMA,
        "bundle_id": "",
        "bundle_files": files,
        "prerequisite_pins": pins,
        "classification_proof_sha256": CLASSIFICATION_PROOF_SHA256,
        "target_table": {"bytes": len(TARGET_TABLE_BYTES), "sha256": TARGET_TABLE_SHA256, "text": TARGET_TABLE_BYTES.decode()},
        "resident_authorities": [
            {"role": role, "path": path, "size_bytes": size, "sha256": digest}
            for role, path, size, digest in RESIDENT_AUTHORITIES
        ],
        "assembly_authorities": _assembly_manifest_rows(),
        "expected_provider_objects": _provider_manifest_rows(),
        "targets": [{"sample_id": s, "locus": LOCUS, "part_count": c} for s, c in zip(TARGETS, PART_COUNTS)],
        "dependencies": [{"path": "cluster_workflows/cnv_v2/src/cnv_plan.py", "sha256": generic_sha, "role": "test_only_parity_oracle"}],
        "output_schemas": dict(OUTPUT_SCHEMAS),
        "resources": json.loads(json.dumps(RESOURCES)),
        "scientific_laws": {
            "diagnostic_only": True,
            "read_inferred_copy_number": "not_estimated",
            "combined_reference": True,
            "complete_unfiltered_alignment": True,
            "primary_exclusion_mask_hex": "0x900",
            "primary_exclusion_mask_decimal": 2304,
            "depth_include_flags_hex": "0x600",
            "depth_include_flags_decimal": 1536,
            "raw_depth_minimum_mapq": 0,
            "mapq10_depth_minimum_mapq": 10,
            "minimum_callable_flank_mean": {"numerator": 5, "denominator": 1},
            "core_expansion_bp": 1000,
            "full_window_additional_expansion_bp": 10000,
            "minimap2_version": "2.26-r1175",
            "samtools_version": "1.21",
            "decode_argv": ["samtools", "fasta", "-F", "0x900"],
            "alignment_argv": ["minimap2", "-ax", "map-ont", "--secondary=no"],
            "kcon_argv": ["minimap2", "-c", "-p", "0.1", "-N", "10"],
            "kcon_query_name": "type2_KCON",
            "kcon_query_length": 9472,
            "kcon_minimum_block_bp": 500,
        },
        "coverage_ratio_schema_identity": copy.deepcopy(COVERAGE_RATIO_SCHEMA_IDENTITY),
        "coverage_ratio_schema_sha256": COVERAGE_RATIO_SCHEMA_SHA256,
        **CONTRACT_BINDINGS,
        "action_gates": dict(ACTION_GATES),
    }
    manifest["bundle_id"] = "exact4-bundle-" + digest_value(_manifest_identity(manifest))[:24]
    return _validate_deployment_manifest(manifest)


def _read_fai(path: Path) -> dict[str, int]:
    records: dict[str, int] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, 1):
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 5 or not fields[0] or not fields[1].isdigit() or int(fields[1]) < 1:
                raise ExactFourPlanError(f"malformed FAI row {path}:{line_number}")
            if fields[0] in records:
                raise ExactFourPlanError(f"duplicate FAI contig {fields[0]!r}: {path}")
            records[fields[0]] = int(fields[1])
    if not records:
        raise ExactFourPlanError(f"empty FAI: {path}")
    return records


def _strict_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t", restkey="__extra__", restval=None, strict=True)
            header = list(reader.fieldnames or [])
            if not header or any(not field for field in header) or len(header) != len(set(header)):
                raise ExactFourPlanError(f"malformed or duplicate TSV header: {path}")
            rows: list[dict[str, str]] = []
            for number, row in enumerate(reader, 2):
                if row.get("__extra__") is not None or any(row.get(field) is None for field in header):
                    raise ExactFourPlanError(f"malformed TSV row {path}:{number}")
                rows.append({field: str(row[field]) for field in header})
            return header, rows
    except (OSError, UnicodeError, csv.Error) as error:
        raise ExactFourPlanError(f"cannot parse strict TSV {path}: {error}") from error


def _merge_intervals(intervals: Sequence[tuple[int, int]]) -> list[tuple[int, int]]:
    ordered = sorted(intervals)
    output: list[list[int]] = []
    for start, end in ordered:
        if type(start) is not int or type(end) is not int or start < 0 or start >= end:
            raise ExactFourPlanError("invalid half-open interval")
        if output and start <= output[-1][1]:
            output[-1][1] = max(output[-1][1], end)
        else:
            output.append([start, end])
    return [(start, end) for start, end in output]


def _subtract_intervals(full: Sequence[tuple[int, int]], core: Sequence[tuple[int, int]]) -> list[tuple[int, int]]:
    result: list[tuple[int, int]] = []
    merged_core = _merge_intervals(core)
    for full_start, full_end in _merge_intervals(full):
        cursor = full_start
        for core_start, core_end in merged_core:
            if core_end <= cursor or core_start >= full_end:
                continue
            if cursor < core_start:
                result.append((cursor, min(core_start, full_end)))
            cursor = max(cursor, core_end)
            if cursor >= full_end:
                break
        if cursor < full_end:
            result.append((cursor, full_end))
    return result


def _tool_command(path: str, argv: Sequence[str], label: str) -> dict[str, object]:
    executable = Path(path)
    if not executable.is_absolute() or not executable.is_file() or executable.is_symlink() or not os.access(executable, os.X_OK):
        raise ExactFourPlanError(f"{label} must be an absolute executable regular file")
    command = [str(executable), *argv]
    try:
        completed = subprocess.run(command, check=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=False, timeout=30)
    except (OSError, subprocess.SubprocessError) as error:
        raise ExactFourPlanError(f"cannot capture {label} tool evidence: {error}") from error
    if completed.returncode not in {0, 1, 2}:
        raise ExactFourPlanError(f"{label} evidence command returned unexpected status {completed.returncode}")
    output = completed.stdout.decode("utf-8", errors="strict")
    if not output:
        raise ExactFourPlanError(f"{label} evidence command returned empty output")
    return {"argv": command, "output": output, "output_sha256": hashlib.sha256(completed.stdout).hexdigest(), "returncode": completed.returncode}


def _validate_tool_receipts(tools: object) -> list[dict[str, object]]:
    if type(tools) is not list or len(tools) != 6:
        raise ExactFourPlanError("resident tool receipt must contain exactly six tools")
    expected = ("bash", "python", "aws", "minimap2", "samtools", "sbatch")
    version_argv = {"bash": ["--version"], "python": ["--version"], "aws": ["--version"], "minimap2": ["--version"], "samtools": ["--version"], "sbatch": ["--version"]}
    help_argv = {"bash": ["--help"], "python": ["--help"], "aws": ["--help"], "minimap2": ["--help"], "samtools": ["depth", "--help"], "sbatch": ["--help"]}
    output: list[dict[str, object]] = []
    for name, value in zip(expected, tools):
        row = _require_dict(value, f"tool receipt {name}")
        _exact_keys(row, ("name", "path", "size_bytes", "sha256", "version", "help"), f"tool receipt {name}")
        if row["name"] != name or not Path(_require_str(row["path"], f"{name} path")).is_absolute():
            raise ExactFourPlanError(f"tool order/path drift for {name}")
        _require_int(row["size_bytes"], f"{name} size", minimum=1)
        _require_sha(row["sha256"], f"{name} SHA-256")
        for evidence_name in ("version", "help"):
            evidence = _require_dict(row[evidence_name], f"{name} {evidence_name}")
            _exact_keys(evidence, ("argv", "output", "output_sha256", "returncode"), f"{name} {evidence_name}")
            output_text = _require_str(evidence["output"], f"{name} {evidence_name} output")
            if _require_sha(evidence["output_sha256"], f"{name} {evidence_name} digest") != hashlib.sha256(output_text.encode("utf-8")).hexdigest():
                raise ExactFourPlanError(f"{name} {evidence_name} output digest drift")
            expected_argv = [row["path"], *(version_argv[name] if evidence_name == "version" else help_argv[name])]
            if evidence["argv"] != expected_argv:
                raise ExactFourPlanError(f"{name} {evidence_name} argv drift")
            if type(evidence["returncode"]) is not int or evidence["returncode"] not in {0, 1, 2}:
                raise ExactFourPlanError(f"{name} {evidence_name} returncode drift")
        output.append(row)
    minimap = output[3]["version"]["output"]
    samtools_version = output[4]["version"]["output"]
    samtools_help = output[4]["help"]["output"]
    if "2.26-r1175" not in minimap or output[3]["sha256"] != MINIMAP2_RESIDENT_SHA256:
        raise ExactFourPlanError("minimap2 must be exact resident 2.26-r1175 bytes")
    if not re.search(r"(?:samtools\s+)?1\.21(?:\s|$)", samtools_version) or output[4]["sha256"] != SAMTOOLS_RESIDENT_SHA256:
        raise ExactFourPlanError("samtools must be exact resident 1.21 bytes")
    if output[3]["version"]["returncode"] != 0 or output[3]["help"]["returncode"] != 0 or output[4]["version"]["returncode"] != 0 or output[4]["help"]["returncode"] != 0:
        raise ExactFourPlanError("resident minimap2/samtools version/help commands must return zero")
    entries = [line.strip() for line in samtools_help.splitlines() if line.strip()]
    if not any(re.search(r"(?:^|\s)-Q(?:[,\s]|$)", line) and any(token in line.lower() for token in ("min-mq", "min mapping quality", "minimum mapping quality")) for line in entries):
        raise ExactFourPlanError("samtools depth help does not bind uppercase -Q to minimum mapping quality")
    if not any(re.search(r"(?:^|\s)-g(?:[,\s]|$)", line) and any(token in line.lower() for token in ("include", "add", "flag")) for line in entries):
        raise ExactFourPlanError("samtools depth help does not bind -g to explicit flag inclusion")
    return output


FIXTURE_IDS = (
    "combined_reference_preparation",
    "primary_source_decoding",
    "competitive_ont_alignment",
    "full_window_primary_bam",
    "part_merge_and_ownership",
    "raw_and_mapq10_depth_laws",
    "kcon_annotation_only_paf",
)


def _validate_fixture_authorization_files(
    authorization_path: object,
    authorization_size_bytes: int,
    authorization_sha256: str,
    acceptance_path: object,
    acceptance_size_bytes: int,
    acceptance_sha256: str,
    tools: object,
    manifest: object,
) -> tuple[dict[str, object], dict[str, object]]:
    """Authenticate the two external one-shot fixture records without cycles."""
    auth_path = _require_lexical_normalized_absolute(authorization_path, "fixture authorization")
    accept_path = _require_lexical_normalized_absolute(acceptance_path, "fixture authorization acceptance")
    if auth_path != _project_path_from_relative(FIXTURE_AUTHORIZATION_PATH) or accept_path != _project_path_from_relative(FIXTURE_AUTHORIZATION_ACCEPTANCE_PATH):
        raise ExactFourPlanError("fixture authorization path-role drift")
    auth_payload = _stable_reopen_absolute_pinned_file(str(auth_path), authorization_size_bytes, authorization_sha256, "fixture authorization")
    acceptance_payload = _stable_reopen_absolute_pinned_file(str(accept_path), acceptance_size_bytes, acceptance_sha256, "fixture authorization acceptance")
    auth = _require_dict(strict_json_bytes(auth_payload, "fixture authorization"), "fixture authorization")
    acceptance = _require_dict(strict_json_bytes(acceptance_payload, "fixture authorization acceptance"), "fixture authorization acceptance")
    if canonical_json_bytes(auth) != auth_payload or canonical_json_bytes(acceptance) != acceptance_payload:
        raise ExactFourPlanError("fixture authorization records must be canonical JSON bytes")
    _exact_keys(acceptance, ("schema_version", "decision", "authorization"), "fixture authorization acceptance")
    if acceptance["schema_version"] != TOOL_FIXTURE_AUTH_ACCEPTANCE_SCHEMA or acceptance["decision"] != "ACCEPT":
        raise ExactFourPlanError("fixture authorization is not independently accepted")
    accepted_binding = _require_dict(acceptance["authorization"], "accepted fixture authorization binding")
    _exact_keys(accepted_binding, ("path", "size_bytes", "sha256"), "accepted fixture authorization binding")
    if accepted_binding != {"path": str(auth_path), "size_bytes": authorization_size_bytes, "sha256": authorization_sha256}:
        raise ExactFourPlanError("fixture authorization acceptance binds different bytes")

    keys = (
        "schema_version", "deployment_bundle_id", "deployment_manifest",
        "correction_addendum", "correction_acceptance", "host_class",
        "tool_bindings", "accepted_roots", "fixture_order", "scratch_root",
        "read_only_action_scope", "artifact_policy", "action_gates",
    )
    _exact_keys(auth, keys, "fixture authorization")
    if auth["schema_version"] != TOOL_FIXTURE_AUTH_SCHEMA or auth["host_class"] != "login.pax.tufts.edu":
        raise ExactFourPlanError("fixture authorization schema/host drift")
    accepted_manifest = _validate_deployment_manifest(manifest)
    manifest_binding = _require_dict(auth["deployment_manifest"], "fixture manifest binding")
    _exact_keys(manifest_binding, ("path", "size_bytes", "sha256"), "fixture manifest binding")
    manifest_path = _require_lexical_normalized_absolute(manifest_binding["path"], "fixture manifest")
    if manifest_path != _project_path_from_relative(DEPLOYMENT_MANIFEST_PATH):
        raise ExactFourPlanError("fixture manifest path-role drift")
    manifest_payload = _stable_reopen_absolute_pinned_file(str(manifest_path), _require_int(manifest_binding["size_bytes"], "fixture manifest size", minimum=1), _require_sha(manifest_binding["sha256"], "fixture manifest SHA-256"), "fixture manifest")
    if manifest_payload != canonical_json_bytes(accepted_manifest) or auth["deployment_bundle_id"] != accepted_manifest["bundle_id"]:
        raise ExactFourPlanError("fixture authorization manifest/bundle drift")

    for key, path, size, digest in (
        ("correction_addendum", V2_CORRECTION_ADDENDUM_PATH, V2_CORRECTION_ADDENDUM_SIZE, V2_CORRECTION_ADDENDUM_SHA256),
        ("correction_acceptance", V2_CORRECTION_ACCEPTANCE_PATH, V2_CORRECTION_ACCEPTANCE_SIZE, V2_CORRECTION_ACCEPTANCE_SHA256),
    ):
        binding = _require_dict(auth[key], f"fixture {key} binding")
        _exact_keys(binding, ("path", "size_bytes", "sha256"), f"fixture {key} binding")
        if binding != {"path": path, "size_bytes": size, "sha256": digest}:
            raise ExactFourPlanError(f"fixture {key} identity drift")
        _stable_reopen_pinned_file(path, size, digest, f"fixture {key}")

    accepted_tools = _validate_tool_receipts(tools)
    by_name = {row["name"]: row for row in accepted_tools}
    tool_bindings = _require_dict(auth["tool_bindings"], "fixture authorization tool bindings")
    _exact_keys(tool_bindings, ("minimap2", "samtools"), "fixture authorization tool bindings")
    for name in ("minimap2", "samtools"):
        expected = {key: by_name[name][key] for key in ("path", "size_bytes", "sha256", "version", "help")}
        if tool_bindings[name] != expected:
            raise ExactFourPlanError(f"fixture authorization {name} binding drift")

    roots = _require_dict(auth["accepted_roots"], "fixture accepted roots")
    expected_roots = {
        "resident_rehash_receipt": {"path": RESIDENT_REHASH_RECEIPT_PATH, "size_bytes": 4571, "sha256": RESIDENT_REHASH_RECEIPT_SHA256},
        "resident_rehash_acceptance": {"path": RESIDENT_REHASH_ACCEPTANCE_PATH, "size_bytes": 1929, "sha256": RESIDENT_REHASH_ACCEPTANCE_SHA256},
        "aws_tool_authority": {"path": AWS_TOOL_AUTHORITY_PATH, "size_bytes": 1165, "sha256": AWS_TOOL_AUTHORITY_SHA256},
        "aws_tool_acceptance": {"path": AWS_TOOL_ACCEPTANCE_PATH, "size_bytes": 1855, "sha256": AWS_TOOL_ACCEPTANCE_SHA256},
    }
    if roots != expected_roots:
        raise ExactFourPlanError("fixture authorization accepted-root identity drift")
    for role, binding in expected_roots.items():
        _stable_reopen_pinned_file(binding["path"], binding["size_bytes"], binding["sha256"], f"fixture {role}")
    if auth["fixture_order"] != list(FIXTURE_IDS) or auth["read_only_action_scope"] != FIXTURE_READ_ONLY_SCOPE:
        raise ExactFourPlanError("fixture authorization order/scope drift")
    _require_lexical_normalized_absolute(auth["scratch_root"], "fixture scratch root")
    if auth["artifact_policy"] != {
        "create_only_local_fixture_artifacts": True,
        "reuse_authorized": False,
        "retry_loop_authorized": False,
    }:
        raise ExactFourPlanError("fixture authorization artifact policy drift")
    gates = _require_dict(auth["action_gates"], "fixture authorization action gates")
    _exact_keys(gates, tuple(ACTION_GATES), "fixture authorization action gates")
    if gates != ACTION_GATES or any(type(gates[key]) is not bool or gates[key] is not False for key in ACTION_GATES):
        raise ExactFourPlanError("fixture authorization external/use gates are not all false")
    host = os.uname().nodename
    if type(host) is not str or re.fullmatch(r"login-[0-9]+\.cluster\.example", host) is None:
        raise ExactFourPlanError("fixture runtime is not an accepted login host")
    return auth, {
        "authorization_path": str(auth_path), "authorization_size_bytes": authorization_size_bytes,
        "authorization_sha256": authorization_sha256, "acceptance_path": str(accept_path),
        "acceptance_size_bytes": acceptance_size_bytes, "acceptance_sha256": acceptance_sha256,
    }


def _fixture_manifest_and_authorization(receipt: Mapping[str, object], tools: object) -> tuple[dict[str, object], dict[str, object]]:
    manifest_binding = _require_dict(receipt["deployment_manifest"], "fixture receipt manifest binding")
    _exact_keys(manifest_binding, ("path", "size_bytes", "sha256"), "fixture receipt manifest binding")
    manifest_path = _require_lexical_normalized_absolute(manifest_binding["path"], "fixture receipt manifest")
    if manifest_path != _project_path_from_relative(DEPLOYMENT_MANIFEST_PATH):
        raise ExactFourPlanError("fixture receipt manifest path-role drift")
    payload = _stable_reopen_absolute_pinned_file(str(manifest_path), _require_int(manifest_binding["size_bytes"], "fixture receipt manifest size", minimum=1), _require_sha(manifest_binding["sha256"], "fixture receipt manifest SHA-256"), "fixture receipt manifest")
    manifest = _validate_deployment_manifest(strict_json_bytes(payload, "fixture receipt manifest"))
    if canonical_json_bytes(manifest) != payload:
        raise ExactFourPlanError("fixture receipt manifest is not canonical")
    binding = _require_dict(receipt["authorization_binding"], "fixture receipt authorization binding")
    _exact_keys(binding, ("authorization_path", "authorization_size_bytes", "authorization_sha256", "acceptance_path", "acceptance_size_bytes", "acceptance_sha256"), "fixture receipt authorization binding")
    _, normalized = _validate_fixture_authorization_files(
        binding["authorization_path"], binding["authorization_size_bytes"], binding["authorization_sha256"],
        binding["acceptance_path"], binding["acceptance_size_bytes"], binding["acceptance_sha256"],
        tools, manifest,
    )
    if binding != normalized:
        raise ExactFourPlanError("fixture receipt authorization binding drift")
    return manifest, normalized


def _fixture_path_map(fixture: Mapping[str, object]) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for collection in ("input_files", "output_files"):
        for row in fixture[collection]:
            path = Path(str(row["path"]))
            if path.name in result:
                raise ExactFourPlanError("fixture artifact basenames must be unique")
            result[path.name] = path
    return result


def _parse_toy_fasta(path: Path) -> list[tuple[bytes, bytes]]:
    payload = path.read_bytes()
    if not payload or not payload.endswith(b"\n"):
        raise ExactFourPlanError("fixture FASTA must be nonempty and LF terminated")
    records: list[tuple[bytes, bytes]] = []
    name: bytes | None = None
    sequence = bytearray()
    for line in payload.splitlines():
        if line.startswith(b">"):
            if name is not None:
                records.append((name, bytes(sequence)))
            name, sequence = line[1:], bytearray()
            if not name or any(byte <= 32 or byte >= 127 for byte in name):
                raise ExactFourPlanError("fixture FASTA name is malformed")
        else:
            if name is None or not line or any(base not in b"ACGTN" for base in line.upper()):
                raise ExactFourPlanError("fixture FASTA sequence is malformed")
            sequence.extend(line.upper())
    if name is None:
        raise ExactFourPlanError("fixture FASTA has no records")
    records.append((name, bytes(sequence)))
    if len({item[0] for item in records}) != len(records):
        raise ExactFourPlanError("fixture FASTA names are duplicated")
    return records


def _parse_toy_sam(path: Path) -> list[list[bytes]]:
    payload = path.read_bytes()
    if not payload or not payload.endswith(b"\n"):
        raise ExactFourPlanError("fixture SAM must be nonempty and LF terminated")
    records: list[list[bytes]] = []
    for line in payload.splitlines():
        if line.startswith(b"@"):
            continue
        fields = line.split(b"\t")
        if len(fields) < 11 or not fields[0]:
            raise ExactFourPlanError("fixture SAM row is malformed")
        try:
            int(fields[1]); int(fields[3]); int(fields[4])
        except ValueError as error:
            raise ExactFourPlanError("fixture SAM numeric field is malformed") from error
        records.append(fields)
    return records


def _fixture_depth_values(path: Path) -> list[int]:
    rows = [line.split(b"\t") for line in path.read_bytes().splitlines()]
    if len(rows) != 20 or any(len(row) != 3 for row in rows) or [(row[0], int(row[1])) for row in rows] != [(b"hapA", index) for index in range(1, 21)]:
        raise ExactFourPlanError("fixture depth coordinates are not exact")
    return [int(row[2]) for row in rows]


def _load_manifest_pinned_plotter(manifest: Mapping[str, object]) -> object:
    rows = [row for row in manifest["bundle_files"] if row["path"] == "src/plot_exact4_depth.py"]
    if len(rows) != 1:
        raise ExactFourPlanError("manifest lacks one exact plotter identity")
    path = _project_path_from_relative("cluster_workflows/cnv_7p22_exact4/src/plot_exact4_depth.py")
    observed = _fingerprint(path)
    if observed["size_bytes"] != rows[0]["size_bytes"] or observed["sha256"] != rows[0]["sha256"]:
        raise ExactFourPlanError("loaded plotter differs from manifest-pinned bytes")
    spec = importlib.util.spec_from_file_location("_exact4_manifest_pinned_plotter", path)
    if spec is None or spec.loader is None or spec.origin != str(path):
        raise ExactFourPlanError("cannot create manifest-pinned plotter import")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if Path(module.__file__) != path:
        raise ExactFourPlanError("manifest-pinned plotter origin drift")
    return module


def _independent_kcon_oracle(authority: Mapping[str, object], paf_rows: list[tuple[list[str], bytes]], paf_sha256: str) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    for element in authority["element_identities"]:
        windows = [
            (row["contig"], start, end)
            for row in authority["bed_authorities"]
            if row["haplotype"] == element["haplotype"] and row["contig"] == element["contig"]
            for start, end in row["full_intervals"]
            if start <= element["source_start"] and element["source_end"] <= end
        ]
        if len(windows) != 1:
            raise ExactFourPlanError("independent KCON oracle lacks one containing window")
        contig, window_start, window_end = windows[0]
        target = f"{contig}:{window_start + 1}-{window_end}"
        qualified = []
        for fields, raw in paf_rows:
            if fields[0] != "type2_KCON" or int(fields[1]) != 9472 or fields[5] != target or int(fields[6]) != window_end - window_start:
                continue
            absolute_start, absolute_end = window_start + int(fields[7]), window_start + int(fields[8])
            if int(fields[10]) >= 500 and absolute_start < element["source_end"] and absolute_end > element["source_start"]:
                qualified.append((fields, raw))
        if not qualified:
            raise ExactFourPlanError("independent KCON oracle has no qualifying row")
        fields, raw = min(qualified, key=lambda item: (-int(item[0][10]), -int(item[0][9]), -int(item[0][11]), int(item[0][7]), int(item[0][2]), item[1]))
        result.append({
            "sample_id": authority["sample_id"], "haplotype_id": element["haplotype"], "contig": contig,
            "source_start_0based": element["source_start"], "source_end_0based": element["source_end"],
            "query_name": fields[0], "query_length": int(fields[1]), "query_start_0based": int(fields[2]),
            "query_end_0based": int(fields[3]), "strand": fields[4], "target_name": fields[5],
            "target_length": int(fields[6]), "target_start_0based": int(fields[7]), "target_end_0based": int(fields[8]),
            "matches": int(fields[9]), "aligned_block_length": int(fields[10]), "mapq": int(fields[11]),
            "paf_row_sha256": hashlib.sha256(raw).hexdigest(), "paf_file_sha256": paf_sha256,
            "annotation_only": True, "read_inferred_copy_number": "not_estimated",
        })
    return result


def _recompute_fixture_oracle(fixture: Mapping[str, object], manifest: Mapping[str, object]) -> tuple[dict[str, int], dict[str, bool]]:
    paths = _fixture_path_map(fixture)
    fixture_id = str(fixture["fixture_id"])
    if fixture_id == FIXTURE_IDS[0]:
        reference = _parse_toy_fasta(paths["combined.fa"])
        bait = _parse_toy_fasta(paths["bait.fa"])
        assertions = {"fai_nonempty": paths["combined.fa.fai"].stat().st_size > 0, "mmi_nonempty": paths["combined.mmi"].stat().st_size > 0, "bait_exact_bounds": bait == [(b"hapA:2-9", b"CGTACGTA")]}
        counts = {"haplotypes": len(reference), "bait_bases": sum(len(seq) for _, seq in bait)}
    elif fixture_id == FIXTURE_IDS[1]:
        source = _parse_toy_sam(paths["flags.sam"]); view = _parse_toy_sam(paths["primary.sam"]); fasta = _parse_toy_fasta(paths["primary.fa"])
        expected = sorted(row[0] for row in source if int(row[1]) & 0x900 == 0)
        assertions = {"view_primary_inventory": sorted(row[0] for row in view) == expected, "fasta_primary_inventory": sorted(name for name, _ in fasta) == expected, "decoded_base_total": sum(len(seq) for _, seq in fasta) == 24}
        counts = {"input_records": len(source), "primary_molecules": len(expected), "primary_bases": sum(len(seq) for _, seq in fasta)}
    elif fixture_id == FIXTURE_IDS[2]:
        reference_names = {name for name, _ in _parse_toy_fasta(paths["combined.fa"])}; reads = _parse_toy_fasta(paths["reads.fa"]); records = _parse_toy_sam(paths["aligned.sam"])
        primaries = [row[0] for row in records if int(row[1]) & 0x900 == 0]
        assertions = {"no_secondary": all(int(row[1]) & 0x100 == 0 for row in records), "one_primary_per_molecule": Counter(primaries) == Counter({name: 1 for name, _ in reads}), "combined_reference_only": all(row[2] == b"*" or row[2] in reference_names for row in records)}
        counts = {"decoded_molecules": len(reads), "primary_records": len(primaries)}
    elif fixture_id == FIXTURE_IDS[3]:
        source = _parse_toy_sam(paths["window.sam"]); selected = _parse_toy_sam(paths["selected.sam"])
        expected = sorted(row[0] for row in source if int(row[1]) & 0x900 == 0)
        assertions = {"interval_primary_only": all(int(row[1]) & 0x900 == 0 for row in selected), "qname_inventory_exact": sorted(row[0] for row in selected) == expected, "quickcheck_passed": any(command["argv"][1:2] == ["quickcheck"] for command in fixture["commands"]), "no_downsampling": Counter(row[0] for row in selected) == Counter(expected)}
        counts = {"primary_records": len(selected)}
    elif fixture_id == FIXTURE_IDS[4]:
        part_a = paths["partA.qnames"].read_bytes().splitlines(); part_b = paths["partB.qnames"].read_bytes().splitlines(); duplicate = paths["duplicate.qnames"].read_bytes().splitlines(); merged = _parse_toy_sam(paths["merged.sam"])
        union = part_a + part_b
        assertions = {"disjoint_union_preserved": len(set(union)) == len(union) and Counter(row[0] for row in merged) == Counter(union), "duplicate_qname_rejected": bool(set(part_a) & set(duplicate)), "combined_reference_contig_owned": all(row[2] == b"hapA" for row in merged), "quickcheck_passed": any(command["argv"][1:2] == ["quickcheck"] for command in fixture["commands"])}
        counts = {"part_records": len(union), "merged_records": len(merged)}
    elif fixture_id == FIXTURE_IDS[5]:
        raw, q10 = _fixture_depth_values(paths["raw.depth"]), _fixture_depth_values(paths["mapq10.depth"])
        assertions = {"all_bed_positions_emitted": len(raw) == len(q10) == 20, "duplicate_qcfail_included": raw == [4] * 10 + [0] * 10, "uppercase_q10_exact": q10 == [3] * 10 + [0] * 10, "raw_not_changed_by_q10": all(left >= right for left, right in zip(raw, q10))}
        counts = {"positions": 20, "raw_depth_sum": sum(raw), "mapq10_depth_sum": sum(q10)}
    elif fixture_id == FIXTURE_IDS[6]:
        plotter = _load_manifest_pinned_plotter(manifest)
        target_records = _parse_toy_fasta(paths["full-window.fa"])
        query_records = _parse_toy_fasta(paths["type2_KCON.fa"])
        authority_payload = paths["selector_authority.json"].read_bytes()
        authority = plotter._validate_kcon_selector_authority(strict_json_bytes(authority_payload, "fixture KCON authority"))
        resident_paf_rows, resident_paf_sha = plotter._parse_paf(paths["resident_minimap.paf"])
        resident_selected = plotter._select_kcon_annotations_core(authority, resident_paf_rows, resident_paf_sha)
        resident_independent = _independent_kcon_oracle(authority, resident_paf_rows, resident_paf_sha)
        paf_rows, paf_sha = plotter._parse_paf(paths["kcon.paf"])
        selected = plotter._select_kcon_annotations_core(authority, paf_rows, paf_sha)
        independent = _independent_kcon_oracle(authority, paf_rows, paf_sha)
        before = strict_json_bytes(paths["boundary_before.json"].read_bytes(), "fixture KCON boundary before")
        after = strict_json_bytes(paths["boundary_after.json"].read_bytes(), "fixture KCON boundary after")
        boundary_observed = {name: {"size_bytes": paths[name].stat().st_size, "sha256": sha256_file(paths[name])} for name in ("boundary.bam", "boundary.depth", "boundary.callability")}
        fields = [row[0] for row in paf_rows]
        exact_seams = (
            any(int(row[10]) == 499 for row in fields)
            and any(int(row[10]) == 500 for row in fields)
            and any(int(row[7]) >= 2500 for row in fields if row[5] == "toy:1-12000")
            and any(row[0] != "type2_KCON" for row in fields)
            and any(int(row[1]) != 9472 for row in fields)
            and any(row[5] != "toy:1-12000" for row in fields)
        )
        selected_numeric = (selected[0]["aligned_block_length"], selected[0]["matches"], selected[0]["mapq"], selected[0]["target_start_0based"], selected[0]["query_start_0based"])
        raw_tie_rows = [raw for row, raw in paf_rows if (int(row[10]), int(row[9]), int(row[11]), int(row[7]), int(row[2])) == selected_numeric]
        ranking_seams = (
            any((int(row[10]), int(row[9]), int(row[11])) == (700, 680, 50) for row in fields)
            and any((int(row[10]), int(row[9]), int(row[11])) == (700, 690, 40) for row in fields)
            and sum((int(row[10]), int(row[9]), int(row[11])) == (700, 690, 50) for row in fields) >= 2
            and selected[0]["target_start_0based"] == 1600
            and selected[0]["query_start_0based"] == 10
            and len(raw_tie_rows) >= 2
            and selected[0]["paf_row_sha256"] == hashlib.sha256(min(raw_tie_rows)).hexdigest()
        )
        exact_sequences = target_records == [(b"toy:1-12000", target_records[0][1])] and len(target_records[0][1]) == 12000 and query_records == [(b"type2_KCON", query_records[0][1])] and len(query_records[0][1]) == 9472
        resident_authority_bound = bool(resident_paf_rows) and resident_selected == resident_independent and bool(resident_selected) and all(row["query_name"] == "type2_KCON" and row["query_length"] == 9472 and row["target_name"] == "toy:1-12000" and row["target_length"] == 12000 for row in resident_selected)
        assertions = {"paf_bounds_valid": exact_sequences and exact_seams and resident_authority_bound, "deterministic_ranking_nonempty": selected == independent and bool(selected) and ranking_seams and resident_selected == resident_independent, "annotation_only": all(row["annotation_only"] is True and row["read_inferred_copy_number"] == "not_estimated" for row in selected + resident_selected), "does_not_change_bam_depth_callability": before == after == boundary_observed}
        counts = {"paf_rows": len(paf_rows), "selected_annotations": len(selected), "resident_paf_rows": len(resident_paf_rows), "resident_selected_annotations": len(resident_selected)}
    else:
        raise ExactFourPlanError("unknown fixture oracle")
    if not all(assertions.values()):
        raise ExactFourPlanError(f"retained scientific oracle failed: {fixture_id}")
    return counts, assertions


def validate_resident_tool_fixture_receipt(value: object, tools: object) -> dict[str, object]:
    """Validate the retained receipt from exactly seven offline real-binary fixtures."""
    accepted_tools = _validate_tool_receipts(tools)
    receipt = _require_dict(value, "resident tool fixture receipt")
    _exact_keys(receipt, ("schema_version", "host_class", "deployment_manifest", "authorization_binding", "tool_bindings", "fixture_order", "fixtures", "action_gates"), "resident tool fixture receipt")
    if receipt["schema_version"] != TOOL_FIXTURE_SCHEMA or receipt["host_class"] != "login.pax.tufts.edu":
        raise ExactFourPlanError("resident tool fixture schema/host drift")
    manifest, _ = _fixture_manifest_and_authorization(receipt, accepted_tools)
    bindings = _require_dict(receipt["tool_bindings"], "fixture tool bindings")
    _exact_keys(bindings, ("minimap2", "samtools"), "fixture tool bindings")
    by_name = {row["name"]: row for row in accepted_tools}
    for name in ("minimap2", "samtools"):
        binding = _require_dict(bindings[name], f"fixture {name} binding")
        _exact_keys(binding, ("path", "size_bytes", "sha256", "version", "help"), f"fixture {name} binding")
        if binding != {key: by_name[name][key] for key in ("path", "size_bytes", "sha256", "version", "help")}:
            raise ExactFourPlanError(f"fixture {name} binding differs from resident tool receipt")
    if receipt["fixture_order"] != list(FIXTURE_IDS) or type(receipt["fixtures"]) is not list or len(receipt["fixtures"]) != 7:
        raise ExactFourPlanError("resident fixture set/order is not exact seven")
    for index, (fixture_value, fixture_id) in enumerate(zip(receipt["fixtures"], FIXTURE_IDS)):
        fixture = _require_dict(fixture_value, f"resident fixture {index}")
        _exact_keys(fixture, ("fixture_id", "scratch_directory", "input_files", "commands", "output_files", "record_counts", "semantic_assertions", "passed"), f"resident fixture {index}")
        if fixture["fixture_id"] != fixture_id:
            raise ExactFourPlanError("resident fixture identity/order drift")
        scratch = Path(_require_str(fixture["scratch_directory"], "fixture scratch directory"))
        if not scratch.is_absolute() or not scratch.is_dir() or scratch.is_symlink():
            raise ExactFourPlanError("fixture scratch directory must be absolute")
        for collection_name in ("input_files", "output_files"):
            rows = fixture[collection_name]
            if type(rows) is not list:
                raise ExactFourPlanError(f"fixture {collection_name} must be a list")
            for row_value in rows:
                row = _require_dict(row_value, f"fixture {collection_name} row")
                _exact_keys(row, ("path", "size_bytes", "sha256"), f"fixture {collection_name} row")
                _require_int(row["size_bytes"], f"fixture {collection_name} size", minimum=0)
                _require_sha(row["sha256"], f"fixture {collection_name} SHA-256")
                artifact = Path(_require_str(row["path"], f"fixture {collection_name} path"))
                try:
                    if os.path.commonpath((str(scratch), str(artifact))) != str(scratch) or artifact == scratch:
                        raise ExactFourPlanError("fixture artifact escaped its private scratch directory")
                except ValueError as error:
                    raise ExactFourPlanError("fixture artifact path domain drift") from error
                observed = _fingerprint(artifact, allow_empty=True)
                if observed["size_bytes"] != row["size_bytes"] or observed["sha256"] != row["sha256"]:
                    raise ExactFourPlanError("retained fixture artifact bytes drifted on replay")
        commands = fixture["commands"]
        if type(commands) is not list or not commands:
            raise ExactFourPlanError("fixture must retain at least one exact command receipt")
        for command_value in commands:
            command = _require_dict(command_value, "fixture command")
            _exact_keys(command, ("argv", "returncode", "stdout_sha256", "stderr_sha256"), "fixture command")
            if type(command["argv"]) is not list or not command["argv"] or any(type(token) is not str for token in command["argv"]):
                raise ExactFourPlanError("fixture argv must be a nonempty exact string list")
            if command["argv"][0] not in {by_name["minimap2"]["path"], by_name["samtools"]["path"]} or command["returncode"] != 0 or type(command["returncode"]) is not int:
                raise ExactFourPlanError("fixture command escaped accepted binaries or failed")
            _require_sha(command["stdout_sha256"], "fixture stdout SHA-256")
            _require_sha(command["stderr_sha256"], "fixture stderr SHA-256")
            forbidden = " ".join(command["argv"]).lower()
            if any(token in forbidden for token in ("aws", "s3://", "sbatch", "squeue", "scancel", "ssh", "curl", "wget")):
                raise ExactFourPlanError("fixture command exposes provider/scheduler/network action")
        argvs = [command["argv"] for command in commands]
        minimap_path, samtools_path = by_name["minimap2"]["path"], by_name["samtools"]["path"]
        required_patterns: dict[str, tuple[tuple[str, ...], ...]] = {
            FIXTURE_IDS[0]: ((samtools_path, "faidx"), (minimap_path, "-d")),
            FIXTURE_IDS[1]: ((samtools_path, "view", "-b", "-o"), (samtools_path, "view", "-F", "0x900"), (samtools_path, "fasta", "-F", "0x900")),
            FIXTURE_IDS[2]: ((minimap_path, "-ax", "map-ont", "--secondary=no"),),
            FIXTURE_IDS[3]: ((samtools_path, "view", "-u", "-F", "0x900", "-L"), (samtools_path, "sort", "-@", "4"), (samtools_path, "index", "-@", "4"), (samtools_path, "quickcheck")),
            FIXTURE_IDS[4]: ((samtools_path, "merge", "-@", "4"), (samtools_path, "index", "-@", "4"), (samtools_path, "quickcheck")),
            FIXTURE_IDS[5]: ((samtools_path, "depth", "-a", "-g", "0x600", "-b"), (samtools_path, "depth", "-a", "-g", "0x600", "-Q", "10", "-b")),
            FIXTURE_IDS[6]: ((minimap_path, "-c", "-p", "0.1", "-N", "10"),),
        }
        for prefix in required_patterns[fixture_id]:
            if not any(tuple(argv[:len(prefix)]) == prefix for argv in argvs):
                raise ExactFourPlanError(f"fixture {fixture_id} lacks exact workflow argv {prefix}")
        counts = _require_dict(fixture["record_counts"], "fixture record counts")
        if any(type(key) is not str or type(count) is not int or count < 0 for key, count in counts.items()):
            raise ExactFourPlanError("fixture record counts must be exact nonnegative integers")
        assertions = _require_dict(fixture["semantic_assertions"], "fixture semantic assertions")
        required_assertions = {
            FIXTURE_IDS[0]: {"fai_nonempty", "mmi_nonempty", "bait_exact_bounds"},
            FIXTURE_IDS[1]: {"view_primary_inventory", "fasta_primary_inventory", "decoded_base_total"},
            FIXTURE_IDS[2]: {"no_secondary", "one_primary_per_molecule", "combined_reference_only"},
            FIXTURE_IDS[3]: {"interval_primary_only", "qname_inventory_exact", "quickcheck_passed", "no_downsampling"},
            FIXTURE_IDS[4]: {"disjoint_union_preserved", "duplicate_qname_rejected", "combined_reference_contig_owned", "quickcheck_passed"},
            FIXTURE_IDS[5]: {"all_bed_positions_emitted", "duplicate_qcfail_included", "uppercase_q10_exact", "raw_not_changed_by_q10"},
            FIXTURE_IDS[6]: {"paf_bounds_valid", "deterministic_ranking_nonempty", "annotation_only", "does_not_change_bam_depth_callability"},
        }[fixture_id]
        if set(assertions) != required_assertions or any(type(key) is not str or type(passed) is not bool or passed is not True for key, passed in assertions.items()):
            raise ExactFourPlanError("every fixture semantic assertion must be exact true")
        recomputed_counts, recomputed_assertions = _recompute_fixture_oracle(fixture, manifest)
        if counts != recomputed_counts or assertions != recomputed_assertions:
            raise ExactFourPlanError("fixture recorded science differs from retained-artifact replay")
        _require_bool(fixture["passed"], True, "fixture passed")
    gates = _require_dict(receipt["action_gates"], "fixture action gates")
    _exact_keys(gates, tuple(ACTION_GATES), "fixture action gates")
    if gates != ACTION_GATES or any(type(gates[key]) is not bool or gates[key] is not False for key in ACTION_GATES):
        raise ExactFourPlanError("fixture action gates must remain exactly false")
    return receipt


def run_resident_tool_fixtures(
    tools: object,
    manifest: object,
    authorization_path: object,
    authorization_size_bytes: int,
    authorization_sha256: str,
    acceptance_path: object,
    acceptance_size_bytes: int,
    acceptance_sha256: str,
) -> dict[str, object]:
    """Run the seven separately authorized offline toy fixtures once.

    This interface has no provider/scheduler executable and retains every
    artifact. It is intentionally not called by resident capture.
    """
    accepted_tools = _validate_tool_receipts(tools)
    by_name = {row["name"]: row for row in accepted_tools}
    accepted_manifest = _validate_deployment_manifest(manifest)
    auth, authorization_binding = _validate_fixture_authorization_files(
        authorization_path, authorization_size_bytes, authorization_sha256,
        acceptance_path, acceptance_size_bytes, acceptance_sha256,
        accepted_tools, accepted_manifest,
    )
    root = _require_lexical_normalized_absolute(auth["scratch_root"], "fixture scratch root")
    if root.exists() or not root.parent.is_dir() or root.parent.is_symlink():
        raise ExactFourPlanError("fixture scratch root must be an absent authorized absolute child")

    minimap = by_name["minimap2"]
    samtools = by_name["samtools"]
    root.mkdir(mode=0o700)

    def write(path: Path, payload: bytes) -> None:
        if not path.parent.is_dir() or path.exists():
            raise ExactFourPlanError("fixture input/output path is not create-only")
        path.write_bytes(payload)

    def fp(path: Path) -> dict[str, object]:
        return _fingerprint(path, allow_empty=True)

    def command(tool: Mapping[str, object], argv: list[str], *, stdout_path: Path | None = None, stdin_path: Path | None = None) -> tuple[dict[str, object], bytes, bytes]:
        executable = Path(str(tool["path"]))
        before = _fingerprint(executable)
        if before["size_bytes"] != tool["size_bytes"] or before["sha256"] != tool["sha256"]:
            raise ExactFourPlanError("resident fixture executable drift before subprocess")
        stdin_handle = stdin_path.open("rb") if stdin_path is not None else subprocess.DEVNULL
        try:
            completed = subprocess.run([str(executable), *argv], stdin=stdin_handle, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False, env={"LANG": "C", "LC_ALL": "C", "PATH": "/dev/null"})
        finally:
            if stdin_path is not None:
                stdin_handle.close()
        after = _fingerprint(executable)
        if after != before:
            raise ExactFourPlanError("resident fixture executable drift after subprocess")
        if completed.returncode != 0:
            raise ExactFourPlanError(f"resident fixture command failed: {[str(executable), *argv]}")
        if stdout_path is not None:
            write(stdout_path, completed.stdout)
        row = {"argv": [str(executable), *argv], "returncode": completed.returncode, "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(), "stderr_sha256": hashlib.sha256(completed.stderr).hexdigest()}
        return row, completed.stdout, completed.stderr

    def start_fixture(fixture_id: str) -> tuple[Path, list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
        directory = root / fixture_id
        directory.mkdir(mode=0o700)
        return directory, [], [], []

    fixtures: list[dict[str, object]] = []
    reference = b">hapA\nACGTACGTACGTACGTACGT\n>hapB\nTTGCTTGCTTGCTTGCTTGC\n"

    # 1. Combined reference, exact regions, FAI and MMI.
    directory, inputs, commands, outputs = start_fixture(FIXTURE_IDS[0])
    ref = directory / "combined.fa"; write(ref, reference); inputs.append(fp(ref))
    row, _, _ = command(samtools, ["faidx", str(ref)]); commands.append(row); outputs.append(fp(Path(str(ref) + ".fai")))
    bait = directory / "bait.fa"; row, bait_bytes, _ = command(samtools, ["faidx", str(ref), "hapA:2-9"], stdout_path=bait); commands.append(row); outputs.append(fp(bait))
    mmi = directory / "combined.mmi"; row, _, _ = command(minimap, ["-d", str(mmi), str(ref)]); commands.append(row); outputs.append(fp(mmi))
    assertions = {"fai_nonempty": outputs[0]["size_bytes"] > 0, "mmi_nonempty": outputs[2]["size_bytes"] > 0, "bait_exact_bounds": b"CGTACGTA" in bait_bytes}
    fixtures.append({"fixture_id": FIXTURE_IDS[0], "scratch_directory": str(directory), "input_files": inputs, "commands": commands, "output_files": outputs, "record_counts": {"haplotypes": 2, "bait_bases": 8}, "semantic_assertions": assertions, "passed": all(assertions.values())})

    # 2. Primary-source decoding: -F 0x900 excludes only secondary/supplementary.
    directory, inputs, commands, outputs = start_fixture(FIXTURE_IDS[1])
    sam = directory / "flags.sam"
    sam_rows = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:hapA\tLN:20"]
    for name, flag in (("primary", 0), ("secondary", 256), ("supplementary", 2048), ("duplicate", 1024), ("qcfail", 512)):
        sam_rows.append(f"{name}\t{flag}\thapA\t1\t20\t8M\t*\t0\t0\tACGTACGT\tFFFFFFFF")
    write(sam, ("\n".join(sam_rows) + "\n").encode("ascii")); inputs.append(fp(sam))
    bam = directory / "flags.bam"; row, _, _ = command(samtools, ["view", "-b", "-o", str(bam), str(sam)]); commands.append(row); outputs.append(fp(bam))
    view = directory / "primary.sam"; row, view_bytes, _ = command(samtools, ["view", "-F", "0x900", str(bam)], stdout_path=view); commands.append(row); outputs.append(fp(view))
    fasta = directory / "primary.fa"; row, fasta_bytes, _ = command(samtools, ["fasta", "-F", "0x900", str(bam)], stdout_path=fasta); commands.append(row); outputs.append(fp(fasta))
    view_names = [line.split(b"\t", 1)[0] for line in view_bytes.splitlines()]
    fasta_names = [line[1:] for line in fasta_bytes.splitlines() if line.startswith(b">")]
    expected_names = [b"primary", b"duplicate", b"qcfail"]
    assertions = {"view_primary_inventory": sorted(view_names) == sorted(expected_names), "fasta_primary_inventory": sorted(fasta_names) == sorted(expected_names), "decoded_base_total": sum(len(line) for line in fasta_bytes.splitlines() if not line.startswith(b">")) == 24}
    fixtures.append({"fixture_id": FIXTURE_IDS[1], "scratch_directory": str(directory), "input_files": inputs, "commands": commands, "output_files": outputs, "record_counts": {"input_records": 5, "primary_molecules": 3, "primary_bases": 24}, "semantic_assertions": assertions, "passed": all(assertions.values())})

    # 3. One competitive reference, one primary per decoded molecule, no secondary.
    directory, inputs, commands, outputs = start_fixture(FIXTURE_IDS[2])
    ref = directory / "combined.fa"; write(ref, reference); inputs.append(fp(ref))
    reads = directory / "reads.fa"; write(reads, b">readA\nACGTACGTACGT\n>readB\nTTGCTTGCTTGC\n"); inputs.append(fp(reads))
    aligned = directory / "aligned.sam"; row, aligned_bytes, _ = command(minimap, ["-ax", "map-ont", "--secondary=no", str(ref), str(reads)], stdout_path=aligned); commands.append(row); outputs.append(fp(aligned))
    records = [line.split(b"\t") for line in aligned_bytes.splitlines() if line and not line.startswith(b"@")]
    primary_names = [fields[0] for fields in records if int(fields[1]) & 0x900 == 0]
    assertions = {"no_secondary": all(int(fields[1]) & 0x100 == 0 for fields in records), "one_primary_per_molecule": sorted(primary_names) == [b"readA", b"readB"], "combined_reference_only": len(inputs) == 2}
    fixtures.append({"fixture_id": FIXTURE_IDS[2], "scratch_directory": str(directory), "input_files": inputs, "commands": commands, "output_files": outputs, "record_counts": {"decoded_molecules": 2, "primary_records": len(primary_names)}, "semantic_assertions": assertions, "passed": all(assertions.values())})

    # 4. Exact primary full-window view, sort, index and quickcheck.
    directory, inputs, commands, outputs = start_fixture(FIXTURE_IDS[3])
    local_sam = directory / "window.sam"
    local_rows = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:hapA\tLN:20"]
    for name, flag in (("primary", 0), ("secondary", 256), ("supplementary", 2048), ("duplicate", 1024), ("qcfail", 512)):
        local_rows.append(f"{name}\t{flag}\thapA\t1\t20\t8M\t*\t0\t0\tACGTACGT\tFFFFFFFF")
    write(local_sam, ("\n".join(local_rows) + "\n").encode("ascii")); inputs.append(fp(local_sam))
    local_bam = directory / "window.bam"
    row, _, _ = command(samtools, ["view", "-b", "-o", str(local_bam), str(local_sam)]); commands.append(row); inputs.append(fp(local_bam))
    bed = directory / "window.bed"; write(bed, b"hapA\t0\t12\n"); inputs.append(fp(bed))
    uncompressed = directory / "primary.unsorted.bam"
    row, _, _ = command(samtools, ["view", "-u", "-F", "0x900", "-L", str(bed), str(local_bam)], stdout_path=uncompressed); commands.append(row); outputs.append(fp(uncompressed))
    sorted_bam = directory / "primary.bam"
    row, _, _ = command(samtools, ["sort", "-@", "4", "-o", str(sorted_bam), str(uncompressed)]); commands.append(row); outputs.append(fp(sorted_bam))
    row, _, _ = command(samtools, ["index", "-@", "4", str(sorted_bam)]); commands.append(row); outputs.append(fp(Path(str(sorted_bam) + ".bai")))
    row, _, _ = command(samtools, ["quickcheck", str(sorted_bam)]); commands.append(row)
    selected_path = directory / "selected.sam"
    row, selected, _ = command(samtools, ["view", str(sorted_bam)], stdout_path=selected_path); commands.append(row); outputs.append(fp(selected_path))
    selected_fields = [line.split(b"\t") for line in selected.splitlines()]
    assertions = {"interval_primary_only": all(int(fields[1]) & 0x900 == 0 for fields in selected_fields), "qname_inventory_exact": sorted(fields[0] for fields in selected_fields) == [b"duplicate", b"primary", b"qcfail"], "quickcheck_passed": True, "no_downsampling": True}
    fixtures.append({"fixture_id": FIXTURE_IDS[3], "scratch_directory": str(directory), "input_files": inputs, "commands": commands, "output_files": outputs, "record_counts": {"primary_records": len(selected_fields)}, "semantic_assertions": assertions, "passed": all(assertions.values())})

    # 5. Exact part merge; Python ownership oracle rejects duplicate QNAMEs.
    directory, inputs, commands, outputs = start_fixture(FIXTURE_IDS[4])
    part_bams = []
    part_names = (("partA", "a1"), ("partB", "b1"))
    for label, qname in part_names:
        part_sam = directory / f"{label}.sam"
        write(part_sam, ("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:hapA\tLN:20\n" + f"{qname}\t0\thapA\t1\t20\t8M\t*\t0\t0\tACGTACGT\tFFFFFFFF\n").encode("ascii")); inputs.append(fp(part_sam))
        part_bam = directory / f"{label}.bam"
        row, _, _ = command(samtools, ["view", "-b", "-o", str(part_bam), str(part_sam)]); commands.append(row); outputs.append(fp(part_bam)); part_bams.append(part_bam)
    merged = directory / "merged.bam"
    row, _, _ = command(samtools, ["merge", "-@", "4", "-o", str(merged), *map(str, part_bams)]); commands.append(row); outputs.append(fp(merged))
    row, _, _ = command(samtools, ["index", "-@", "4", str(merged)]); commands.append(row); outputs.append(fp(Path(str(merged) + ".bai")))
    row, _, _ = command(samtools, ["quickcheck", str(merged)]); commands.append(row)
    merged_sam_path = directory / "merged.sam"
    row, merged_sam, _ = command(samtools, ["view", str(merged)], stdout_path=merged_sam_path); commands.append(row); outputs.append(fp(merged_sam_path))
    merged_fields = [line.split(b"\t") for line in merged_sam.splitlines()]
    inventory_a = directory / "partA.qnames"; write(inventory_a, b"a1\n"); inputs.append(fp(inventory_a))
    inventory_b = directory / "partB.qnames"; write(inventory_b, b"b1\n"); inputs.append(fp(inventory_b))
    union = verify_qname_inventory_union([inventory_a, inventory_b])
    duplicated = directory / "duplicate.qnames"; write(duplicated, b"a1\n"); inputs.append(fp(duplicated))
    duplicate_rejected = False
    try:
        verify_qname_inventory_union([inventory_a, duplicated])
    except ExactFourPlanError:
        duplicate_rejected = True
    assertions = {"disjoint_union_preserved": sorted(fields[0] for fields in merged_fields) == [b"a1", b"b1"] and union["count"] == 2, "duplicate_qname_rejected": duplicate_rejected, "combined_reference_contig_owned": all(fields[2] == b"hapA" for fields in merged_fields), "quickcheck_passed": True}
    fixtures.append({"fixture_id": FIXTURE_IDS[4], "scratch_directory": str(directory), "input_files": inputs, "commands": commands, "output_files": outputs, "record_counts": {"part_records": 2, "merged_records": len(merged_fields)}, "semantic_assertions": assertions, "passed": all(assertions.values())})

    # 6. Raw and uppercase-Q10 depth with explicit -g 0x600 inclusion.
    directory, inputs, commands, outputs = start_fixture(FIXTURE_IDS[5])
    depth_sam = directory / "depth.sam"
    depth_rows = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:hapA\tLN:20"]
    for name, flag, mapq in (("low", 0, 5), ("high", 0, 20), ("duplicate", 1024, 20), ("qcfail", 512, 20)):
        depth_rows.append(f"{name}\t{flag}\thapA\t1\t{mapq}\t10M\t*\t0\t0\tACGTACGTAC\tFFFFFFFFFF")
    write(depth_sam, ("\n".join(depth_rows) + "\n").encode("ascii")); inputs.append(fp(depth_sam))
    depth_unsorted = directory / "depth.unsorted.bam"; row, _, _ = command(samtools, ["view", "-b", "-o", str(depth_unsorted), str(depth_sam)]); commands.append(row)
    depth_bam = directory / "depth.bam"; row, _, _ = command(samtools, ["sort", "-@", "4", "-o", str(depth_bam), str(depth_unsorted)]); commands.append(row); outputs.append(fp(depth_bam))
    row, _, _ = command(samtools, ["index", "-@", "4", str(depth_bam)]); commands.append(row); outputs.append(fp(Path(str(depth_bam) + ".bai")))
    depth_bed = directory / "depth.bed"; write(depth_bed, b"hapA\t0\t20\n"); inputs.append(fp(depth_bed))
    raw_path = directory / "raw.depth"; row, raw_bytes, _ = command(samtools, ["depth", "-a", "-g", "0x600", "-b", str(depth_bed), str(depth_bam)], stdout_path=raw_path); commands.append(row); outputs.append(fp(raw_path))
    q10_path = directory / "mapq10.depth"; row, q10_bytes, _ = command(samtools, ["depth", "-a", "-g", "0x600", "-Q", "10", "-b", str(depth_bed), str(depth_bam)], stdout_path=q10_path); commands.append(row); outputs.append(fp(q10_path))
    def depth_values(payload: bytes) -> list[int]:
        rows = [line.split(b"\t") for line in payload.splitlines()]
        if [(fields[0], int(fields[1])) for fields in rows] != [(b"hapA", index) for index in range(1, 21)]:
            return []
        return [int(fields[2]) for fields in rows]
    raw_values, q10_values = depth_values(raw_bytes), depth_values(q10_bytes)
    assertions = {"all_bed_positions_emitted": len(raw_values) == len(q10_values) == 20, "duplicate_qcfail_included": raw_values == [4] * 10 + [0] * 10, "uppercase_q10_exact": q10_values == [3] * 10 + [0] * 10, "raw_not_changed_by_q10": all(raw >= q10 for raw, q10 in zip(raw_values, q10_values))}
    fixtures.append({"fixture_id": FIXTURE_IDS[5], "scratch_directory": str(directory), "input_files": inputs, "commands": commands, "output_files": outputs, "record_counts": {"positions": 20, "raw_depth_sum": sum(raw_values), "mapq10_depth_sum": sum(q10_values)}, "semantic_assertions": assertions, "passed": all(assertions.values())})

    # 7. Production-core annotation-only PAF selection under exact seams.
    directory, inputs, commands, outputs = start_fixture(FIXTURE_IDS[6])
    target = directory / "full-window.fa"; query = directory / "type2_KCON.fa"
    target_sequence = bytes(b"ACGT"[((index * 1103515245 + 12345) >> 16) & 3] for index in range(12000))
    query_sequence = target_sequence[500:9972]
    write(target, b">toy:1-12000\n" + target_sequence + b"\n"); write(query, b">type2_KCON\n" + query_sequence + b"\n")
    inputs.extend((fp(target), fp(query)))
    authority_path = directory / "selector_authority.json"
    authority_value = {
        "sample_id": "FIXTURE", "haplotypes": ["hapA", "hapB"],
        "bed_authorities": [{"haplotype": "hapA", "contig": "toy", "full_intervals": [[0, 12000]]}],
        "element_identities": [{"haplotype": "hapA", "contig": "toy", "source_start": 1500, "source_end": 2500}],
    }
    write(authority_path, canonical_json_bytes(authority_value)); inputs.append(fp(authority_path))
    boundary_files = {}
    for name, payload in (("boundary.bam", b"synthetic-bam-boundary\n"), ("boundary.depth", b"synthetic-depth-boundary\n"), ("boundary.callability", b"synthetic-callability-boundary\n")):
        path = directory / name; write(path, payload); inputs.append(fp(path)); boundary_files[name] = fp(path)
    boundary_value = {name: {"size_bytes": row["size_bytes"], "sha256": row["sha256"]} for name, row in boundary_files.items()}
    before_path = directory / "boundary_before.json"; write(before_path, canonical_json_bytes(boundary_value)); inputs.append(fp(before_path))
    resident_paf = directory / "resident_minimap.paf"
    row, _, _ = command(minimap, ["-c", "-p", "0.1", "-N", "10", str(target), str(query)], stdout_path=resident_paf)
    commands.append(row); outputs.append(fp(resident_paf))
    # The selector seam is a deterministic synthetic PAF over the same exact
    # authority/query.  It retains 499/500, overlap, and every ranking tier;
    # the real resident invocation above separately proves the exact argv.
    paf = directory / "kcon.paf"
    write(paf, b"".join((
        b"type2_KCON\t9472\t0\t499\t+\ttoy:1-12000\t12000\t1500\t1999\t499\t499\t60\n",
        b"type2_KCON\t9472\t0\t500\t+\ttoy:1-12000\t12000\t1500\t2000\t500\t500\t60\n",
        b"type2_KCON\t9472\t50\t750\t+\ttoy:1-12000\t12000\t1700\t2400\t680\t700\t50\n",
        b"type2_KCON\t9472\t40\t740\t+\ttoy:1-12000\t12000\t1700\t2400\t690\t700\t40\n",
        b"type2_KCON\t9472\t30\t730\t+\ttoy:1-12000\t12000\t1700\t2400\t690\t700\t50\n",
        b"type2_KCON\t9472\t20\t720\t+\ttoy:1-12000\t12000\t1600\t2300\t690\t700\t50\n",
        b"type2_KCON\t9472\t10\t710\t+\ttoy:1-12000\t12000\t1600\t2300\t690\t700\t50\tzz:Z:b\n",
        b"type2_KCON\t9472\t10\t710\t+\ttoy:1-12000\t12000\t1600\t2300\t690\t700\t50\tzz:Z:a\n",
        b"type2_KCON\t9472\t0\t700\t+\ttoy:1-12000\t12000\t2600\t3300\t700\t700\t60\n",
        b"wrong_query\t9472\t0\t700\t+\ttoy:1-12000\t12000\t1600\t2300\t700\t700\t60\n",
        b"type2_KCON\t9471\t0\t700\t+\ttoy:1-12000\t12000\t1600\t2300\t700\t700\t60\n",
        b"type2_KCON\t9472\t0\t700\t+\twrong:1-12000\t12000\t1600\t2300\t700\t700\t60\n",
    )))
    outputs.append(fp(paf))
    after_value = {name: {"size_bytes": fp(directory / name)["size_bytes"], "sha256": fp(directory / name)["sha256"]} for name in boundary_files}
    after_path = directory / "boundary_after.json"; write(after_path, canonical_json_bytes(after_value)); outputs.append(fp(after_path))
    provisional = {"fixture_id": FIXTURE_IDS[6], "scratch_directory": str(directory), "input_files": inputs, "commands": commands, "output_files": outputs, "record_counts": {}, "semantic_assertions": {}, "passed": True}
    counts, assertions = _recompute_fixture_oracle(provisional, accepted_manifest)
    provisional["record_counts"] = counts; provisional["semantic_assertions"] = assertions; provisional["passed"] = all(assertions.values()); fixtures.append(provisional)

    receipt = {
        "schema_version": TOOL_FIXTURE_SCHEMA,
        "host_class": "login.pax.tufts.edu",
        "deployment_manifest": {"path": str(_project_path_from_relative(DEPLOYMENT_MANIFEST_PATH)), "size_bytes": len(canonical_json_bytes(accepted_manifest)), "sha256": digest_value(accepted_manifest)},
        "authorization_binding": authorization_binding,
        "tool_bindings": {name: {key: by_name[name][key] for key in ("path", "size_bytes", "sha256", "version", "help")} for name in ("minimap2", "samtools")},
        "fixture_order": list(FIXTURE_IDS),
        "fixtures": fixtures,
        "action_gates": dict(ACTION_GATES),
    }
    if not all(row["passed"] is True for row in fixtures):
        raise ExactFourPlanError("one or more real-binary fixture oracles failed")
    return validate_resident_tool_fixture_receipt(receipt, tools)


def _bed_bytes(beds: Sequence[Mapping[str, object]], sample: str, interval_key: str, suffix: str) -> bytes:
    rows: list[tuple[str, int, int, str]] = []
    for bed in beds:
        if bed["sample_id"] != sample:
            continue
        for start, end in bed[interval_key]:
            rows.append((str(bed["contig"]), int(start), int(end), str(bed["haplotype"]) + suffix))
    rows.sort(key=lambda row: (row[0].encode("utf-8"), row[1], row[2], row[3].encode("utf-8")))
    return b"".join(("\t".join(map(str, row)) + "\n").encode("ascii") for row in rows)


def _derive_resident_biology(
    manifest: Mapping[str, object],
    resident_rows: Sequence[Mapping[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], dict[str, object], list[dict[str, object]]]:
    """Freshly derive all coordinate-bearing receipt material from pinned bytes."""
    authority_path = Path(next(row["path"] for row in resident_rows if row["role"] == "assembly_authority"))
    header, authority = _strict_tsv(authority_path)
    if header != ["sample_id", "haplotype", "fasta_path", "fai_path", "fasta_sha256", "fai_sha256"]:
        raise ExactFourPlanError("assembly authority header drift")
    expected_assemblies = manifest["assembly_authorities"]
    assemblies: list[dict[str, object]] = []
    for expected in expected_assemblies:
        matches = [row for row in authority if row["sample_id"] == expected["sample_id"] and row["haplotype"] == expected["haplotype"]]
        if len(matches) != 1:
            raise ExactFourPlanError("assembly authority must contain each exact sample/haplotype once")
        source = matches[0]
        if any(source[key] != expected[key] for key in ("fasta_path", "fai_path", "fasta_sha256", "fai_sha256")):
            raise ExactFourPlanError("assembly authority identity differs from deployment manifest")
        fasta, fai = Path(source["fasta_path"]), Path(source["fai_path"])
        fasta_fp, fai_fp = _fingerprint(fasta), _fingerprint(fai)
        if fasta_fp["sha256"] != source["fasta_sha256"] or fai_fp["sha256"] != source["fai_sha256"]:
            raise ExactFourPlanError("FASTA/FAI bytes drifted from authority")
        assemblies.append({
            "sample_id": source["sample_id"], "haplotype": source["haplotype"],
            "fasta": fasta_fp, "fai": fai_fp,
            "contigs": _read_fai(fai),
        })
    for sample in TARGETS:
        sample_assemblies = [row for row in assemblies if row["sample_id"] == sample]
        if len(sample_assemblies) != 2:
            raise ExactFourPlanError(f"sample {sample} must have exactly two declared assemblies")
        duplicate_contigs = set(sample_assemblies[0]["contigs"]) & set(sample_assemblies[1]["contigs"])
        if duplicate_contigs:
            raise ExactFourPlanError(f"duplicate contig identifier across {sample} haplotypes")

    combined_path = Path(next(row["path"] for row in resident_rows if row["role"] == "combined_table"))
    combined_header, combined_rows = _strict_tsv(combined_path)
    needed = {"Locus", "ID", "Haplotype", "Source_Identifier"}
    if not needed.issubset(combined_header):
        raise ExactFourPlanError("combined table lacks Locus/ID/Haplotype/Source_Identifier")
    elements: list[dict[str, object]] = []
    combined_row_sets: list[dict[str, object]] = []
    seen: set[tuple[object, ...]] = set()
    assembly_lookup = {(row["sample_id"], row["haplotype"]): row for row in assemblies}
    for sample in TARGETS:
        selected = [row for row in combined_rows if row["Locus"] == LOCUS and row["ID"] == sample]
        combined_row_sets.append({"sample_id": sample, "locus": LOCUS, "header": combined_header, "rows": selected})
    for row_set in combined_row_sets:
        for row in row_set["rows"]:
            key = (row["ID"], row["Haplotype"])
            if key not in assembly_lookup:
                raise ExactFourPlanError(f"combined row has unauthenticated haplotype: {key}")
            match = SOURCE_RE.fullmatch(row["Source_Identifier"])
            if match is None:
                raise ExactFourPlanError("Source_Identifier must be one strict contig:start-end value")
            contig, start_text, end_text = match.groups()
            start, end = int(start_text), int(end_text)
            contigs = assembly_lookup[key]["contigs"]
            if contig not in contigs or not 0 <= start < end <= contigs[contig]:
                raise ExactFourPlanError("source interval absent from authenticated FAI or out of bounds")
            identity = (row["ID"], row["Haplotype"], contig, start, end)
            if identity in seen:
                raise ExactFourPlanError("duplicate assembly-encoded element identity")
            seen.add(identity)
            core_start, core_end = max(0, start - 1000), min(contigs[contig], end + 1000)
            full_start, full_end = max(0, core_start - 10000), min(contigs[contig], core_end + 10000)
            elements.append({
                "sample_id": row["ID"], "haplotype": row["Haplotype"], "contig": contig,
                "source_start": start, "source_end": end,
                "core_start": core_start, "core_end": core_end,
                "full_start": full_start, "full_end": full_end,
            })
    assembly_order = {(row["sample_id"], row["haplotype"]): index for index, row in enumerate(assemblies)}
    elements.sort(key=lambda row: (TARGETS.index(row["sample_id"]), assembly_order[(row["sample_id"], row["haplotype"])], row["contig"].encode("utf-8"), row["source_start"], row["source_end"]))
    if {row["sample_id"] for row in elements} != set(TARGETS):
        raise ExactFourPlanError("combined table must contain at least one exact locus row for every target")
    beds = []
    for sample in TARGETS:
        for assembly in [row for row in assemblies if row["sample_id"] == sample]:
            subset = [row for row in elements if row["sample_id"] == sample and row["haplotype"] == assembly["haplotype"]]
            by_contig: dict[str, list[dict[str, object]]] = {}
            for row in subset:
                by_contig.setdefault(row["contig"], []).append(row)
            for contig, members in sorted(by_contig.items()):
                core = _merge_intervals([(row["core_start"], row["core_end"]) for row in members])
                full = _merge_intervals([(row["full_start"], row["full_end"]) for row in members])
                flanks = _subtract_intervals(full, core)
                beds.append({"sample_id": sample, "haplotype": assembly["haplotype"], "contig": contig,
                             "core_intervals": [list(item) for item in core], "full_intervals": [list(item) for item in full],
                             "flank_intervals": [list(item) for item in flanks]})

    count_rows = [{"sample_id": row["sample_id"], "haplotype": row["haplotype"], "count": sum(1 for element in elements if element["sample_id"] == row["sample_id"] and element["haplotype"] == row["haplotype"])} for row in assemblies]
    sample_totals = [{"sample_id": sample, "count": sum(1 for element in elements if element["sample_id"] == sample)} for sample in TARGETS]
    if any(row["count"] < 1 for row in count_rows + sample_totals):
        raise ExactFourPlanError("each declared sample/haplotype must contain assembly-encoded element evidence")
    assembly_encoded_counts = {"by_assembly": count_rows, "by_sample": sample_totals}
    sample_beds = []
    for sample in TARGETS:
        core = _bed_bytes(beds, sample, "core_intervals", "_core")
        full = _bed_bytes(beds, sample, "full_intervals", "_full")
        sample_beds.append({
            "sample_id": sample,
            "core_bed": {"size_bytes": len(core), "sha256": hashlib.sha256(core).hexdigest(), "text": core.decode("ascii")},
            "full_window_bed": {"size_bytes": len(full), "sha256": hashlib.sha256(full).hexdigest(), "text": full.decode("ascii")},
        })
    return assemblies, elements, beds, combined_row_sets, assembly_encoded_counts, sample_beds


def capture_resident_receipt(
    deployment_manifest: dict[str, object],
    tool_paths: dict[str, str],
    tool_fixture_receipt: dict[str, object],
) -> dict[str, object]:
    manifest = _validate_deployment_manifest(deployment_manifest)
    if type(tool_paths) is not dict or set(tool_paths) != {"bash", "python", "aws", "minimap2", "samtools", "sbatch"}:
        raise ExactFourPlanError("tool_paths must contain exactly bash/python/aws/minimap2/samtools/sbatch")
    resident_rows = []
    for expected in manifest["resident_authorities"]:
        row = _require_dict(expected, "resident authority")
        path = Path(_require_str(row["path"], "resident path"))
        fingerprint = _fingerprint(path)
        if fingerprint["size_bytes"] != row["size_bytes"] or fingerprint["sha256"] != row["sha256"]:
            raise ExactFourPlanError(f"resident authority drift: {path}")
        resident_rows.append({"role": row["role"], **fingerprint})
    assemblies, elements, beds, combined_row_sets, assembly_encoded_counts, sample_beds = _derive_resident_biology(manifest, resident_rows)
    authority_path = Path(next(row["path"] for row in resident_rows if row["role"] == "assembly_authority"))

    tools = []
    version_argv = {"bash": ["--version"], "python": ["--version"], "aws": ["--version"], "minimap2": ["--version"], "samtools": ["--version"], "sbatch": ["--version"]}
    help_argv = {"bash": ["--help"], "python": ["--help"], "aws": ["--help"], "minimap2": ["--help"], "samtools": ["depth", "--help"], "sbatch": ["--help"]}
    for name in ("bash", "python", "aws", "minimap2", "samtools", "sbatch"):
        path = tool_paths[name]
        fp = _fingerprint(Path(path))
        tools.append({"name": name, **fp, "version": _tool_command(path, version_argv[name], f"{name} version"), "help": _tool_command(path, help_argv[name], f"{name} help")})
    _validate_tool_receipts(tools)
    fixtures = validate_resident_tool_fixture_receipt(tool_fixture_receipt, tools)
    statvfs = os.statvfs(authority_path)
    receipt: dict[str, object] = {
        "schema_version": RESIDENT_SCHEMA,
        "deployment_bundle_id": manifest["bundle_id"],
        "deployment_manifest_sha256": digest_value(manifest),
        "resident_files": resident_rows,
        "assemblies": assemblies,
        "element_identities": elements,
        "bed_authorities": beds,
        "combined_row_sets": combined_row_sets,
        "assembly_encoded_counts": assembly_encoded_counts,
        "sample_bed_authorities": sample_beds,
        "tools": tools,
        "tool_fixture_receipt": fixtures,
        "storage": {"path": str(authority_path.parent), "block_size": statvfs.f_frsize, "blocks_available": statvfs.f_bavail, "inodes_available": statvfs.f_favail},
        "resources": manifest["resources"],
        "capture_metadata": {"observation_label": "read_only_local_capture", "network_contact_performed": False, "scheduler_query_performed": False, "filesystem_write_performed": False},
    }
    return _validate_resident_receipt(receipt, manifest)


def _validate_resident_receipt(value: object, manifest: Mapping[str, object]) -> dict[str, object]:
    receipt = _require_dict(value, "resident receipt")
    _exact_keys(receipt, ("schema_version", "deployment_bundle_id", "deployment_manifest_sha256", "resident_files", "assemblies", "element_identities", "bed_authorities", "combined_row_sets", "assembly_encoded_counts", "sample_bed_authorities", "tools", "tool_fixture_receipt", "storage", "resources", "capture_metadata"), "resident receipt")
    if receipt["schema_version"] != RESIDENT_SCHEMA or receipt["deployment_bundle_id"] != manifest["bundle_id"]:
        raise ExactFourPlanError("resident receipt schema/bundle drift")
    if receipt["deployment_manifest_sha256"] != digest_value(manifest):
        raise ExactFourPlanError("resident receipt deployment manifest digest drift")
    tools = _validate_tool_receipts(receipt["tools"])
    if tools != receipt["tools"]:
        raise ExactFourPlanError("resident tool receipt normalization drift")
    validate_resident_tool_fixture_receipt(receipt["tool_fixture_receipt"], tools)
    resident_files = receipt["resident_files"]
    if type(resident_files) is not list or len(resident_files) != len(manifest["resident_authorities"]):
        raise ExactFourPlanError("resident receipt must contain the exact four resident files")
    for index, (row_value, expected) in enumerate(zip(resident_files, manifest["resident_authorities"])):
        row = _require_dict(row_value, f"resident file {index}")
        _exact_keys(row, ("role", "path", "size_bytes", "sha256"), f"resident file {index}")
        if row != expected:
            raise ExactFourPlanError(f"resident file identity/order drift at index {index}")
    derived = _derive_resident_biology(manifest, resident_files)
    for key, expected in zip(("assemblies", "element_identities", "bed_authorities", "combined_row_sets", "assembly_encoded_counts", "sample_bed_authorities"), derived):
        if receipt[key] != expected:
            raise ExactFourPlanError(f"resident coordinate/biology authority drift: {key}")
    if receipt["resources"] != manifest["resources"]:
        raise ExactFourPlanError("resident resource receipt drift")
    capture = _require_dict(receipt["capture_metadata"], "capture metadata")
    _exact_keys(capture, ("observation_label", "network_contact_performed", "scheduler_query_performed", "filesystem_write_performed"), "capture metadata")
    for key in ("network_contact_performed", "scheduler_query_performed", "filesystem_write_performed"):
        _require_bool(capture[key], False, key)
    _require_str(capture["observation_label"], "resident observation label")
    storage = _require_dict(receipt["storage"], "resident storage")
    _exact_keys(storage, ("path", "block_size", "blocks_available", "inodes_available"), "resident storage")
    if not Path(_require_str(storage["path"], "resident storage path")).is_absolute():
        raise ExactFourPlanError("resident storage path must be absolute")
    for key in ("block_size", "blocks_available", "inodes_available"):
        _require_int(storage[key], f"resident storage {key}", minimum=1)
    return receipt


def _validate_amended_live_capture_authorization(
    value: object,
    manifest: Mapping[str, object],
    *,
    revalidate_files: bool = True,
) -> dict[str, object]:
    binding = _require_dict(value, "amended live-capture authorization")
    keys = (
        "authorization_path", "authorization_size_bytes", "authorization_sha256",
        "acceptance_path", "acceptance_size_bytes", "acceptance_sha256",
        "deployment_bundle_id", "deployment_manifest_path",
        "deployment_manifest_size_bytes", "deployment_manifest_sha256",
        "provider_receipt_schema", "biological_target", "ordered_samples",
        "provider_part_counts", "read_only_action_scope",
    )
    _exact_keys(binding, keys, "amended live-capture authorization")
    if binding["authorization_path"] != AMENDED_AUTHORIZATION_PATH or binding["acceptance_path"] != AMENDED_AUTHORIZATION_ACCEPTANCE_PATH:
        raise ExactFourPlanError("amended live-capture authorization path-role drift")
    if binding["deployment_manifest_path"] != DEPLOYMENT_MANIFEST_PATH:
        raise ExactFourPlanError("amended live-capture manifest path drift")
    for prefix in ("authorization", "acceptance", "deployment_manifest"):
        _require_int(binding[f"{prefix}_size_bytes"], f"{prefix} size", minimum=1)
        _require_sha(binding[f"{prefix}_sha256"], f"{prefix} SHA-256")
    if binding["deployment_bundle_id"] != manifest["bundle_id"]:
        raise ExactFourPlanError("amended live-capture bundle drift")
    if binding["deployment_manifest_sha256"] != digest_value(manifest):
        raise ExactFourPlanError("amended live-capture manifest digest drift")
    if binding["provider_receipt_schema"] != PROVIDER_SCHEMA or binding["biological_target"] != LOCUS:
        raise ExactFourPlanError("amended live-capture schema/target drift")
    if binding["ordered_samples"] != list(TARGETS):
        raise ExactFourPlanError("amended live-capture target order drift")
    if binding["provider_part_counts"] != {sample: count for sample, count in zip(TARGETS, PART_COUNTS)}:
        raise ExactFourPlanError("amended live-capture part-count drift")
    if binding["read_only_action_scope"] != AMENDED_READ_ONLY_SCOPE:
        raise ExactFourPlanError("amended live-capture action-scope drift")
    if revalidate_files:
        _stable_reopen_pinned_file(binding["authorization_path"], binding["authorization_size_bytes"], binding["authorization_sha256"], "amended live-capture authorization")
        _stable_reopen_pinned_file(binding["acceptance_path"], binding["acceptance_size_bytes"], binding["acceptance_sha256"], "amended live-capture authorization acceptance")
        manifest_payload = _stable_reopen_pinned_file(binding["deployment_manifest_path"], binding["deployment_manifest_size_bytes"], binding["deployment_manifest_sha256"], "amended deployment manifest")
        if canonical_json_bytes(manifest) != manifest_payload:
            raise ExactFourPlanError("amended deployment manifest bytes differ from validated manifest")
    return copy.deepcopy(binding)


def _validate_provider_authority_bindings(
    value: object,
    manifest: Mapping[str, object],
    *,
    revalidate_external_files: bool,
) -> dict[str, object]:
    bindings = _require_dict(value, "provider authority bindings")
    _exact_keys(bindings, ("expected_provider_metadata_authority", "aws_tool_authority", "amended_live_capture_authorization"), "provider authority bindings")
    metadata = _require_dict(bindings["expected_provider_metadata_authority"], "expected-provider authority binding")
    _exact_keys(metadata, ("path", "size_bytes", "sha256", "acceptance_path", "acceptance_size_bytes", "acceptance_sha256", "deployment_bundle_id", "capture_authorization_sha256"), "expected-provider authority binding")
    if (metadata["size_bytes"], metadata["sha256"], metadata["acceptance_size_bytes"], metadata["acceptance_sha256"]) != (5801, EXPECTED_PROVIDER_AUTHORITY_SHA256, 2089, EXPECTED_PROVIDER_ACCEPTANCE_SHA256):
        raise ExactFourPlanError("expected-provider authority identity drift")
    if metadata["deployment_bundle_id"] != "exact4-bundle-cee202a9f2a520062573fe50" or metadata["capture_authorization_sha256"] != LEGACY_CAPTURE_AUTHORIZATION_SHA256:
        raise ExactFourPlanError("expected-provider legacy provenance drift")
    for path_key in ("path", "acceptance_path"):
        path = Path(_require_str(metadata[path_key], f"expected-provider {path_key}"))
        if not path.is_absolute() or str(path) != os.path.normpath(str(path)):
            raise ExactFourPlanError("expected-provider authority path is not normalized absolute")

    aws = _require_dict(bindings["aws_tool_authority"], "AWS authority binding")
    _exact_keys(aws, ("path", "size_bytes", "sha256", "acceptance_path", "acceptance_size_bytes", "acceptance_sha256", "normalized_absolute_path", "executable_sha256", "version_output", "module_name", "cluster_host_observation_class", "read_only_resident_receipt_sha256", "capture_authorization_sha256"), "AWS authority binding")
    if (aws["size_bytes"], aws["sha256"], aws["acceptance_size_bytes"], aws["acceptance_sha256"]) != (1165, AWS_TOOL_AUTHORITY_SHA256, 1855, AWS_TOOL_ACCEPTANCE_SHA256):
        raise ExactFourPlanError("AWS authority file identity drift")
    if aws["executable_sha256"] != "226b41b30cf5707f262ec4087a8eda83283a3477b59ab5bd6786059f6b1d0662" or aws["version_output"] != "aws-cli/1.16.308 Python/3.9.9 Linux/3.10.0-862.el7.x86_64 botocore/1.13.44" or aws["module_name"] != "awscli/1.16" or aws["cluster_host_observation_class"] != "login.pax.tufts.edu":
        raise ExactFourPlanError("AWS executable/version/module/host authority drift")
    if aws["read_only_resident_receipt_sha256"] != RESIDENT_REHASH_RECEIPT_SHA256 or aws["capture_authorization_sha256"] != LEGACY_CAPTURE_AUTHORIZATION_SHA256:
        raise ExactFourPlanError("AWS authority provenance drift")
    executable_path = Path(_require_str(aws["normalized_absolute_path"], "AWS normalized path"))
    if not executable_path.is_absolute() or str(executable_path) != os.path.normpath(str(executable_path)):
        raise ExactFourPlanError("AWS executable authority path is not normalized absolute")
    for path_key in ("path", "acceptance_path"):
        path = Path(_require_str(aws[path_key], f"AWS {path_key}"))
        if not path.is_absolute() or str(path) != os.path.normpath(str(path)):
            raise ExactFourPlanError("AWS authority path is not normalized absolute")
    amended = _validate_amended_live_capture_authorization(bindings["amended_live_capture_authorization"], manifest, revalidate_files=revalidate_external_files)
    return {"expected_provider_metadata_authority": copy.deepcopy(metadata), "aws_tool_authority": copy.deepcopy(aws), "amended_live_capture_authorization": amended}


def _validate_provider_receipt(
    value: object,
    manifest: Mapping[str, object],
    *,
    revalidate_external_files: bool = True,
) -> dict[str, object]:
    receipt = _require_dict(value, "provider receipt")
    _exact_keys(receipt, ("schema_version", "deployment_bundle_id", "objects", "provider_tool_receipt", "capture_metadata", "authority_bindings", "action_gates"), "provider receipt")
    if receipt["schema_version"] != PROVIDER_SCHEMA or receipt["deployment_bundle_id"] != manifest["bundle_id"]:
        raise ExactFourPlanError("provider receipt schema/bundle drift")
    tool = _require_dict(receipt["provider_tool_receipt"], "provider tool receipt")
    _exact_keys(tool, ("path", "size_bytes", "sha256", "version_output", "version_output_sha256"), "provider tool receipt")
    if not Path(_require_str(tool["path"], "provider tool path")).is_absolute():
        raise ExactFourPlanError("provider tool path must be absolute")
    _require_int(tool["size_bytes"], "provider tool size", minimum=1)
    _require_sha(tool["sha256"], "provider tool SHA-256")
    _require_sha(tool["version_output_sha256"], "provider version output digest")
    if tool["version_output_sha256"] != hashlib.sha256((_require_str(tool["version_output"], "provider version output") + "\n").encode("utf-8")).hexdigest():
        raise ExactFourPlanError("provider version output digest drift")
    if type(receipt["objects"]) is not list or len(receipt["objects"]) != 13:
        raise ExactFourPlanError("provider receipt must contain exactly 13 objects")
    expected = manifest["expected_provider_objects"]
    for index, (row_value, expected_row) in enumerate(zip(receipt["objects"], expected)):
        row = _require_dict(row_value, f"provider object {index}")
        _exact_keys(row, ("sample_id", "part_index", "uri", "version_id", "etag", "content_length_bytes", "last_modified", "observation_time", "provider_command", "response_digest"), f"provider object {index}")
        for key in ("sample_id", "part_index", "uri", "version_id", "etag", "content_length_bytes", "last_modified"):
            if row[key] != expected_row[key]:
                raise ExactFourPlanError(f"provider object identity/order drift at global index {index}: {key}")
        _require_str(row["observation_time"], "provider observation_time")
        uri = str(row["uri"])
        if not uri.startswith("s3://") or "/" not in uri[5:]:
            raise ExactFourPlanError("provider URI cannot form exact HEAD command")
        bucket, key = uri[5:].split("/", 1)
        expected_command = [tool["path"], "s3api", "head-object", "--bucket", bucket, "--key", key, "--version-id", row["version_id"], "--no-sign-request"]
        if row["provider_command"] != expected_command:
            raise ExactFourPlanError("provider command is not the exact ordered versioned HEAD-only command")
        _require_sha(row["response_digest"], "provider response digest")
    capture = _require_dict(receipt["capture_metadata"], "provider capture metadata")
    _exact_keys(capture, ("observation_host", "capture_label", "object_body_downloaded"), "provider capture metadata")
    if re.fullmatch(r"login-[0-9]+\.cluster\.example", _require_str(capture["observation_host"], "provider observation host")) is None:
        raise ExactFourPlanError("provider observation host is not an accepted kernel login FQDN")
    if capture["capture_label"] != PROVIDER_CAPTURE_LABEL:
        raise ExactFourPlanError("provider capture label drift")
    _require_bool(capture["object_body_downloaded"], False, "object_body_downloaded")
    bindings = _validate_provider_authority_bindings(receipt["authority_bindings"], manifest, revalidate_external_files=revalidate_external_files)
    if tool["path"] != bindings["aws_tool_authority"]["normalized_absolute_path"] or tool["sha256"] != bindings["aws_tool_authority"]["executable_sha256"] or tool["version_output"] != bindings["aws_tool_authority"]["version_output"]:
        raise ExactFourPlanError("provider tool receipt differs from AWS authority binding")
    gates = _require_dict(receipt["action_gates"], "provider action gates")
    _exact_keys(gates, tuple(ACTION_GATES), "provider action gates")
    if gates != ACTION_GATES or any(type(gates[key]) is not bool or gates[key] is not False for key in ACTION_GATES):
        raise ExactFourPlanError("provider action gates must remain exactly false")
    return receipt


def _resident_identity(receipt: Mapping[str, object]) -> dict[str, object]:
    return {key: value for key, value in receipt.items() if key not in {"capture_metadata", "storage"}}


def _provider_identity(receipt: Mapping[str, object]) -> dict[str, object]:
    objects = []
    for row in receipt["objects"]:
        objects.append({key: value for key, value in row.items() if key != "observation_time"})
    return {"schema_version": receipt["schema_version"], "deployment_bundle_id": receipt["deployment_bundle_id"], "objects": objects, "provider_tool_receipt": receipt["provider_tool_receipt"], "authority_bindings": receipt["authority_bindings"], "action_gates": receipt["action_gates"]}


def _sample_haplotypes(resident: Mapping[str, object], sample: str) -> list[str]:
    return [row["haplotype"] for row in resident["assemblies"] if row["sample_id"] == sample]


def _run_paths(run_id: str, resident: Mapping[str, object]) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    root = f"{RUN_ROOT}/{run_id}"
    preparation: list[dict[str, object]] = []
    alignment: list[dict[str, object]] = []
    analysis: list[dict[str, object]] = []
    global_part = 0
    for sample_index, (sample, part_count) in enumerate(zip(TARGETS, PART_COUNTS)):
        sample_root = f"{root}/samples/{sample}"
        prep_root = f"{sample_root}/preparation"
        haps = _sample_haplotypes(resident, sample)
        if len(haps) != 2 or len(set(haps)) != 2:
            raise ExactFourPlanError(f"sample {sample} must have exactly two ordered haplotypes")
        prep_destinations = [
            f"{prep_root}/combined_reference.fa", f"{prep_root}/combined_reference.fa.fai",
            f"{prep_root}/combined_reference.mmi", f"{prep_root}/core.bed",
            f"{prep_root}/full_window.bed", f"{prep_root}/bait.fa",
            f"{prep_root}/PREPARED.json",
        ]
        preparation.append({
            "index": sample_index, "sample_id": sample, "locus": LOCUS, "haplotypes": haps,
            "claim": f"{root}/claims/prepare/{sample_index}.claim.json",
            "staging_root": f"{root}/staging/prepare/{sample_index}",
            "bundle_destination": prep_root,
            "destinations": prep_destinations,
        })
        for local_part in range(part_count):
            object_row = next(row for row in _provider_manifest_rows() if row["sample_id"] == sample and row["part_index"] == local_part)
            part_root = f"{sample_root}/parts/{local_part:02d}"
            alignment.append({
                "index": global_part, "sample_index": sample_index, "sample_id": sample,
                "part_index": local_part, "locus": LOCUS, "provider_object": object_row,
                "claim": f"{root}/claims/part/{global_part}.claim.json",
                "staging_root": f"{root}/staging/part/{global_part}",
                "node_scratch_basename": f"{run_id}.{sample}.{local_part:02d}.source.bam",
                "unfiltered_bam": f"{part_root}/complete_unfiltered_alignment.bam",
                "primary_bam": f"{part_root}/primary_full_window.bam",
                "primary_bai": f"{part_root}/primary_full_window.bam.bai",
                "qname_inventory": f"{part_root}/aligned_primary_qnames.txt",
                "source_qname_inventory": f"{part_root}/source_primary_qnames.txt",
                "decoded_qname_inventory": f"{part_root}/decoded_qnames.txt",
                "aligned_primary_qname_inventory": f"{part_root}/aligned_primary_qnames.txt",
                "supplementary_qname_inventory": f"{part_root}/supplementary_qnames.txt",
                "receipt": f"{part_root}/PART.json",
                "destinations": [
                    f"{part_root}/complete_unfiltered_alignment.bam",
                    f"{part_root}/primary_full_window.bam",
                    f"{part_root}/primary_full_window.bam.bai",
                    f"{part_root}/source_primary_qnames.txt",
                    f"{part_root}/decoded_qnames.txt",
                    f"{part_root}/aligned_primary_qnames.txt",
                    f"{part_root}/supplementary_qnames.txt",
                    f"{part_root}/PART.json",
                ],
            })
            global_part += 1
        analysis_root = f"{sample_root}/analysis"
        analysis.append({
            "index": sample_index, "sample_id": sample, "locus": LOCUS,
            "expected_global_part_indices": [row["index"] for row in alignment if row["sample_id"] == sample],
            "claim": f"{root}/claims/analysis/{sample_index}.claim.json",
            "staging_root": f"{root}/staging/analysis/{sample_index}",
            "merged_primary_bam": f"{analysis_root}/merged_primary_full_window.bam",
            "merged_primary_bai": f"{analysis_root}/merged_primary_full_window.bam.bai",
            "raw_depth": f"{analysis_root}/raw.depth.tsv",
            "mapq10_depth": f"{analysis_root}/mapq10.depth.tsv",
            "summary": f"{analysis_root}/depth_summary.json",
            "kcon_paf": f"{analysis_root}/type2_KCON.paf",
            "raw_pdf": f"{analysis_root}/raw.depth.pdf",
            "mapq10_pdf": f"{analysis_root}/mapq10.depth.pdf",
            "receipt": f"{analysis_root}/ANALYSIS.json",
            "done": f"{analysis_root}/DONE.json",
            "destinations": [
                f"{analysis_root}/merged_primary_full_window.bam",
                f"{analysis_root}/merged_primary_full_window.bam.bai",
                f"{analysis_root}/raw.depth.tsv", f"{analysis_root}/mapq10.depth.tsv",
                f"{analysis_root}/depth_summary.json", f"{analysis_root}/type2_KCON.paf",
                f"{analysis_root}/raw.depth.pdf", f"{analysis_root}/mapq10.depth.pdf",
                f"{analysis_root}/ANALYSIS.json", f"{analysis_root}/DONE.json",
            ],
        })
    if len(alignment) != 13:
        raise ExactFourPlanError("internal exact-four part map is not 13")
    return preparation, alignment, analysis


def build_target_only_plan(
    deployment_manifest: dict[str, object],
    resident_receipt: dict[str, object],
    provider_receipt: dict[str, object],
) -> dict[str, object]:
    manifest = _validate_deployment_manifest(deployment_manifest)
    resident = _validate_resident_receipt(resident_receipt, manifest)
    provider = _validate_provider_receipt(provider_receipt, manifest)
    immutable = {
        "deployment_manifest_sha256": digest_value(manifest),
        "bundle_id": manifest["bundle_id"],
        **{key: manifest[key] for key in CONTRACT_BINDINGS},
        "target_table_sha256": TARGET_TABLE_SHA256,
        "resident_identity": _resident_identity(resident),
        "provider_identity": _provider_identity(provider),
        "targets": manifest["targets"],
        "bundle_files": manifest["bundle_files"],
        "resources": manifest["resources"],
        "scientific_laws": manifest["scientific_laws"],
        "coverage_ratio_schema_identity": manifest["coverage_ratio_schema_identity"],
        "coverage_ratio_schema_sha256": manifest["coverage_ratio_schema_sha256"],
        "output_schemas": manifest["output_schemas"],
    }
    run_id = "exact4-" + digest_value(immutable)[:24]
    preparation, alignment, analysis = _run_paths(run_id, resident)
    entrypoint_hashes = {
        row["role"]: row["sha256"] for row in manifest["bundle_files"]
        if row["role"] in {"preparation_entrypoint", "alignment_entrypoint", "analysis_entrypoint", "planner_and_release_gate", "plotter"}
    }
    plan: dict[str, object] = {
        "schema_version": PLAN_SCHEMA,
        "run_id": run_id,
        "run_root": f"{RUN_ROOT}/{run_id}",
        "plan_digest": "",
        "deployment_manifest_sha256": digest_value(manifest),
        "bundle_id": manifest["bundle_id"],
        "resident_receipt_sha256": digest_value(resident),
        "provider_receipt_sha256": digest_value(provider),
        "immutable_run_material": immutable,
        "targets": manifest["targets"],
        "preparation_map": preparation,
        "alignment_map": alignment,
        "analysis_map": analysis,
        "entrypoint_hashes": entrypoint_hashes,
        "resources": manifest["resources"],
        "scientific_laws": manifest["scientific_laws"],
        "coverage_ratio_schema_identity": manifest["coverage_ratio_schema_identity"],
        "coverage_ratio_schema_sha256": manifest["coverage_ratio_schema_sha256"],
        **{key: manifest[key] for key in CONTRACT_BINDINGS},
        "commands": {
            "decode": ["samtools", "fasta", "-F", "0x900"],
            "align": ["minimap2", "-ax", "map-ont", "--secondary=no"],
            "primary_view": ["samtools", "view", "-u", "-F", "0x900", "-L", "<full_window.bed>", "<unfiltered_alignment.bam>"],
            "raw_depth": ["samtools", "depth", "-a", "-g", "0x600", "-b", "<full_window.bed>", "<merged_primary_full_window.bam>"],
            "mapq10_depth": ["samtools", "depth", "-a", "-g", "0x600", "-Q", "10", "-b", "<full_window.bed>", "<merged_primary_full_window.bam>"],
            "kcon": ["minimap2", "-c", "-p", "0.1", "-N", "10", "<full-window.fa>", "<type2_KCON.fa>"],
        },
        "submission_recipe": {
            "inert_only": True,
            "phases": [
                {"phase": "preparation", "array": "0-3%2", "dependency": None},
                {"phase": "alignment", "array": "0-12%2", "dependency": "afterok:preparation"},
                {"phase": "analysis", "array": "0-3%2", "dependency": "afterok:alignment"},
            ],
        },
        "forbidden_features_absent": [
            "other_sample", "other_locus", "legacy_result_path", "stale_ten_target_plan",
            "source_missing_six", "bait_ascertainment", "minimum_overlap_selection",
            "read_cap", "downsampling", "top_n", "denominator_shrinkage",
            "per_haplotype_double_alignment", "fallback", "automatic_retry",
            "cancellation", "resubmission", "cleanup", "read_inferred_copy_number",
            "ratio_category_threshold", "coverage_compatibility_category",
            "coverage_discordance_category", "ratio_derived_status_or_color",
            "legacy_ratio_bound",
        ],
        "output_schemas": manifest["output_schemas"],
        "action_gates": dict(PLAN_ACTION_GATES),
        "dry_run_only": True,
    }
    without_digest = dict(plan)
    without_digest.pop("plan_digest")
    plan["plan_digest"] = digest_value(without_digest)
    return validate_target_only_plan(plan)


def validate_target_only_plan(plan: dict[str, object]) -> dict[str, object]:
    value = _require_dict(plan, "target-only plan")
    keys = (
        "schema_version", "run_id", "run_root", "plan_digest", "deployment_manifest_sha256",
        "bundle_id", "resident_receipt_sha256", "provider_receipt_sha256", "immutable_run_material",
        "targets", "preparation_map", "alignment_map", "analysis_map", "entrypoint_hashes",
        "resources", "scientific_laws", "coverage_ratio_schema_identity", "coverage_ratio_schema_sha256", "commands",
        *CONTRACT_BINDINGS.keys(),
        "submission_recipe", "forbidden_features_absent", "output_schemas", "action_gates", "dry_run_only",
    )
    _exact_keys(value, keys, "target-only plan")
    if value["schema_version"] != PLAN_SCHEMA:
        raise ExactFourPlanError("target-only plan schema drift")
    if not re.fullmatch(r"exact4-[0-9a-f]{24}", _require_str(value["run_id"], "run_id")):
        raise ExactFourPlanError("invalid run_id")
    if value["run_id"] != "exact4-" + digest_value(value["immutable_run_material"])[:24]:
        raise ExactFourPlanError("run_id does not match immutable run material")
    expected_root = f"{RUN_ROOT}/{value['run_id']}"
    if value["run_root"] != expected_root:
        raise ExactFourPlanError("run root is not the exact isolated run path")
    if value["targets"] != [{"sample_id": s, "locus": LOCUS, "part_count": c} for s, c in zip(TARGETS, PART_COUNTS)]:
        raise ExactFourPlanError("plan target scope/order/count drift")
    if type(value["preparation_map"]) is not list or [row["index"] for row in value["preparation_map"]] != list(range(4)):
        raise ExactFourPlanError("preparation mapping drift")
    if type(value["alignment_map"]) is not list or [row["index"] for row in value["alignment_map"]] != list(range(13)):
        raise ExactFourPlanError("alignment mapping drift")
    if type(value["analysis_map"]) is not list or [row["index"] for row in value["analysis_map"]] != list(range(4)):
        raise ExactFourPlanError("analysis mapping drift")
    if [row["sample_id"] for row in value["alignment_map"]] != [sample for sample, count in zip(TARGETS, PART_COUNTS) for _ in range(count)]:
        raise ExactFourPlanError("alignment sample/part cardinality drift")
    for collection in (value["preparation_map"], value["alignment_map"], value["analysis_map"]):
        for row in collection:
            for key in ("claim", "staging_root"):
                _require_plan_path(value, row[key], f"mapping {key}")
            for destination in row["destinations"]:
                _require_plan_path(value, destination, "mapping destination")
    commands = _require_dict(value["commands"], "commands")
    _exact_keys(commands, ("decode", "align", "primary_view", "raw_depth", "mapq10_depth", "kcon"), "commands")
    if commands["align"] != ["minimap2", "-ax", "map-ont", "--secondary=no"]:
        raise ExactFourPlanError("complete-unfiltered competitive alignment command drift")
    if commands["primary_view"][:6] != ["samtools", "view", "-u", "-F", "0x900", "-L"]:
        raise ExactFourPlanError("primary-only view law drift")
    if commands["raw_depth"][:5] != ["samtools", "depth", "-a", "-g", "0x600"]:
        raise ExactFourPlanError("raw depth must explicitly include QCFAIL/DUP primaries with -g 0x600")
    if commands["mapq10_depth"][:7] != ["samtools", "depth", "-a", "-g", "0x600", "-Q", "10"]:
        raise ExactFourPlanError("MAPQ10 must include QCFAIL/DUP primaries and use uppercase -Q 10")
    if value["coverage_ratio_schema_identity"] != COVERAGE_RATIO_SCHEMA_IDENTITY or value["coverage_ratio_schema_sha256"] != COVERAGE_RATIO_SCHEMA_SHA256:
        raise ExactFourPlanError("plan coverage-ratio schema identity drift")
    immutable = _require_dict(value["immutable_run_material"], "immutable run material")
    _exact_keys(immutable, ("deployment_manifest_sha256", "bundle_id", *CONTRACT_BINDINGS.keys(), "target_table_sha256", "resident_identity", "provider_identity", "targets", "bundle_files", "resources", "scientific_laws", "coverage_ratio_schema_identity", "coverage_ratio_schema_sha256", "output_schemas"), "immutable run material")
    for key in ("deployment_manifest_sha256", "bundle_id", "targets", "resources", "scientific_laws", "coverage_ratio_schema_identity", "coverage_ratio_schema_sha256", "output_schemas"):
        if immutable[key] != value[key]:
            raise ExactFourPlanError(f"plan/immutable identity drift: {key}")
    if immutable["target_table_sha256"] != TARGET_TABLE_SHA256:
        raise ExactFourPlanError("immutable target table identity drift")
    resident_identity = _require_dict(immutable["resident_identity"], "immutable resident identity")
    _exact_keys(resident_identity, ("schema_version", "deployment_bundle_id", "deployment_manifest_sha256", "resident_files", "assemblies", "element_identities", "bed_authorities", "combined_row_sets", "assembly_encoded_counts", "sample_bed_authorities", "tools", "tool_fixture_receipt", "resources"), "immutable resident identity")
    _validate_tool_receipts(resident_identity["tools"])
    validate_resident_tool_fixture_receipt(resident_identity["tool_fixture_receipt"], resident_identity["tools"])
    provider_identity = _require_dict(immutable["provider_identity"], "immutable provider identity")
    _exact_keys(provider_identity, ("schema_version", "deployment_bundle_id", "objects", "provider_tool_receipt", "authority_bindings", "action_gates"), "immutable provider identity")
    if provider_identity["schema_version"] != PROVIDER_SCHEMA or provider_identity["deployment_bundle_id"] != value["bundle_id"]:
        raise ExactFourPlanError("immutable provider schema/bundle drift")
    amended_binding = _require_dict(provider_identity["authority_bindings"], "immutable provider authority bindings").get("amended_live_capture_authorization")
    amended_binding = _require_dict(amended_binding, "immutable amended live-capture authorization")
    current_manifest_payload = _stable_reopen_pinned_file(
        DEPLOYMENT_MANIFEST_PATH,
        amended_binding.get("deployment_manifest_size_bytes"),
        amended_binding.get("deployment_manifest_sha256"),
        "amended deployment manifest",
    )
    current_manifest = _validate_deployment_manifest(_require_dict(strict_json_bytes(current_manifest_payload, "amended deployment manifest"), "amended deployment manifest"))
    if digest_value(current_manifest) != value["deployment_manifest_sha256"]:
        raise ExactFourPlanError("current amended deployment manifest differs from immutable plan")
    _validate_provider_authority_bindings(provider_identity["authority_bindings"], current_manifest, revalidate_external_files=True)
    if provider_identity["action_gates"] != ACTION_GATES or any(type(provider_identity["action_gates"].get(key)) is not bool or provider_identity["action_gates"].get(key) is not False for key in ACTION_GATES):
        raise ExactFourPlanError("immutable provider action gates drift")
    expected_preparation, expected_alignment, expected_analysis = _run_paths(value["run_id"], resident_identity)
    if value["preparation_map"] != expected_preparation or value["alignment_map"] != expected_alignment or value["analysis_map"] != expected_analysis:
        raise ExactFourPlanError("plan phase maps differ from deterministic authenticated mapping")
    bundle_files = immutable["bundle_files"]
    if type(bundle_files) is not list:
        raise ExactFourPlanError("immutable bundle files must be an exact list")
    derived_entrypoints = {
        row["role"]: row["sha256"] for row in bundle_files
        if type(row) is dict and row.get("role") in {"preparation_entrypoint", "alignment_entrypoint", "analysis_entrypoint", "planner_and_release_gate", "plotter"}
    }
    if value["entrypoint_hashes"] != derived_entrypoints or len(derived_entrypoints) != 5:
        raise ExactFourPlanError("plan entrypoint hashes differ from immutable bundle files")
    for key, expected in CONTRACT_BINDINGS.items():
        if value[key] != expected or immutable.get(key) != expected:
            raise ExactFourPlanError(f"plan/immutable contract binding drift: {key}")
    plan_text = canonical_json_bytes(value).lower()
    for phrase in (b"coverage-compatible", b"coverage-discordant", b"compatibility_tolerance", b"ratio_min", b"ratio_max"):
        if phrase in plan_text:
            raise ExactFourPlanError("forbidden ratio-derived category/tolerance entered plan")
    gates = _require_dict(value["action_gates"], "plan action gates")
    if gates != PLAN_ACTION_GATES or any(type(gates[key]) is not bool or gates[key] is not False for key in PLAN_ACTION_GATES):
        raise ExactFourPlanError("plan action/performed gates must remain exactly false")
    _require_bool(value["dry_run_only"], True, "dry_run_only")
    without_digest = dict(value)
    observed = _require_sha(without_digest.pop("plan_digest"), "plan_digest")
    if observed != digest_value(without_digest):
        raise ExactFourPlanError("plan_digest mismatch")
    return value


def _require_plan_path(plan: Mapping[str, object], path_value: object, label: str) -> Path:
    path_text = _require_str(path_value, label)
    if ".." in Path(path_text).parts:
        raise ExactFourPlanError(f"{label} contains lexical '..'")
    path = Path(path_text)
    root = Path(_require_str(plan["run_root"], "run_root"))
    if not path.is_absolute():
        raise ExactFourPlanError(f"{label} must be absolute")
    try:
        if os.path.commonpath((str(root), str(path))) != str(root) or path == root:
            raise ExactFourPlanError(f"{label} is outside the exact run root")
    except ValueError as error:
        raise ExactFourPlanError(f"{label} is outside the run path domain") from error
    return path


def _phase_maps(plan: Mapping[str, object]) -> dict[str, object]:
    return {
        "preparation": [{"index": row["index"], "sample_id": row["sample_id"]} for row in plan["preparation_map"]],
        "alignment": [{"index": row["index"], "sample_id": row["sample_id"], "part_index": row["part_index"]} for row in plan["alignment_map"]],
        "analysis": [{"index": row["index"], "sample_id": row["sample_id"]} for row in plan["analysis_map"]],
    }


def _allowed_claims(plan: Mapping[str, object]) -> list[dict[str, object]]:
    claims = []
    for kind, key in (("prepare", "preparation_map"), ("part", "alignment_map"), ("analysis", "analysis_map")):
        claims.extend({"kind": kind, "index": row["index"], "path": row["claim"]} for row in plan[key])
    return claims


def _allowed_destinations(plan: Mapping[str, object]) -> list[str]:
    output: list[str] = []
    for key in ("preparation_map", "alignment_map", "analysis_map"):
        for row in plan[key]:
            if key == "preparation_map":
                output.append(row["bundle_destination"])
            output.extend(row["destinations"])
    if len(output) != len(set(output)):
        raise ExactFourPlanError("plan contains duplicate publication destinations")
    return output


def _validate_current_source_bytes(resident_identity: object, kind: str, sample_id: str) -> None:
    """Rehash only the phase-relevant resident bytes immediately before use."""
    identity = _require_dict(resident_identity, "immutable resident identity")
    sample = _require_str(sample_id, "current-source sample_id")
    if sample not in TARGETS:
        raise ExactFourPlanError("current-source sample_id is outside the exact four targets")
    if kind == "prepare":
        assemblies = identity.get("assemblies")
        if type(assemblies) is not list:
            raise ExactFourPlanError("immutable resident assemblies must be an exact list")
        selected = [row for row in assemblies if type(row) is dict and row.get("sample_id") == sample]
        expected_pairs = [(candidate, haplotype) for candidate, haplotype, *_ in ASSEMBLY_AUTHORITIES if candidate == sample]
        if len(selected) != 2:
            raise ExactFourPlanError(f"current-source preparation requires exactly two assemblies for {sample}")
        for index, (row_value, expected_pair) in enumerate(zip(selected, expected_pairs)):
            row = _require_dict(row_value, f"current-source assembly {index}")
            _exact_keys(row, ("sample_id", "haplotype", "fasta", "fai", "contigs"), f"current-source assembly {index}")
            if (row["sample_id"], row["haplotype"]) != expected_pair:
                raise ExactFourPlanError(f"current-source assembly identity/order drift for {sample}")
            for label in ("fasta", "fai"):
                record = _require_dict(row[label], f"current-source assembly {index} {label}")
                _exact_keys(record, ("path", "size_bytes", "sha256"), f"current-source assembly {index} {label}")
                path = Path(_require_str(record["path"], f"current-source assembly {index} {label} path"))
                expected_size = _require_int(record["size_bytes"], f"current-source assembly {index} {label} size", minimum=1)
                expected_sha = _require_sha(record["sha256"], f"current-source assembly {index} {label} SHA-256")
                observed = _fingerprint(path)
                if observed["size_bytes"] != expected_size or observed["sha256"] != expected_sha:
                    raise ExactFourPlanError(f"current assembly {label} bytes drifted: {sample}/{row['haplotype']}")
        return
    if kind == "analysis":
        resident_files = identity.get("resident_files")
        if type(resident_files) is not list:
            raise ExactFourPlanError("immutable resident files must be an exact list")
        matches = [row for row in resident_files if type(row) is dict and row.get("role") == "type2_kcon"]
        if len(matches) != 1:
            raise ExactFourPlanError("current-source analysis requires exactly one type2_kcon authority")
        row = _require_dict(matches[0], "current-source type2_kcon")
        _exact_keys(row, ("role", "path", "size_bytes", "sha256"), "current-source type2_kcon")
        path = Path(_require_str(row["path"], "current-source type2_kcon path"))
        expected_size = _require_int(row["size_bytes"], "current-source type2_kcon size", minimum=1)
        expected_sha = _require_sha(row["sha256"], "current-source type2_kcon SHA-256")
        observed = _fingerprint(path)
        if observed["size_bytes"] != expected_size or observed["sha256"] != expected_sha:
            raise ExactFourPlanError("current resident file bytes drifted: type2_kcon")
        return
    raise ExactFourPlanError("current-source kind must be prepare or analysis")


def _validate_current_execution_identity(plan: Mapping[str, object], release_tools: object) -> None:
    """Freshly rehash every executable byte admitted to a write-capable path."""
    tools = _validate_tool_receipts(release_tools)
    expected_tools = plan["immutable_run_material"]["resident_identity"]["tools"]
    if tools != expected_tools:
        raise ExactFourPlanError("execution release tool receipts do not equal immutable resident tools")
    for row in tools:
        observed = _fingerprint(Path(row["path"]))
        if observed["size_bytes"] != row["size_bytes"] or observed["sha256"] != row["sha256"]:
            raise ExactFourPlanError(f"current tool bytes drifted: {row['name']}")

    module_path = Path(__file__)
    if module_path.is_symlink():
        raise ExactFourPlanError("planner entrypoint may not be a symlink")
    bundle_root = module_path.resolve(strict=True).parent.parent
    roles = {"preparation_entrypoint", "alignment_entrypoint", "analysis_entrypoint", "planner_and_release_gate", "plotter"}
    observed_hashes: dict[str, str] = {}
    for row_value in plan["immutable_run_material"]["bundle_files"]:
        row = _require_dict(row_value, "immutable bundle file")
        if row["role"] not in roles:
            continue
        path = bundle_root / row["path"]
        fingerprint = _fingerprint(path)
        mode = f"{stat.S_IMODE(path.stat().st_mode):04o}"
        if fingerprint["size_bytes"] != row["size_bytes"] or fingerprint["sha256"] != row["sha256"] or mode != row["mode"]:
            raise ExactFourPlanError(f"current entrypoint bytes/mode drifted: {row['role']}")
        observed_hashes[row["role"]] = fingerprint["sha256"]
    if set(observed_hashes) != roles or observed_hashes != plan["entrypoint_hashes"]:
        raise ExactFourPlanError("current entrypoint identity set differs from plan")


def validate_execution_release(
    release: dict[str, object],
    plan: dict[str, object],
    plan_file_sha256: str,
) -> dict[str, object]:
    plan_value = validate_target_only_plan(plan)
    value = _require_dict(release, "execution release")
    keys = (
        "schema_version", "plan_file_sha256", "plan_digest", "run_id", "bundle_id",
        "deployment_manifest_sha256", "resident_receipt_sha256", "provider_receipt_sha256",
        *CONTRACT_BINDINGS.keys(), "coverage_ratio_schema_identity", "coverage_ratio_schema_sha256",
        "no_drift_receipt_sha256", "no_drift_identity_sha256", "entrypoint_hashes",
        "tool_receipts", "phase_maps", "resources", "allowed_phases", "one_attempt",
        "attempt_number", "allowed_claims", "allowed_destinations", "execution_authorized",
        "production_authorized", "copy_number_inference_authorized", "model_use_authorized",
        "manuscript_use_authorized",
    )
    _exact_keys(value, keys, "execution release")
    if value["schema_version"] != RELEASE_SCHEMA:
        raise ExactFourPlanError("execution release schema drift")
    expected = {
        "plan_file_sha256": _require_sha(plan_file_sha256, "current plan-file SHA-256"),
        "plan_digest": plan_value["plan_digest"],
        "run_id": plan_value["run_id"],
        "bundle_id": plan_value["bundle_id"],
        "deployment_manifest_sha256": plan_value["deployment_manifest_sha256"],
        "resident_receipt_sha256": plan_value["resident_receipt_sha256"],
        "provider_receipt_sha256": plan_value["provider_receipt_sha256"],
        **{key: plan_value[key] for key in CONTRACT_BINDINGS},
        "coverage_ratio_schema_identity": plan_value["coverage_ratio_schema_identity"],
        "coverage_ratio_schema_sha256": plan_value["coverage_ratio_schema_sha256"],
        "entrypoint_hashes": plan_value["entrypoint_hashes"],
        "phase_maps": _phase_maps(plan_value),
        "resources": plan_value["resources"],
        "allowed_phases": ["preparation", "alignment", "analysis"],
        "allowed_claims": _allowed_claims(plan_value),
        "allowed_destinations": _allowed_destinations(plan_value),
    }
    for key, expected_value in expected.items():
        if value[key] != expected_value:
            raise ExactFourPlanError(f"execution release does not bind current {key}")
    _require_sha(value["no_drift_receipt_sha256"], "no-drift receipt SHA-256")
    _require_sha(value["no_drift_identity_sha256"], "no-drift identity SHA-256")
    _validate_current_execution_identity(plan_value, value["tool_receipts"])
    _require_bool(value["one_attempt"], True, "one_attempt")
    if value["attempt_number"] != 1 or type(value["attempt_number"]) is not int:
        raise ExactFourPlanError("execution release attempt_number must be exact integer 1")
    _require_bool(value["execution_authorized"], True, "execution_authorized")
    for key in ("production_authorized", "copy_number_inference_authorized", "model_use_authorized", "manuscript_use_authorized"):
        _require_bool(value[key], False, key)
    return value


def _load_context(plan_path: Path, release_path: Path) -> tuple[dict[str, object], dict[str, object], str]:
    plan = _load_json_file(plan_path, "target-only plan")
    release = _load_json_file(release_path, "execution release")
    plan_sha = sha256_file(plan_path)
    validate_execution_release(release, plan, plan_sha)
    return plan, release, plan_sha


def _mapping(plan: Mapping[str, object], kind: str, index: int) -> dict[str, object]:
    if kind == "prepare":
        collection = plan["preparation_map"]
    elif kind == "part":
        collection = plan["alignment_map"]
    elif kind == "analysis":
        collection = plan["analysis_map"]
    else:
        raise ExactFourPlanError("claim kind must be prepare, part, or analysis")
    if type(index) is not int or index < 0:
        raise ExactFourPlanError("claim index must be a nonnegative exact integer")
    matches = [row for row in collection if row["index"] == index]
    if len(matches) != 1:
        raise ExactFourPlanError(f"claim mapping is not exact for {kind}/{index}")
    return matches[0]


def _assert_existing_components_safe(root: Path, candidate: Path) -> None:
    if ".." in candidate.parts:
        raise ExactFourPlanError("path contains lexical '..'")
    try:
        if os.path.commonpath((str(root), str(candidate))) != str(root):
            raise ExactFourPlanError("path is outside exact run root")
    except ValueError as error:
        raise ExactFourPlanError("path domain differs from run root") from error
    current = Path(root.anchor)
    for component in root.parts[1:]:
        current /= component
        if os.path.lexists(current):
            info = os.lstat(current)
            if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
                raise ExactFourPlanError(f"unsafe existing run-root component: {current}")
    relative = candidate.relative_to(root)
    current = root
    for component in relative.parts[:-1]:
        current /= component
        if os.path.lexists(current):
            info = os.lstat(current)
            if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
                raise ExactFourPlanError(f"unsafe existing parent component: {current}")


def _mkdir_parents_create_only(
    plan_path: Path,
    release_path: Path,
    destination: Path,
) -> None:
    plan, _, _ = _load_context(plan_path, release_path)
    root = Path(plan["run_root"])
    _assert_existing_components_safe(root, destination)
    chain: list[Path] = []
    current = destination.parent
    while current != current.parent and not os.path.lexists(current):
        chain.append(current)
        current = current.parent
    for directory in reversed(chain):
        fresh_plan, _, _ = _load_context(plan_path, release_path)
        _require_plan_path(fresh_plan, str(destination), "mutation destination")
        _assert_existing_components_safe(Path(fresh_plan["run_root"]), destination)
        try:
            os.mkdir(directory, 0o750)
        except FileExistsError:
            info = os.lstat(directory)
            if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
                raise ExactFourPlanError(f"concurrent unsafe parent appeared: {directory}")


def _create_claim(plan_path: Path, release_path: Path, kind: str, index: int, token: str) -> dict[str, object]:
    if not SAFE_TOKEN_RE.fullmatch(_require_str(token, "claim token")):
        raise ExactFourPlanError("claim token is not a safe exact token")
    plan, _, plan_sha = _load_context(plan_path, release_path)
    mapping = _mapping(plan, kind, index)
    if kind in ("prepare", "analysis"):
        _validate_current_source_bytes(
            plan["immutable_run_material"]["resident_identity"], kind, mapping["sample_id"]
        )
    claim_path = _require_plan_path(plan, mapping["claim"], "claim path")
    _mkdir_parents_create_only(plan_path, release_path, claim_path)
    claim: dict[str, object] = {
        "schema_version": CLAIM_SCHEMA, "run_id": plan["run_id"], "plan_file_sha256": plan_sha,
        "plan_digest": plan["plan_digest"], "bundle_id": plan["bundle_id"], "kind": kind,
        "index": index, "token": token, "claim_digest": "",
    }
    material = dict(claim)
    material.pop("claim_digest")
    claim["claim_digest"] = digest_value(material)
    payload = canonical_json_bytes(claim)
    fresh_plan, _, fresh_sha = _load_context(plan_path, release_path)
    if fresh_sha != plan_sha or _mapping(fresh_plan, kind, index)["claim"] != str(claim_path):
        raise ExactFourPlanError("claim context drifted immediately before creation")
    _assert_existing_components_safe(Path(fresh_plan["run_root"]), claim_path)
    if os.path.lexists(claim_path):
        raise ExactFourPlanError(f"claim path already exists and blocks one attempt: {claim_path}")
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(claim_path, flags, 0o440)
    try:
        os.write(descriptor, payload)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    staging_root = _require_plan_path(plan, mapping["staging_root"], "claim staging root")
    _mkdir_parents_create_only(plan_path, release_path, staging_root)
    fresh_plan, _, fresh_sha = _load_context(plan_path, release_path)
    if fresh_sha != plan_sha or _mapping(fresh_plan, kind, index)["staging_root"] != str(staging_root):
        raise ExactFourPlanError("claim staging context drifted before creation")
    _assert_existing_components_safe(Path(fresh_plan["run_root"]), staging_root)
    if os.path.lexists(staging_root):
        raise ExactFourPlanError(f"staging root already exists and blocks one attempt: {staging_root}")
    try:
        os.mkdir(staging_root, 0o750)
    except FileExistsError as error:
        raise ExactFourPlanError(f"concurrent staging root appeared and blocks one attempt: {staging_root}") from error
    return claim


def _rename_noreplace(source: Path, destination: Path) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    if renameat2 is None:
        raise ExactFourPlanError("Linux renameat2(RENAME_NOREPLACE) is unavailable")
    renameat2.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
    renameat2.restype = ctypes.c_int
    result = renameat2(-100, os.fsencode(source), -100, os.fsencode(destination), 1)
    if result != 0:
        error = ctypes.get_errno()
        raise ExactFourPlanError(f"no-replace directory publication failed: {os.strerror(error)}")


def _publish(plan_path: Path, release_path: Path, staged: Path, destination: Path) -> None:
    plan, release, _ = _load_context(plan_path, release_path)
    destination = _require_plan_path(plan, str(destination), "publication destination")
    if str(destination) not in release["allowed_destinations"]:
        raise ExactFourPlanError("publication destination is not in exact release mapping")
    allowed_staging = [Path(row["staging_root"]) for key in ("preparation_map", "alignment_map", "analysis_map") for row in plan[key]]
    if not staged.is_absolute() or not any(os.path.commonpath((str(root), str(staged))) == str(root) and staged != root for root in allowed_staging):
        raise ExactFourPlanError("staged path is outside exact plan staging roots")
    if not os.path.lexists(staged) or staged.is_symlink():
        raise ExactFourPlanError("staged artifact is missing or symlinked")
    staged_info = os.lstat(staged)
    if not (stat.S_ISREG(staged_info.st_mode) or stat.S_ISDIR(staged_info.st_mode)):
        raise ExactFourPlanError("staged artifact must be a regular file or directory")
    if stat.S_ISREG(staged_info.st_mode) and staged_info.st_nlink != 1:
        raise ExactFourPlanError("staged regular file has a hard-link alias")
    _mkdir_parents_create_only(plan_path, release_path, destination)
    fresh_plan, fresh_release, _ = _load_context(plan_path, release_path)
    if str(destination) not in fresh_release["allowed_destinations"]:
        raise ExactFourPlanError("destination release mapping drifted before publication")
    _assert_existing_components_safe(Path(fresh_plan["run_root"]), destination)
    if os.path.lexists(destination):
        raise ExactFourPlanError("destination already exists, including dangling symlink")
    if os.lstat(staged).st_dev != os.lstat(destination.parent).st_dev:
        raise ExactFourPlanError("cross-filesystem publication is forbidden")
    if staged.is_dir():
        _rename_noreplace(staged, destination)
    else:
        os.link(staged, destination, follow_symlinks=False)
        _load_context(plan_path, release_path)
        staged.unlink()


def inspect_qname_inventory(path: os.PathLike[str] | str) -> dict[str, object]:
    item = Path(path)
    if not item.is_file() or item.is_symlink():
        raise ExactFourPlanError(f"QNAME inventory is not a regular file: {item}")
    payload = item.read_bytes()
    if payload and not payload.endswith(b"\n"):
        raise ExactFourPlanError("QNAME inventory must be LF terminated")
    rows = payload.splitlines()
    if any(not row or b"\t" in row or b"\r" in row or any(byte < 33 or byte > 126 for byte in row) for row in rows):
        raise ExactFourPlanError("QNAME inventory contains an invalid exact ASCII name")
    if rows != sorted(rows) or len(rows) != len(set(rows)):
        raise ExactFourPlanError("QNAME inventory must be bytewise sorted and unique")
    return {"path": str(item), "size_bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest(), "count": len(rows), "encoding": "SAM-QNAME ASCII bytes, LF terminated", "order": "LC_ALL=C bytewise ascending, unique"}


def verify_qname_inventory_union(paths: Sequence[os.PathLike[str] | str]) -> dict[str, object]:
    if not paths:
        raise ExactFourPlanError("at least one QNAME inventory is required")
    seen: set[bytes] = set()
    inventories = []
    for path in paths:
        info = inspect_qname_inventory(path)
        inventories.append(info)
        for name in Path(path).read_bytes().splitlines():
            if name in seen:
                raise ExactFourPlanError(f"cross-part duplicate QNAME: {name!r}")
            seen.add(name)
    union_bytes = b"".join(name + b"\n" for name in sorted(seen))
    return {"inputs": inventories, "count": len(seen), "sha256": hashlib.sha256(union_bytes).hexdigest(), "encoding": "SAM-QNAME ASCII bytes, LF terminated", "order": "LC_ALL=C bytewise ascending, unique"}


def validate_complete_qname_ownership(
    source_path: os.PathLike[str] | str,
    decoded_path: os.PathLike[str] | str,
    aligned_primary_path: os.PathLike[str] | str,
    supplementary_path: os.PathLike[str] | str,
) -> dict[str, object]:
    """Validate the exact production-worker QNAME bijection/subset seam."""
    paths = [Path(value) for value in (source_path, decoded_path, aligned_primary_path, supplementary_path)]
    rows = [inspect_qname_inventory(path) for path in paths]
    payloads = [path.read_bytes() for path in paths]
    if payloads[0] != payloads[1] or payloads[0] != payloads[2]:
        raise ExactFourPlanError("source, decoded, and aligned-primary QNAME inventories differ")
    source_names = set(payloads[0].splitlines())
    supplementary_names = set(payloads[3].splitlines())
    if not supplementary_names.issubset(source_names):
        raise ExactFourPlanError("supplementary QNAME inventory is not a subset of complete molecules")
    return {
        "source": rows[0], "decoded": rows[1], "aligned_primary": rows[2],
        "supplementary": rows[3], "qname_set_cardinality": rows[0]["count"],
        "exact_qname_set_and_multiplicity_equal": True,
        "supplementary_qnames_subset": True,
    }


def _validate_artifact_rows(rows: object, allowed: Sequence[str], label: str) -> list[dict[str, object]]:
    if type(rows) is not list or not rows:
        raise ExactFourPlanError(f"{label} artifacts must be a nonempty exact list")
    output = []
    seen = set()
    for index, value in enumerate(rows):
        row = _require_dict(value, f"{label} artifact {index}")
        _exact_keys(row, ("path", "size_bytes", "sha256"), f"{label} artifact {index}")
        path = _require_str(row["path"], f"{label} artifact path")
        if path not in allowed or path in seen:
            raise ExactFourPlanError(f"{label} has foreign or duplicate artifact path")
        seen.add(path)
        allow_empty = path.endswith("/supplementary_qnames.txt")
        expected_size = _require_int(row["size_bytes"], f"{label} artifact size", minimum=0 if allow_empty else 1)
        expected_sha = _require_sha(row["sha256"], f"{label} artifact SHA-256")
        observed = _fingerprint(Path(path), allow_empty=allow_empty)
        if observed["size_bytes"] != expected_size or observed["sha256"] != expected_sha:
            raise ExactFourPlanError(f"{label} artifact fingerprint drift: {path}")
        output.append(row)
    return output


def _validate_closure_receipt(
    plan: Mapping[str, object],
    mapping: Mapping[str, object],
    path: str,
    schema: str,
    kind: str,
) -> dict[str, object]:
    receipt = _load_json_file(path, f"{kind} receipt")
    keys = (
        "schema_version", "kind", "run_id", "plan_digest", "bundle_id", "index",
        "sample_id", "locus", "claim_token", "artifacts", "scientific_assertions",
        "details", "read_inferred_copy_number", "receipt_digest",
    )
    _exact_keys(receipt, keys, f"{kind} receipt")
    if receipt["schema_version"] != schema or receipt["kind"] != kind:
        raise ExactFourPlanError(f"{kind} receipt schema/kind drift")
    bindings = {"run_id": plan["run_id"], "plan_digest": plan["plan_digest"], "bundle_id": plan["bundle_id"], "index": mapping["index"], "sample_id": mapping["sample_id"], "locus": LOCUS}
    for key, expected in bindings.items():
        if receipt[key] != expected:
            raise ExactFourPlanError(f"{kind} receipt binding drift: {key}")
    if receipt["read_inferred_copy_number"] != "not_estimated":
        raise ExactFourPlanError("read_inferred_copy_number must remain literal not_estimated")
    assertions = _require_dict(receipt["scientific_assertions"], f"{kind} scientific assertions")
    required_assertions = {
        "diagnostic_only": True,
        "complete_unfiltered_alignment_required": True,
        "competitive_combined_reference_required": True,
        "primary_only_depth_required": True,
        "read_cap_or_downsampling_allowed": False,
    }
    if assertions != required_assertions:
        raise ExactFourPlanError(f"{kind} scientific assertions drift")
    details = _require_dict(receipt["details"], f"{kind} details")
    if kind == "prepare":
        _exact_keys(details, ("haplotypes", "core_bed_sha256", "full_window_bed_sha256", "combined_reference_sha256", "combined_fai_sha256", "combined_mmi_sha256", "bait_sha256", "publication_mode"), "prepare details")
        if details["haplotypes"] != mapping["haplotypes"] or details["publication_mode"] != "same_filesystem_directory_noreplace":
            raise ExactFourPlanError("preparation details drift")
        bed_rows = [row for row in plan["immutable_run_material"]["resident_identity"]["sample_bed_authorities"] if row["sample_id"] == mapping["sample_id"]]
        if len(bed_rows) != 1 or details["core_bed_sha256"] != bed_rows[0]["core_bed"]["sha256"] or details["full_window_bed_sha256"] != bed_rows[0]["full_window_bed"]["sha256"]:
            raise ExactFourPlanError("prepared BED bytes differ from authenticated coordinate authority")
    elif kind == "part":
        _exact_keys(details, ("provider_object", "download_size_bytes", "download_sha256", "decoded_molecules", "decoded_bases", "primary_records", "supplementary_records", "qname_inventory_sha256", "source_qname_inventory_sha256", "decoded_qname_inventory_sha256", "aligned_primary_qname_inventory_sha256", "supplementary_qname_inventory_sha256", "qname_set_cardinality", "exact_qname_set_and_multiplicity_equal", "supplementary_qnames_subset", "exclusion_mask_hex", "exclusion_mask_decimal", "decode_argv", "alignment_argv", "primary_view_argv"), "part details")
        if details["provider_object"] != mapping["provider_object"] or details["exclusion_mask_hex"] != "0x900" or details["exclusion_mask_decimal"] != 2304:
            raise ExactFourPlanError("part provider/mask details drift")
        for count_key in ("download_size_bytes", "decoded_molecules", "decoded_bases", "primary_records", "supplementary_records", "qname_set_cardinality"):
            _require_int(details[count_key], f"part {count_key}")
        _require_sha(details["download_sha256"], "download SHA-256")
        for key in ("qname_inventory_sha256", "source_qname_inventory_sha256", "decoded_qname_inventory_sha256", "aligned_primary_qname_inventory_sha256", "supplementary_qname_inventory_sha256"):
            _require_sha(details[key], key)
        _require_bool(details["exact_qname_set_and_multiplicity_equal"], True, "exact QNAME set and multiplicity equality")
        _require_bool(details["supplementary_qnames_subset"], True, "supplementary QNAME subset")
        inventories = {name: inspect_qname_inventory(mapping[name]) for name in ("source_qname_inventory", "decoded_qname_inventory", "aligned_primary_qname_inventory", "supplementary_qname_inventory")}
        ownership = validate_complete_qname_ownership(
            mapping["source_qname_inventory"], mapping["decoded_qname_inventory"],
            mapping["aligned_primary_qname_inventory"], mapping["supplementary_qname_inventory"],
        )
        source_rows = Path(mapping["source_qname_inventory"]).read_bytes().splitlines()
        decoded_rows = Path(mapping["decoded_qname_inventory"]).read_bytes().splitlines()
        aligned_rows = Path(mapping["aligned_primary_qname_inventory"]).read_bytes().splitlines()
        supplementary_rows = Path(mapping["supplementary_qname_inventory"]).read_bytes().splitlines()
        if not source_rows or ownership["exact_qname_set_and_multiplicity_equal"] is not True or ownership["supplementary_qnames_subset"] is not True:
            raise ExactFourPlanError("part retained QNAME inventories violate exact ownership")
        if details["qname_set_cardinality"] != len(source_rows) or details["decoded_molecules"] != len(decoded_rows) or details["primary_records"] != len(aligned_rows) or details["supplementary_records"] < len(supplementary_rows):
            raise ExactFourPlanError("part QNAME counts differ from retained exact inventories")
        if details["qname_inventory_sha256"] != inventories["aligned_primary_qname_inventory"]["sha256"]:
            raise ExactFourPlanError("part legacy QNAME digest differs from aligned inventory")
        for name in inventories:
            if details[f"{name}_sha256"] != inventories[name]["sha256"]:
                raise ExactFourPlanError(f"part retained QNAME digest drift: {name}")
        if details["alignment_argv"][:4] != ["minimap2", "-ax", "map-ont", "--secondary=no"] or details["decode_argv"][:4] != ["samtools", "fasta", "-F", "0x900"]:
            raise ExactFourPlanError("part decode/alignment argv drift")
    elif kind == "analysis":
        _exact_keys(details, ("expected_global_part_indices", "qname_union", "merged_primary_bam_sha256", "raw_depth_sha256", "mapq10_depth_sha256", "summary_sha256", "kcon_paf_sha256", "raw_pdf_sha256", "mapq10_pdf_sha256", "raw_depth_argv", "mapq10_depth_argv", "kcon_argv", "kcon_selected_annotations", "deployment_manifest_sha256", "entrypoint_hashes", *CONTRACT_BINDINGS.keys(), "coverage_ratio_schema_identity", "coverage_ratio_schema_sha256", "done_written_last"), "analysis details")
        if details["expected_global_part_indices"] != mapping["expected_global_part_indices"] or details["done_written_last"] is not True:
            raise ExactFourPlanError("analysis part map/DONE ordering detail drift")
        if details["raw_depth_argv"][:5] != ["samtools", "depth", "-a", "-g", "0x600"] or details["mapq10_depth_argv"][:7] != ["samtools", "depth", "-a", "-g", "0x600", "-Q", "10"]:
            raise ExactFourPlanError("analysis raw/MAPQ10 argv or explicit QCFAIL/DUP inclusion drift")
        if details["deployment_manifest_sha256"] != plan["deployment_manifest_sha256"] or details["entrypoint_hashes"] != plan["entrypoint_hashes"]:
            raise ExactFourPlanError("analysis deployment/entrypoint identity drift")
        try:
            import plot_exact4_depth
        except ImportError as error:
            raise ExactFourPlanError(f"cannot import frozen KCON validator: {error}") from error
        try:
            expected_annotations = plot_exact4_depth._select_kcon_annotations_from_plan(dict(plan), int(mapping["index"]), Path(str(mapping["kcon_paf"])))
        except plot_exact4_depth.ExactFourPlotError as error:
            raise ExactFourPlanError(f"cannot revalidate plan-bound KCON annotations: {error}") from error
        if details["kcon_selected_annotations"] != expected_annotations:
            raise ExactFourPlanError("analysis KCON selected annotations drift from current PAF/plan authority")
        if details["coverage_ratio_schema_identity"] != COVERAGE_RATIO_SCHEMA_IDENTITY or details["coverage_ratio_schema_sha256"] != COVERAGE_RATIO_SCHEMA_SHA256:
            raise ExactFourPlanError("analysis coverage-ratio schema identity drift")
        for key, expected in CONTRACT_BINDINGS.items():
            if details[key] != expected:
                raise ExactFourPlanError(f"analysis contract binding drift: {key}")
    _validate_artifact_rows(receipt["artifacts"], mapping["destinations"], kind)
    material = dict(receipt)
    observed_digest = _require_sha(material.pop("receipt_digest"), f"{kind} receipt digest")
    if observed_digest != digest_value(material):
        raise ExactFourPlanError(f"{kind} receipt digest mismatch")
    return receipt


def validate_prepared_sample(plan: dict[str, object], sample_index: int) -> dict[str, object]:
    value = validate_target_only_plan(plan)
    mapping = _mapping(value, "prepare", sample_index)
    return _validate_closure_receipt(value, mapping, mapping["destinations"][-1], PREPARED_SCHEMA, "prepare")


def validate_alignment_part(plan: dict[str, object], part_index: int) -> dict[str, object]:
    value = validate_target_only_plan(plan)
    mapping = _mapping(value, "part", part_index)
    receipt = _validate_closure_receipt(value, mapping, mapping["receipt"], PART_SCHEMA, "part")
    for key in ("source_qname_inventory", "decoded_qname_inventory", "aligned_primary_qname_inventory", "supplementary_qname_inventory"):
        qname_path = mapping[key]
        info = inspect_qname_inventory(qname_path)
        if not any(row["path"] == qname_path and row["sha256"] == info["sha256"] for row in receipt["artifacts"]):
            raise ExactFourPlanError(f"part receipt does not bind current {key}")
    return receipt


def validate_analysis_sample(plan: dict[str, object], sample_index: int) -> dict[str, object]:
    value = validate_target_only_plan(plan)
    mapping = _mapping(value, "analysis", sample_index)
    for part_index in mapping["expected_global_part_indices"]:
        validate_alignment_part(value, part_index)
    qnames = [value["alignment_map"][part_index]["qname_inventory"] for part_index in mapping["expected_global_part_indices"]]
    verify_qname_inventory_union(qnames)
    receipt = _validate_closure_receipt(value, mapping, mapping["receipt"], ANALYSIS_SCHEMA, "analysis")
    try:
        import plot_exact4_depth
    except ImportError as error:
        raise ExactFourPlanError(f"cannot import frozen analysis validator: {error}") from error
    try:
        expected_summary, expected_annotations = plot_exact4_depth._build_summary_and_kcon_annotations_from_plan(
            value,
            sample_index,
            Path(mapping["raw_depth"]),
            Path(mapping["mapq10_depth"]),
            Path(mapping["kcon_paf"]),
        )
        observed_summary = _load_json_file(mapping["summary"], "coverage-ratio summary")
        if observed_summary != expected_summary or receipt["details"]["kcon_selected_annotations"] != expected_annotations:
            raise ExactFourPlanError("analysis summary/KCON closure differs from current plan-bound evidence")
        plot_exact4_depth.validate_pdf(Path(mapping["raw_pdf"]), {"summary": observed_summary, "metric": "raw", "kcon_annotations": expected_annotations})
        plot_exact4_depth.validate_pdf(Path(mapping["mapq10_pdf"]), {"summary": observed_summary, "metric": "mapq10", "kcon_annotations": expected_annotations})
    except plot_exact4_depth.ExactFourPlotError as error:
        raise ExactFourPlanError(f"analysis scientific/PDF closure failed: {error}") from error
    done = _load_json_file(mapping["done"], "DONE receipt")
    _exact_keys(done, ("schema_version", "run_id", "plan_digest", "bundle_id", "deployment_manifest_sha256", "entrypoint_hashes", "sample_id", "analysis_receipt_sha256", *CONTRACT_BINDINGS.keys(), "coverage_ratio_schema_identity", "coverage_ratio_schema_sha256", "read_inferred_copy_number", "done_digest"), "DONE receipt")
    if done["schema_version"] != "hml2_7p22_exact4_done_1" or done["run_id"] != value["run_id"] or done["plan_digest"] != value["plan_digest"] or done["bundle_id"] != value["bundle_id"] or done["deployment_manifest_sha256"] != value["deployment_manifest_sha256"] or done["entrypoint_hashes"] != value["entrypoint_hashes"] or done["sample_id"] != mapping["sample_id"]:
        raise ExactFourPlanError("DONE receipt identity drift")
    if done["analysis_receipt_sha256"] != sha256_file(mapping["receipt"]) or done["read_inferred_copy_number"] != "not_estimated":
        raise ExactFourPlanError("DONE receipt does not bind current diagnostic analysis")
    if done["coverage_ratio_schema_identity"] != COVERAGE_RATIO_SCHEMA_IDENTITY or done["coverage_ratio_schema_sha256"] != COVERAGE_RATIO_SCHEMA_SHA256:
        raise ExactFourPlanError("DONE coverage-ratio schema identity drift")
    for key, expected in CONTRACT_BINDINGS.items():
        if done[key] != expected:
            raise ExactFourPlanError(f"DONE contract binding drift: {key}")
    material = dict(done)
    if _require_sha(material.pop("done_digest"), "DONE digest") != digest_value(material):
        raise ExactFourPlanError("DONE digest mismatch")
    return receipt


def validate_run_closure(plan: dict[str, object]) -> dict[str, object]:
    value = validate_target_only_plan(plan)
    done_rows = []
    for index, sample in enumerate(TARGETS):
        validate_analysis_sample(value, index)
        path = value["analysis_map"][index]["done"]
        done_rows.append({"sample_id": sample, "path": path, "sha256": sha256_file(path)})
    closure: dict[str, object] = {
        "schema_version": RUN_CLOSURE_SCHEMA, "run_id": value["run_id"],
        "plan_digest": value["plan_digest"], "bundle_id": value["bundle_id"],
        "deployment_manifest_sha256": value["deployment_manifest_sha256"],
        "entrypoint_hashes": value["entrypoint_hashes"],
        "targets": list(TARGETS), "done_receipts": done_rows,
        **CONTRACT_BINDINGS,
        "coverage_ratio_schema_identity": copy.deepcopy(COVERAGE_RATIO_SCHEMA_IDENTITY),
        "coverage_ratio_schema_sha256": COVERAGE_RATIO_SCHEMA_SHA256,
        "read_inferred_copy_number": "not_estimated", "closure_digest": "",
    }
    material = dict(closure)
    material.pop("closure_digest")
    closure["closure_digest"] = digest_value(material)
    return closure


def _parse_cli(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    sub = parser.add_subparsers(dest="command", required=True)
    freeze = sub.add_parser("freeze-bundle", allow_abbrev=False)
    freeze.add_argument("--bundle-root", type=Path, required=True)
    capture = sub.add_parser("capture-resident", allow_abbrev=False)
    capture.add_argument("--deployment-manifest", required=True)
    capture.add_argument("--tool-fixture-receipt", type=Path, required=True)
    for name in ("bash", "python", "aws", "minimap2", "samtools", "sbatch"):
        capture.add_argument(f"--{name}", required=True)
    fixtures = sub.add_parser("run-fixtures", allow_abbrev=False)
    fixtures.add_argument("--deployment-manifest", required=True)
    fixtures.add_argument("--tools-receipt", type=Path, required=True)
    fixtures.add_argument("--fixture-authorization", required=True)
    fixtures.add_argument("--fixture-authorization-size-bytes", type=int, required=True)
    fixtures.add_argument("--fixture-authorization-sha256", required=True)
    fixtures.add_argument("--fixture-authorization-acceptance", required=True)
    fixtures.add_argument("--fixture-authorization-acceptance-size-bytes", type=int, required=True)
    fixtures.add_argument("--fixture-authorization-acceptance-sha256", required=True)
    validate_fixtures = sub.add_parser("validate-fixtures", allow_abbrev=False)
    validate_fixtures.add_argument("--fixture-receipt", type=Path, required=True)
    validate_fixtures.add_argument("--tools-receipt", type=Path, required=True)
    qnames = sub.add_parser("validate-qname-ownership", allow_abbrev=False)
    qnames.add_argument("--source", type=Path, required=True)
    qnames.add_argument("--decoded", type=Path, required=True)
    qnames.add_argument("--aligned-primary", type=Path, required=True)
    qnames.add_argument("--supplementary", type=Path, required=True)
    render = sub.add_parser("render-plan", allow_abbrev=False)
    render.add_argument("--deployment-manifest", type=Path, required=True)
    render.add_argument("--resident-receipt", type=Path, required=True)
    render.add_argument("--provider-receipt", type=Path, required=True)
    validate = sub.add_parser("validate-plan", allow_abbrev=False)
    validate.add_argument("--plan", type=Path, required=True)
    validate_release = sub.add_parser("validate-release", allow_abbrev=False)
    validate_release.add_argument("--release", type=Path, required=True)
    validate_release.add_argument("--plan", type=Path, required=True)
    claim = sub.add_parser("claim", allow_abbrev=False)
    claim.add_argument("--plan", type=Path, required=True)
    claim.add_argument("--execution-release", type=Path, required=True)
    claim.add_argument("--kind", choices=("prepare", "part", "analysis"), required=True)
    claim.add_argument("--index", type=int, required=True)
    claim.add_argument("--token", required=True)
    prepared = sub.add_parser("validate-prepared", allow_abbrev=False)
    prepared.add_argument("--plan", type=Path, required=True)
    prepared.add_argument("--sample-index", type=int, required=True)
    part = sub.add_parser("validate-part", allow_abbrev=False)
    part.add_argument("--plan", type=Path, required=True)
    part.add_argument("--part-index", type=int, required=True)
    analysis = sub.add_parser("validate-analysis", allow_abbrev=False)
    analysis.add_argument("--plan", type=Path, required=True)
    analysis.add_argument("--sample-index", type=int, required=True)
    run = sub.add_parser("validate-run", allow_abbrev=False)
    run.add_argument("--plan", type=Path, required=True)
    publish = sub.add_parser("publish", allow_abbrev=False)
    publish.add_argument("--plan", type=Path, required=True)
    publish.add_argument("--execution-release", type=Path, required=True)
    publish.add_argument("--staged", type=Path, required=True)
    publish.add_argument("--destination", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_cli(argv)
    try:
        if args.command == "freeze-bundle":
            output = build_deployment_manifest(args.bundle_root)
        elif args.command == "run-fixtures":
            manifest = _load_fixed_deployment_manifest_argument(args.deployment_manifest, "run-fixtures deployment manifest")
            tools_wrapper = _load_json_file(args.tools_receipt, "resident tools receipt")
            _exact_keys(tools_wrapper, ("tools",), "resident tools receipt")
            output = run_resident_tool_fixtures(
                tools_wrapper["tools"], manifest,
                args.fixture_authorization, args.fixture_authorization_size_bytes, args.fixture_authorization_sha256,
                args.fixture_authorization_acceptance, args.fixture_authorization_acceptance_size_bytes,
                args.fixture_authorization_acceptance_sha256,
            )
        elif args.command == "validate-fixtures":
            tools_wrapper = _load_json_file(args.tools_receipt, "resident tools receipt")
            _exact_keys(tools_wrapper, ("tools",), "resident tools receipt")
            output = validate_resident_tool_fixture_receipt(_load_json_file(args.fixture_receipt, "resident fixture receipt"), tools_wrapper["tools"])
        elif args.command == "capture-resident":
            manifest = _load_fixed_deployment_manifest_argument(args.deployment_manifest, "capture-resident deployment manifest")
            output = capture_resident_receipt(manifest, {name: getattr(args, name) for name in ("bash", "python", "aws", "minimap2", "samtools", "sbatch")}, _load_json_file(args.tool_fixture_receipt, "resident tool fixture receipt"))
        elif args.command == "validate-qname-ownership":
            output = validate_complete_qname_ownership(args.source, args.decoded, args.aligned_primary, args.supplementary)
        elif args.command == "render-plan":
            output = build_target_only_plan(_load_json_file(args.deployment_manifest, "deployment manifest"), _load_json_file(args.resident_receipt, "resident receipt"), _load_json_file(args.provider_receipt, "provider receipt"))
        elif args.command == "validate-plan":
            output = validate_target_only_plan(_load_json_file(args.plan, "target-only plan"))
        elif args.command == "validate-release":
            plan = _load_json_file(args.plan, "target-only plan")
            output = validate_execution_release(_load_json_file(args.release, "execution release"), plan, sha256_file(args.plan))
        elif args.command == "claim":
            output = _create_claim(args.plan, args.execution_release, args.kind, args.index, args.token)
        elif args.command == "validate-prepared":
            output = validate_prepared_sample(_load_json_file(args.plan, "target-only plan"), args.sample_index)
        elif args.command == "validate-part":
            output = validate_alignment_part(_load_json_file(args.plan, "target-only plan"), args.part_index)
        elif args.command == "validate-analysis":
            output = validate_analysis_sample(_load_json_file(args.plan, "target-only plan"), args.sample_index)
        elif args.command == "validate-run":
            output = validate_run_closure(_load_json_file(args.plan, "target-only plan"))
        elif args.command == "publish":
            _publish(args.plan, args.execution_release, args.staged, args.destination)
            output = {"published": str(args.destination)}
        else:
            raise ExactFourPlanError("unknown exact-four command")
        sys.stdout.buffer.write(canonical_json_bytes(output))
        return 0
    except ExactFourPlanError as error:
        print(f"exact4_plan: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
