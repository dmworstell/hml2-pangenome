# September 2026 manuscript revision analyses

This directory reproduces the revised nucleotide comparisons, Type-I cassette analyses, Figure 3 bootstrap/sequence comparisons, and corrected terminal-boundary measurements. It does not include ORF-break dating.

From a fresh repository checkout, create a separate Python 3.12 environment:

```bash
python3.12 -m venv .venv-september-revision
source .venv-september-revision/bin/activate
python -m pip install -r analysis/september2026_revision/requirements.txt
python analysis/september2026_revision/reproduce.py --figures --out /tmp/hml2-revision-reproduction
```

The revision requirements pin the tested NumPy 2.3.5 environment. The root requirements describe the separate compact-panel workflow and pin NumPy 2.5.3. Do not install both requirements files in the same environment.

Choose a new output directory. The command never overwrites the published supplementary tables. It compares 41 regenerated tables or sequence files with their retained manuscript versions. Every table field is checked, including counts, denominators and exclusions. Floating-point values use a relative tolerance of 1e-12 and absolute tolerance of 1e-14. Row order is ignored because some original tables followed filesystem traversal order. The four complete cassette FASTA files must match byte for byte. The final report is `reproduction_checks.json` in the output directory.

Use `--only arrays`, `--only solo`, `--only type1`, `--only phylogeny`, or `--only tsd` for a single analysis. Omit `--figures` to skip the additional array, Type-I and Figure 3 rendering. The solo-LTR and TSD builders also render their figures as part of their original analysis scripts.

## Analysis and figure inventory

| Analysis | Portable scripts | Published numerical inputs/results | Final figure |
| --- | --- | --- | --- |
| Tandem-array nucleotide differences and conditional clock sensitivity | `arrays_analyze.py`, `arrays_sites.py`, `arrays_clock.py`, `arrays_figure.py` | `Supplementary_Data/Table_S15` | S5 |
| Within-locus solo-LTR diversity and illustrative substitution expectations | `solo_compare.py`, `solo_figure.py` | `Supplementary_Data/Table_S16` | S18 |
| Human cassette comparisons, LTR5Hs-matched controls, and primate alignments | `type1_analyze.py`, `type1_figures.py` | `Supplementary_Data/Table_S17` | Main 5 (panel C revised), S13, Table S11 alignment data |
| LTR/Pol bootstrap and exact-sharing network nucleotide differences | `phylogeny_analyze.py`, `phylogeny_tables.py`, `phylogeny_figure.py` | `Supplementary_Data/Table_S14/phylogeny_review` | Main 3 |
| Corrected paired terminal-boundary observations | `tsd_rebuild.py` | `Supplementary_Data/Figure_S7_boundary_evidence` | S6 |

The portable scripts output `Figure_S5_array_nucleotide_variation.*` and `Figure_S18_8q11_LTR_diversity.*`. The archived source scripts retain earlier S20/S21 filenames. Those names predate the final manuscript order. Table S15 is the source of final Figure S5, and Table S16 is the source of final Figure S18. Main Figure 5 and supplemental Figure S13 retain their final numbers. The corrected boundary tables in `Supplementary_Data/Figure_S7_boundary_evidence/` supply final Figure S6. The older `Figure_S10_TSD_per_locus.tsv` table does not supply that figure.

The complete cassette aligned FASTA files and sequence/frequency tables are retained as Table S11 supporting data under the historical path `Supplementary_Data/Table_S17/`. The separate 24-page alignment PDF is omitted from the submission. The renderer retains that historical output for reproducibility.

## Retained inputs and limits

- Arrays start from the 617 extracted copy sequences, the 12 supplied KCON alignments and copy provenance. They cover 289 arrays. They reproduce the comparisons without the original genome extraction directories or MAFFT installation. The catalog-to-local 500-base padding correction is retained in the provenance and validation tables. Repeating the upstream extraction requires the original source FASTAs identified by hashes in that provenance. The 14q11.2 NA20282 h2 source was unavailable and is still excluded.
- Solo-LTR analysis includes the sequence-bearing observation map for all 3,052 locally available source observations. It rebuilds the same sorted unique-sequence identifiers and consumes the retained minimap2 PAF to remove host flanks. It does not require the original local filenames. To independently regenerate the alignment, run `minimap2 -x asm20 -c --secondary=no inputs/solo/ltr_detection_ref.fa <output>/solo/solo_comparison_unique.fa` and replace the retained mapping in an isolated output copy. Alignment software changes can change mappings, so this is separate from reproducing the retained numerical result.
- Type-I analysis includes all seven original analysis inputs, including the 45-human representative alignment, the authoritative LTR subfamily table, and the 83-candidate primate KCON projection and ledger. It does not infer unobserved sources, expression, copackaging, or insertion ages.
- Figure 3 includes the retained sequence clusters and unrooted trees, alignment admissions and source observations, twelve KCON alignments and the fallback projected rows. It repeats all 1,000 LTR and 1,000 Pol bootstrap replicates with the original seeds. Xq28 comparisons begin at the retained, source-hashed Env sequences and MAFFT alignment in the supplementary data. Repeating the earlier native sequence extraction requires the source FASTAs identified in `Xq28_native_source_files.tsv`.
- TSD analysis repeats the exact catalog-to-boundary joins. The retained 61,936-row catalog is projected to the 11 columns this analysis reads. The paired-boundary table, its metadata, direct terminal evidence, discrepancy table, population metadata and a Fiber-seq donor-ID-only roster are included. The public Fiber-seq input contains only 39 unique `Individual_ID` values. Collaborator locus-level and per-donor actuation/coverage measurements are excluded. The catalog's historical structural/ORF annotations are not replaced by the boundary measurements. All 59,656 public-donor records remain in the denominator. The 42 length-ambiguous 8p23.1a pairs remain separate from the 20 candidates that differ across all tested 4–6-base lengths.

`input_provenance.json` records original source hashes and any column projection or historical path-label normalization. These are immutable biological-input checksums, not code-version gates. `inputs/solo/observations.tsv` also records the hash of each original extraction file. The original workspace-bound scripts in `Supplementary_Data` are retained as analysis history. Use the command above for a portable rerun.

`figure_support/render_support.py` retains only the functions and constants consumed by revised Figures 3 and 5 from the original manuscript figure builder. Rendering uses the same data and plotting logic. The alignment PDF uses Matplotlib's bundled DejaVu fonts instead of macOS-specific font files. Its sequence content and 24-page structure are retained. Font availability and plotting-library versions can alter figure pixels, so figure byte identity is not an acceptance condition. The published figures remain the visually reviewed manuscript assets. Figure 3 displays “Telomeric Type I” and “Telomeric Type II”. Historical `acro_type1` and `acro_type2` input identifiers are preserved. This follows the current terminology and final numbering in the repository’s `FIGURE_TABLE_MAP.md`.

## Validation

`validation.json` records the actual reproduction check and PDF page counts used for this release. This is validation of the bounded revised analyses from the included retained inputs. It is not a claim that the upstream genome extraction, public-read pipelines, or unrelated dating work was rerun.
