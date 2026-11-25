# Methodology

Goal: place the predicted IDRs of human protein isoforms onto exact GRCh38 coordinates, so that any two isoforms of a
gene can be compared residue by residue and every coding exon carries an IDR fraction. All coordinates are GENCODE v50 /
GRCh38. This document describes `pipeline/`; the analyses built on its outputs are in `analysis/` and the paper's Methods.

**Inputs.** UniProtKB/Swiss-Prot human entries including isoforms (release 2026_01); MetaPredict v3.1 IDR intervals
(`IDRome_AllSwissProtHumansProIsos.csv`); the isoform tables built by `notebooks/01–03` (`data/annotated/`: TF family,
gene age, APPRIS); GENCODE v50 comprehensive GTF, `pc_translations` and `pc_transcripts`. Release, paths and thresholds
are set in `tfidr_pipeline/config.py` and stamped onto every output.

**Reading the counts.** Counts are for the main cohort (genes with ≥2 isoforms, `TFIDR_MIN_ISOFORMS=2`; IDRs ≥20
residues) and percentages are of its 32,579 isoforms unless stated otherwise. Nothing is deleted: an isoform that fails
a step gets an eligibility flag and is left out of the analyses that need that step. The same counts for the earlier
≥3-isoform cohort are in `manuscript/supplementary_table1_cohort_comparison.md`.

## Data flow and reproducibility

The numbered steps below are a narrative order. The code runs as stages `s1`–`s9`, four post-modules (`resolve`,
`events`, `exons`, `igv`) and the analyses `s9b`, `s9c`, `s9d` and `s10b`, driven by `pipeline/run_pipeline.py`.

| Role | Modules (step) | What it does |
|---|---|---|
| **Data acquisition** (network) | `fetch_uniprot_ensembl` (1) · `s1_download_gencode` (3, 5) | Streams UniProt cross-references; downloads the GENCODE GTF and FASTAs |
| **Coordinate extraction** | `s2_exon_cds_coords` (3) · `s4b_genomic_projection` (5) | Parses the GTF into exon and CDS tables; builds per-nucleotide genomic position arrays and BED intervals |
| **Validation / filtering** | `s3_validate_translations` (2) · `resolve_isoform_transcripts` (6) · `s8_qc_and_subsets` | Assigns sequence-match classes, resolves transcripts by coding structure, writes QC and subset tables |
| **Analysis / assembly** | `s4` (4) · `s5`, `s6` (7) · `exon_level_table` (8) · `classify_splice_events` (9) · `s9` (10) · `s7`, `make_igv_tracks` | Joins, summaries, statistics and exports |

Only two modules touch the network (`fetch_uniprot_ensembl`, `s1_download_gencode`). With the UniProt release and
GENCODE v50 pinned (both are cached under `data/raw/` and checked against `data/manifest.tsv`), every later stage is
deterministic.

---

## 0 — Starting cohort

32,579 reviewed human UniProt isoforms from 10,532 genes (10,532 canonical, 22,047 alternative), after filtering to
entries with a named gene, a length of at least 20 residues and protein-existence level below 5, restricted to genes
with at least two isoforms. MetaPredict supplies 51,627 IDRs of at least 20 residues in this cohort.

## 1 — Candidate Ensembl transcripts

UniProt's isoform-tagged Ensembl cross-references (for example `ENST…; ENSP…; ENSG…. [P04637-1]`) nominate candidate
transcripts for each isoform. A canonical accession (`P04637`) is matched to UniProt's `-1` tag; an alternative isoform
must match its exact `-N` tag. Version suffixes are kept in raw columns and stripped in normalised columns.
_Result:_ 25,340 isoforms (77.8%) have at least one isoform-specific candidate transcript; the rest have only a gene-level
link or none and are left out of transcript, exon and genomic analyses.
_Out:_ `isoform_transcript_map.csv`, `isoform_transcript_candidates_long.csv`.

## 2 — Sequence-validate each candidate

Each candidate's GENCODE v50 translation is compared with the UniProt sequence. Classes: **exact** · near_exact ·
length_match_only · mismatch · no_translation. **Only exact matches are used**; the weaker classes are kept for QC.
When an isoform has several candidates, all are scored and the best becomes the resolved transcript.
_Result:_ 24,983 isoforms (76.7%; 98.6% of those with a candidate) have an exact match.
_Out:_ `mapping_validation_scaled.parquet`.

## 3 — Exon and coding coordinates

From the comprehensive GTF, the pipeline takes genomic exon and CDS coordinates of each mapped transcript, orders coding
exons 5′→3′ with cumulative CDS offsets, and keeps CDS phase and strand so split codons and minus-strand transcripts
project correctly. The comprehensive annotation is used so that every sequence-valid transcript receives coordinates;
the extra low-support transcripts it contains only survive if they pass the exact-match (step 2) and translation
(step 5) checks.
_Out:_ `ensembl_transcript_exon_coordinates.parquet`, `transcript_cds_to_protein_coordinate_map.parquet`.

## 4 — IDR and exon overlap in protein coordinates

MetaPredict IDR intervals and their ordered complements are overlaid on the amino-acid ranges of the coding exons, giving
per-segment length, number of overlapping exons and exon boundaries. This is approximate (a codon can be split across
two exons); step 5 replaces it with the exact genomic projection.
_Out:_ `idr_exon_overlap_table.parquet`, `nonidr_exon_overlap_table.parquet`.

## 5 — Project to GRCh38 and check by translation

For each validated isoform, the coding exons of its transcript are walked in order to build a strand-aware array with
one genomic position per coding nucleotide. Residues `[start, end]` map to coding nucleotides `[(start−1)×3 : end×3]`;
slicing the array and collapsing runs gives genomic (BED) intervals. Grouping the array in threes gives each codon, keyed
in `analysis/reference_free/` by the genomic coordinate of its first base.
_Two checks:_ the full coding sequence must translate back to the UniProt protein, and every projected segment must
translate back to its protein subsequence. Failures are logged and excluded.
_Result:_ 24,817 isoforms pass the full-sequence check; all 76,781 projected IDR and ordered segments translate back
exactly.
_Out:_ `genomic_scaffold_validation.parquet`, `genomic_projection_segments.parquet`, `idr/nonidr_genomic_intervals.parquet`,
`idrs.bed`, `nonidrs.bed`, `exons.bed`, `genomic_projection_failures.csv`.

## 6 — Resolve several valid transcripts by coding structure

For each isoform with more than one exact-match transcript, a coding-structure signature (chromosome, strand, ordered
CDS-exon intervals) is built. Transcripts that differ only in untranslated regions are collapsed; isoforms whose valid
transcripts have different coding structures are flagged as ambiguous and never picked automatically.
_Result:_ **24,803 isoforms (76.1%) form the MAIN set** (one coding structure) used by every exon- and codon-level
analysis; 125 are ambiguous. The MAIN isoforms come from 10,388 genes.
_Out:_ `isoform_transcript_mapping_resolved.{csv,parquet}`.

## 7 — Isoform master table and gene summary

Mapping, validation, projection, TF annotation, gene age, APPRIS and IDR measures are joined onto all 32,579 isoforms,
with eligibility flags (`analysis_include_isoform_level`, `_exon_level`, `_genomic_projection`) and an `exclude_reason`.
A gene-level table summarises isoform counts, mapping coverage and the range of %IDR across each gene's isoforms.
_Out:_ `tf_idr_isoform_master.{csv,parquet}`, `tf_idr_gene_summary.{csv,parquet}`.

## 8 — Exon table

One row per MAIN isoform × exon of its transcript: Ensembl gene, transcript, protein and exon identifiers; genomic, CDS
and amino-acid coordinates; CDS phase; per-exon IDR fraction (the fraction of the exon's encoded residues in IDRs; ≥0.5
is mostly disordered); MANE, Ensembl-canonical and APPRIS flags; and `exon_usage`, the definition every analysis uses:
among a gene's mapped transcripts (genes with ≥2), an exon in all of them is **constitutive** and one in some but not all
is **alternative**. `is_constitutive` and `transcript_inclusion_fraction` count all GENCODE protein-coding transcripts of
the gene and are kept for the sensitivity analysis only.
_Result:_ 294,826 coding-exon rows; 124,613 unique (gene, exon) pairs. The exon analyses use the 8,196 genes with at
least two mapped transcripts: 124,984 unique exons, 62,054 alternative and 62,930 constitutive.
_Out:_ `exon_level_table.{csv,parquet}` (+ data dictionary).

## 9 — Canonical-versus-alternative splice events

For genes in which the canonical and at least one alternative isoform are in the MAIN set, the coding-exon structures are
compared and each differing exon is classified (exon skip, alternative 5′ or 3′ splice site, alternative first or last
exon, intron retention, complex) and labelled by whether it adds, removes or resizes a mostly disordered exon. These
tables describe individual events; the paper's comparison of isoforms uses all isoform pairs by genomic codon instead
(`analysis/reference_free/`).
_Out:_ `splice_events.{csv,parquet}`, `splice_event_idr_summary.csv`.

## 10 — Exon-level splicing statistics

Stage `s9` restricts to validated isoforms and writes the per-isoform × exon table with exon usage
(`exon_idr_annotation.csv`) that the analyses in `analysis/` read, plus descriptive statistics (`splicing_idr_*.csv`).
The paper's estimates (the alternative-minus-constitutive difference with genes as the unit, adjusted and within-gene
models, sensitivity analyses) are computed in `analysis/reanalysis_2026-10-26/` and `analysis/protein_classes/`.
_Out:_ `exon_idr_annotation.csv`, `splicing_idr_*.csv`.

---

## Cohort attrition (isoforms)

Each tier is a filter on the 32,579-isoform cohort; "removed" means left out of that tier, not deleted.

| Tier | Isoforms | % of cohort | Removed |
|---|---:|---:|---:|
| Cohort (genes with ≥2 isoforms) | 32,579 | 100% | — |
| ≥1 candidate transcript (step 1) | 25,340 | 77.8% | 7,239 |
| Exact sequence match (step 2) | 24,983 | 76.7% | 7,596 |
| Full-sequence translation check (step 5) | 24,817 | 76.2% | 7,762 |
| MAIN, one coding structure (step 6) | 24,803 | 76.1% | 7,776 |

## Which population to use

| Claim | Population |
|---|---|
| Protein sequence, %IDR, protein class, family, gene age | All 32,579 isoforms |
| Exon usage and exon IDR fraction | MAIN isoforms of genes with ≥2 mapped transcripts |
| Codon-level comparison of isoforms (shared and isoform-specific positions, IDR outcomes) | MAIN isoforms |
| Exact genomic coordinates | `genomic_scaffold_validated=True` |

## Key caveats

- **Match policy:** only exact sequence matches are used; `near_exact` is a QC label.
- **Comprehensive GTF:** low-support transcripts it adds (NMD, readthrough, partial) are removed by the exact-match and
  translation checks.
- **Pseudoreplication:** exons and isoform pairs are nested within genes, so the analyses use genes as the unit
  (resampling genes, gene-clustered errors, gene fixed effects).
- **Coordinates:** GTF-derived tables are 1-based inclusive; BED files are 0-based half-open.
