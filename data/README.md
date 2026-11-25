# Data

`python run.py` first checks every input against [`manifest.tsv`](manifest.tsv) (path, SHA-256, size, source,
which step reads it). Inputs are not in git (about 610 MB); [`raw/README.md`](raw/README.md) lists every file with its
size and source. Download them as below.

## Inputs (`data/raw/`)

| Path | What | How to obtain |
|---|---|---|
| `uniprot/uniprotkb_allproteins.fasta` | UniProtKB/Swiss-Prot human canonical and isoform sequences (42,562 entries) | uniprot.org, query `reviewed:true AND model_organism:9606`, FASTA including isoforms; UniProt release 2026_01 (files published 2026-01-28; downloaded 2026-06-02); 20,431 canonical and 22,131 isoform accessions |
| `uniprot/uniprotkb_annotated` | The same proteome as TSV with all return fields | uniprot.org, same query, TSV, all columns; release 2026_01 |
| `uniprot/idmapping/upkb_to_ensemble_humanisos.tsv` | UniProt isoform → Ensembl gene | UniProt ID mapping (UniProtKB AC/ID → Ensembl) |
| `uniprot/uniprot_ensembl_xref_raw.tsv` | UniProt → Ensembl transcript/protein cross-references | Downloaded by the pipeline's `fetch` step from the live UniProt REST API if absent. A new download reflects the current UniProt release and will not match the checksum, so keep this copy for exact reproduction |
| `metapredict/IDRome_AllSwissProtHumansProIsos.csv` | MetaPredict v3.1 IDR segments and ALBATROSS properties for the FASTA above (72,333 segments) | Made by the authors with the Holehouse-lab `idrome_constructor` Colab notebook. **Ground truth; not regenerated** |
| `genorigin/Homo_sapiens.csv` | GenOrigin gene ages | chenzxlab.hzau.edu.cn/GenOrigin → download → Homo sapiens |
| `appris/appris_data.appris.txt` | APPRIS per-transcript annotation | appris.bioinfo.cnio.es → downloads → Homo sapiens (GENCODE) |
| `gencode/gencode.v50.{annotation.gtf,pc_translations.fa,pc_transcripts.fa}.gz` | GENCODE release 50, GRCh38 | Downloaded by pipeline stage `s1` from `ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/` if absent |
| `humantfs/DatabaseExtract_v_1.01.xlsx` | HumanTFs v1.01 TF list (Lambert et al. 2018) | humantfs.ccbr.utoronto.ca/download.php |
| `splice_sites_conservation/splice_sites.csv.gz` | Splice-site conservation across 405 mammalian genomes (Minkin & Salzberg 2025) | github.com/iminkin/splice-sites-conservation (`data/processed/splice_sites.csv.gz`, Git LFS); Zenodo 10.5281/zenodo.14893716 |

## Frozen cluster results (`data/external/hpc_snapshot/`, in git)

TF labels (416 genes) and gene ages (5,108 proteins) derived from an Ensembl-gene→UniProt table on the lab's
computing cluster. See [`external/hpc_snapshot/README.md`](external/hpc_snapshot/README.md).

## Derived data (rebuilt by `run.py`, not in git)

| Location | Written by | Main tables |
|---|---|---|
| `data/annotated/` | notebooks 01–03 | `IDRisoforms_df.csv` (isoform master: sequence, TF status, %IDR, IDR properties), `idr_df.csv` (one row per IDR segment), `*_geneage.csv` (+ gene age, isoform counts, %IDR range), `*_appris.csv` |
| `data/genomic/` | `pipeline/run_pipeline.py` | `processed/tf_idr_isoform_master`, `processed/exon_level_table` (+ data dictionary), `analysis/exon_idr_annotation.csv`, `analysis/splicing_idr_*.csv`, `isoform_mapping/splice_events`, `qc/qc_report.md`, BED/IGV tracks |
| `analysis/*/outputs/` | `analysis/` | summary tables and the numbers quoted in the manuscript |

## Cohort

42,562 Swiss-Prot human entries → reviewed, named gene, length ≥20 aa, protein existence <5 (isoforms inherit
the canonical's PE) → genes with ≥2 isoforms (`TFIDR_MIN_ISOFORMS`, default 2): **32,579 isoforms, 10,532 genes**.
IDRs: MetaPredict segments **≥20 aa** (51,627 in the cohort). Genomic analyses use isoforms whose GENCODE v50
translation matches the UniProt sequence exactly and resolves to one coding structure (24,803 isoforms).
The earlier ≥3-isoform cohort (21,801 isoforms, 5,143 genes) is rebuilt with `TFIDR_MIN_ISOFORMS=3`.
