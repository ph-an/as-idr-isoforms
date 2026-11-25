# Alternative splicing and intrinsically disordered regions

Code to reproduce the analyses and figures for the preprint:

**Alternative splicing edits intrinsically disordered regions across human protein isoforms**  
Tommy Phan, Susie Song and Joshua M. Akey<br>
Preprint: bioRxiv (2026) https://doi.org/10.64898/2026.11.23.768762

---

## What this repository contains

Scripts and notebooks to:
- build the cohort of human protein isoforms (UniProtKB/Swiss-Prot) and their MetaPredict IDRs
- map every isoform exactly to its GENCODE transcript and project each residue onto the genome (GRCh38)
- compare isoforms by genomic codon (shared and isoform-specific positions, IDR outcomes)
- analyse disorder by exon usage, protein class, protein family, gene age and splice-site conservation
- run the robustness analyses and generate every main and Extended Data figure

## Repository layout

```text
.
├── analysis/             # manuscript analyses; each folder writes outputs/ with the numbers used in the text
│   └── manuscript_figures/   # one script per main / Extended Data figure
├── data/                 # input manifest and download instructions (large inputs are not in git)
├── manuscript/           # rendered figures (main/, extended_data/, supplementary/) and Supplementary Table 1
├── notebooks/            # cohort, IDR tables, gene age and APPRIS annotation
├── pipeline/             # UniProt -> GENCODE mapping, genomic projection, exon tables
├── scripts/              # helper scripts
├── run.py                # rebuilds every table and figure in order
└── requirements.txt
```

## Getting started

### Requirements

Python 3.13 (other recent versions should work).

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

### Reproducing the paper

1. Download the inputs listed under [Data availability](#data-availability) into `data/raw/` (exact paths and
   checksums are in [`data/manifest.tsv`](data/manifest.tsv); details in [`data/README.md`](data/README.md)).
2. Run everything, or one step at a time:

```bash
python run.py                          # inputs -> proteome -> genomic -> analysis -> paper (~40 min, plus MetaPredict for the sensitivity step)
python run.py --steps analysis paper   # analyses and figures only, once data/ is built
python run.py --list                   # describe the steps
```

| Step | What it does | Writes |
|---|---|---|
| `inputs` | checks the SHA-256 of every input | — |
| `proteome` | runs `notebooks/01–03` (cohort, IDR tables, gene age, APPRIS) | `data/annotated/` |
| `genomic` | `pipeline/run_pipeline.py`: UniProt → Ensembl mapping, GENCODE coordinates, genomic projection, exon tables | `data/genomic/` |
| `analysis` | the manuscript analyses in `analysis/` | `analysis/*/outputs/` |
| `paper` | figure scripts in `analysis/manuscript_figures/` | `manuscript/figures/` |

Small analysis outputs (summary tables and the JSON files with every number quoted in the text) are included in
`analysis/*/outputs/`, and the rendered figures in `manuscript/figures/`. Some figure scripts also read the rebuilt
tables in `data/`, so run the `proteome` and `genomic` steps once before redrawing figures.

## Data availability

Large input files are not included in this repository. Download them to the paths below.

| Input | Download | Description | Source | Version | Path |
|---|---|---|---|---|---|
| Human proteins and isoforms | [UniProt](https://www.uniprot.org/uniprotkb?query=reviewed%3Atrue+AND+model_organism%3A9606) | Swiss-Prot human entries including isoforms (FASTA, and TSV with all columns) | UniProt | 2026_01 (42,562 sequences) | `data/raw/uniprot/uniprotkb_allproteins.fasta`, `data/raw/uniprot/uniprotkb_annotated` |
| UniProt → Ensembl | [UniProt REST](https://rest.uniprot.org/) | isoform → Ensembl transcript and protein cross-references (downloaded by the pipeline if absent) | UniProt | 2026_01 | `data/raw/uniprot/uniprot_ensembl_xref_raw.tsv` |
| UniProt isoform → Ensembl gene | [UniProt ID mapping](https://www.uniprot.org/id-mapping) | used by the APPRIS notebook | UniProt | 2026_01 | `data/raw/uniprot/idmapping/upkb_to_ensemble_humanisos.tsv` |
| Disordered regions | [idrome_constructor](https://github.com/holehouse-lab/idrome_constructor) | MetaPredict v3.1 IDRs for every sequence in the FASTA above | Lotthammer et al. 2024; Emenecker et al. 2021 | MetaPredict v3.1 | `data/raw/metapredict/IDRome_AllSwissProtHumansProIsos.csv` |
| GENCODE | [GENCODE 50](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/) | comprehensive annotation, protein-coding translations and transcripts (downloaded by the pipeline if absent) | GENCODE | v50, GRCh38 | `data/raw/gencode/` |
| APPRIS | [APPRIS](https://appris.bioinfo.cnio.es/) | principal-isoform annotation | APPRIS | GENCODE-based | `data/raw/appris/appris_data.appris.txt` |
| TF list | [HumanTFs](https://humantfs.ccbr.utoronto.ca/download.php) | human transcription factors and DBD families | Lambert et al. 2018 | v1.01 | `data/raw/humantfs/DatabaseExtract_v_1.01.xlsx` |
| GenOrigin | [GenOrigin](http://chenzxlab.hzau.edu.cn/GenOrigin/#!/download) | gene age | Tong et al. 2021 | 2021 | `data/raw/genorigin/Homo_sapiens.csv` |
| Splice-site conservation | [GitHub](https://github.com/iminkin/splice-sites-conservation) / [Zenodo](https://doi.org/10.5281/zenodo.14893716) | GT/AG conservation across 405 mammalian genomes | Minkin & Salzberg 2025 | `data/processed/splice_sites.csv.gz` | `data/raw/splice_sites_conservation/splice_sites.csv.gz` |

Small derived annotation tables (TF labels and gene ages per UniProt accession) are included in `data/external/`.

## Citation

```bibtex
@article{phan_song_akey_2026_splicing_idrs,
  title   = {Alternative splicing edits intrinsically disordered regions across human protein isoforms},
  author  = {Phan, Tommy and Song, Susie and Akey, Joshua M.},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {10.64898/2026.11.23.768762}
}
```

## License

MIT (see [`LICENSE`](LICENSE)).
