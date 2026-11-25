# Exon usage and disorder

Exon-usage analyses used in the manuscript, with robustness checks. They read the pipeline outputs
(`data/genomic/`) and re-process no raw input. `outputs/unique_exons.csv` is read by the other analyses.

```bash
python analysis/reanalysis_2026-10-26/reanalysis.py   # ~1-5 min; writes outputs/ incl. key_numbers.json
```

| ID | Question | Headline (see `outputs/key_numbers.json`) |
|---|---|---|
| R5 | Sensitivity to the exon-usage definition | Alternative exons remain more disordered when usage is defined over all GENCODE transcripts; disorder falls with inclusion fraction |
| R6 | Is IDR "exon modularity" a length artefact? | Partly. Against a same-protein, same-length placement null, IDRs contain 8–12% fewer exon boundaries than expected (ordered segments 3–10% more); single-exon encoding persists in every length bin but the raw "~3x" ratio is inflated by length |
| — | MetaPredict version (`check_metapredict_version.py`) | metapredict 3.1.01 defaults reproduce the stored IDRs for 400/400 sampled isoforms |

Unit conventions: unique exon = (UniProt gene accession, Ensembl exon ID), IDR fraction averaged over isoform rows
(as in pipeline stage s9c). Models use gene-clustered sandwich errors.
