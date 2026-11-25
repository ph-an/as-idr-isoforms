# Sensitivity to the disorder definition

Recomputes the main results (class medians, alternative minus constitutive exon disorder, the within-gene effect of an exon
with a poorly conserved splice site, and the gene-age slopes) under other disorder calls.

| Script | What it does |
|---|---|
| `predict_scores.py` | Per-residue MetaPredict v3.1 scores for every cohort isoform (about 1.5 h on two cores; restartable by chunk) |
| `predict_v1_subset.py` | MetaPredict V1 scores (the earlier MetaPredict version): first the genes with lineage-restricted exons plus as many random genes, then `all` for the rest |
| `sensitivity.py` | Builds each disorder call (stored IDRs; v3 at thresholds 0.4/0.6 and minimum lengths 10/30/50; raw v3 scores; V1; MobiDB-lite from UniProt, canonical only) and recomputes the results. Optional arguments name the calls to run. |

Outputs: `outputs/sensitivity.json`, `outputs/sensitivity_summary.csv`. The per-residue scores (`outputs/scores/`, `outputs/v1*/`)
are not committed. Recomputing the default call from the scores reproduces the stored IDRs at 99.9999% of residues.
