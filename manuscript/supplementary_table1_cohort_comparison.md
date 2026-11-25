# Supplementary Table 1. Cohort comparison: genes with ≥3 vs ≥2 annotated isoforms

The main analysis uses genes with >=2 isoforms (`TFIDR_MIN_ISOFORMS=2`, the default). The >=3 cohort was the original
analysis. Every number is read from the outputs of the same scripts run on each cohort (`TFIDR_MIN_ISOFORMS=3` reruns
the whole workflow on the ≥3 cohort). Fate percentages compare each canonical IDR with its alternative isoforms by sequence alignment.

| Quantity | >=3 isoforms | >=2 isoforms (main) |
|---|---|---|
| Isoforms | 21,801 | 32,579 |
| Genes | 5,143 | 10,532 |
| IDRs (>=20 aa) | 36,862 | 51,627 |
| MAIN (exact single-structure) isoforms | 15,888 | 24,803 |
| Projected segments | 51,955 | 76,781 |
| Unique exons | 80,498 | 124,613 |
| Genes with >=2 mapped transcripts | 4,530 | 8,196 |
| Alt − const exon IDR fraction, genes other than TFs (resampling genes) | 0.141 | 0.1368 |
|   95% CI low | 0.1282 | 0.126 |
|   95% CI high | 0.154 | 0.1466 |
| Adjusted OLS alt coefficient [coef, lo, hi, P] | 0.0891 / 0.0762 / 0.102 / 1.2e-41 | 0.0811 / 0.0713 / 0.0909 / 6.81e-59 |
| Gene fixed-effects within-gene difference [coef, P] | 0.1849 / 1.08e-254 | 0.1749 / 0 |
| Fixed-effects n [genes, exons] | 4,038 / 73,587 | 7,341 / 115,749 |
| Median %IDR canonical | 26.67 | 23.4 |
| Median %IDR alternative | 27.34 | 25.63 |
| Longest IDR vs %IDR, Spearman | 0.7568 / 0 | 0.7595 / 0 |
| Protein length vs %IDR, Spearman | 0.113 / 0 | 0.1279 / 0 |
| High-disorder genes (>=75%) | 525 | 1,014 |
| Very-low-disorder genes (<=15%) | 1,716 | 3,928 |
| Canonical-alternative pairs | 16,658 | 22,047 |
| Pooled disorder: spliced-out | 0.4051 | 0.3972 |
| Pooled disorder: spliced-in | 0.6287 | 0.6199 |
| Pooled disorder: shared | 0.3481 | 0.3399 |
| Canonical IDR fate: intact % | 55.9 | 56.46 |
|   internal deletion % | 21.77 | 21.16 |
|   removed % | 13.69 | 13.6 |
|   truncated % | 8.63 | 8.78 |
| Alternative IDRs: carried over % | 59.11 | 59.01 |
|   remodelled % | 40.51 | 40.59 |
|   new % | 0.38 | 0.41 |
| Pairs with an IDR status flip at shared residues | 0.708 | 0.688 |
