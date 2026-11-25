"""Cohort and pipeline counts for Figure 1 and the Results text (reads the current pipeline outputs)."""
import json, os
from pathlib import Path
import pandas as pd
R = Path(os.environ.get('TFIDR_ROOT', Path(__file__).resolve().parents[2]))
iso = pd.read_csv(R / 'data/annotated/IDRisoforms_df.csv', usecols=['isoform_accession', 'base_accession', 'is_canonical'])
idr = pd.read_csv(R / 'data/annotated/idr_df.csv', usecols=['isoform_accession', 'idr_aaLen'])
m = pd.read_parquet(R / 'data/genomic/processed/tf_idr_isoform_master.parquet', columns=['isoform_accession', 'genomic_scaffold_validated', 'n_candidate_ENST', 'sequence_match_type'])
res = pd.read_parquet(R / 'data/genomic/isoform_mapping/isoform_transcript_mapping_resolved.parquet')
seg = pd.read_parquet(R / 'data/genomic/interim/genomic_projection_segments.parquet')
E = pd.read_parquet(R / 'data/genomic/processed/exon_level_table.parquet', columns=['gene_id', 'exon_id', 'is_coding_exon', 'exon_usage', 'isoform_accession'])
ec = E[E.is_coding_exon]
main = res[res.resolved_mapping_status.isin(['resolved_unique', 'resolved_collapsed_equivalent'])]
N = dict(min_isoforms=int(os.environ.get('TFIDR_MIN_ISOFORMS', '2')), isoforms=len(iso), genes=int(iso.base_accession.nunique()),
         canonical=int(iso.is_canonical.sum()), alternative=int((~iso.is_canonical).sum()),
         idrs=int((idr.idr_aaLen >= 20).sum()), with_candidate_enst=int((m.n_candidate_ENST > 0).sum()),
         exact_match=int((m.sequence_match_type == 'exact').sum()), main=int(len(main)),
         ambiguous=int(res.resolved_mapping_status.str.contains('ambiguous').sum()),
         scaffold_validated=int(m.genomic_scaffold_validated.sum()), segments=int(len(seg)),
         coding_exon_rows=int(len(ec)), mapped_genes=int(ec.gene_id.nunique()),
         unique_exons=int(ec[ec.exon_usage.isin(['alternative', 'constitutive'])].drop_duplicates(['gene_id', 'exon_id']).shape[0]),
         genes_two_plus_mapped=int(ec[ec.exon_usage.isin(['alternative', 'constitutive'])].gene_id.nunique()))
# the exon analyses (reanalysis_2026-10-26) key exons by (UniProt gene accession, exon ID) over the s9 annotation table;
# report those counts too so Figure 1 matches the Results
ux = Path(__file__).resolve().parents[1] / 'reanalysis_2026-10-26' / 'outputs' / 'unique_exons.csv'
if ux.exists():
    u = pd.read_csv(ux, usecols=['gene', 'alt']); N.update(analysis_unique_exons=len(u), analysis_genes=int(u.gene.nunique()),
                                                      analysis_alt=int(u.alt.sum()), analysis_const=int((1 - u.alt).sum()))
out = Path(__file__).resolve().parent / 'outputs' / 'cohort_counts.json'; json.dump(N, open(out, 'w'), indent=1); print(json.dumps(N, indent=1))
