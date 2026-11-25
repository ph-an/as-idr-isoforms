# Reference-free IDR architecture

The canonical-referenced analyses (`analysis/idr_architecture`) compare each alternative isoform with the UniProt
canonical isoform by protein-sequence alignment. Here no isoform is a reference: every coding residue of every mapped
isoform is keyed by the genomic position of the first base of its codon, so a gene's isoforms share one coordinate
system and a residue is "shared" only if the same codon encodes it.

* `reference_free.py` — core (in every isoform) vs variable positions; disorder by inclusion level; composition of
  disordered positions by inclusion and position vs core disorder; IDR fates over all ordered isoform pairs; disorder
  classes from core %IDR and median isoform %IDR. Outputs in `outputs/` (`key_numbers.json`).
* `fate_method_check.py` — the same canonical IDRs classified by alignment and by shared codons, and genomic-codon fates
  in both directions (`outputs/fate_method_check.json`).

Main points (2026-11-04):
* Variable positions are more disordered than core positions (41% vs 32% pooled; gene-paired P = 1e-129), most of all
  positions present in a minority of isoforms (51%).
* IDR fates depend on how shared residues are defined, not on the reference: alignment calls 22% of canonical IDRs
  internally deleted, shared codons call 9.5%. 57% of alignment "internal deletions" are truncations or full removals
  at the codon level (coincidental sequence matches inside IDRs). With codons, partial loss (internal deletion +
  truncation, 25%) still exceeds whole removal (16%), and internal deletion still dominates for IDRs > 200 aa.
* The composition signature (acidic down, histidine/cysteine/aromatic up, net charge up) is concentrated in disorder
  present in a minority of isoforms and at termini; internal variable disorder shows mainly aliphatic depletion and
  proline/histidine enrichment.
