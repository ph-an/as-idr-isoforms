# Manuscript figures

Every figure in the preprint, grouped as the paper uses them. Each file is written by one script in
`analysis/manuscript_figures/` (rebuild all with `python run.py --steps paper`); the numbers on each figure come from
the outputs of the scripts in `analysis/`. PDFs are for the manuscript, PNGs for previews.

| Folder | What it holds |
|---|---|
| [`main/`](main) | Main Figures 1–7 |
| [`extended_data/`](extended_data) | Extended Data Figures 1–10 |
| [`supplementary/`](supplementary) | Exon IDR distributions by disorder class (written with the exon statistics used in Fig. 4e) |

## Main figures

| Figure | File | Script | Shows |
|---|---|---|---|
| 1 | [`main/fig1_graphical_abstract`](main/fig1_graphical_abstract.png) | `main_fig1_graphical_abstract.py` | Graphical abstract: how splicing combines ordered and disordered exons into isoforms; the approach in five steps; the five main findings |
| 2 | [`main/fig2_overview`](main/fig2_overview.png) | `main_fig2_overview.py` | Definitions (constitutive/alternative exons, shared/isoform-specific positions, IDR outcomes); exon view vs codon view; %IDR per class; range across a gene's isoforms; longest IDR vs %IDR |
| 3 | [`main/fig3_families`](main/fig3_families.png) | `main_fig3_families.py` | Most and least disordered families; annotations |
| 4 | [`main/fig4_exons`](main/fig4_exons.png) | `main_fig4_exons.py` | Alternative vs constitutive exon disorder in every class; share alternative; exon boundaries inside IDRs |
| 5 | [`main/fig5_isoforms`](main/fig5_isoforms.png) | `main_fig5_isoforms.py` | Disorder by inclusion; IDR outcomes; outcomes by class; composition of isoform-specific disorder |
| 6 | [`main/fig6_evolution`](main/fig6_evolution.png) | `main_fig6_evolution.py` | All-isoform %IDR and longest IDR through gene age per class; alternative share and alt − const by age; IDR outcomes by age |
| 7 | [`main/fig7_splice_sites`](main/fig7_splice_sites.png) | `main_fig7_splice_sites.py` | Exon disorder and alternative share by splice-site conservation; composition vs young genes |

## Extended Data

Numbered in the order the text first cites them.

| Figure | File | Script | Shows |
|---|---|---|---|
| 1 | [`extended_data/ed1_methods`](extended_data/ed1_methods.png) | `ed1_methods.py` | Methods workflow in eight steps with counts |
| 2 | [`extended_data/ed2_robustness`](extended_data/ed2_robustness.png) | `ed2_robustness.py` | Exonic vs intronic flank conservation; three results under other disorder definitions; MetaPredict vs MobiDB-lite per class; estimators of the exon difference |
| 3 | [`extended_data/ed3_disease`](extended_data/ed3_disease.png) | `ed3_disease.py` | Disease inheritance (dominant, recessive, both, none): canonical disorder, alternative − constitutive and isoform-specific − shared IDR fraction |
| 4 | [`extended_data/ed4_fate_validation`](extended_data/ed4_fate_validation.png) | `ed4_fate_validation.py` | IDR outcomes under eight comparison methods; gene weighting, isoform count, direction; disorder of kept codons |
| 5 | [`extended_data/ed5_examples`](extended_data/ed5_examples.png) | `ed5_examples.py plot PDLIM3 MDM4` | PDLIM3 and MDM4 isoforms |
| 6 | [`extended_data/ed6_composition`](extended_data/ed6_composition.png) | `ed6_composition.py` | Composition of isoform-specific vs shared disorder by inclusion (majority, half, minority) |
| 7 | [`extended_data/ed7_family_age`](extended_data/ed7_family_age.png) | `ed7_family_age.py` | Age–disorder correlation within each of 48 families |
| 8 | [`extended_data/ed8_family_isoform_age`](extended_data/ed8_family_isoform_age.png) | `ed8_family_isoform_age.py` | All-isoform %IDR and longest IDR through gene age in 15 families |
| 9 | [`extended_data/ed9_exon_age`](extended_data/ed9_exon_age.png) | `ed9_exon_age.py` | Splice-site conservation within genes: paired genes; the within-gene result under alternative explanations |
| 10 | [`extended_data/ed10_exon_age_families`](extended_data/ed10_exon_age_families.png) | `ed10_exon_age_families.py` | Exon IDR by splice-site conservation across 30 families |

## Supplementary

| File | Script | Shows |
|---|---|---|
| [`supplementary/s2_exons_by_disorder_class`](supplementary/s2_exons_by_disorder_class.png) | `supp2_exons_by_disorder_class.py` | Exon IDR distributions and the alternative–constitutive difference by disorder class |
