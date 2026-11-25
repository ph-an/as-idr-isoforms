"""Gene-level protein classifications from the UniProt annotation table, joined to disorder and splicing outcomes.

One row per gene (UniProt base accession) in the >=2-isoform cohort. Writes outputs/gene_annotations.csv.

Classifications (from data/raw/uniprot/uniprotkb_annotated, canonical entry)
  family            first UniProt 'Protein families' entry (e.g. 'Krueppel C2H2-type zinc-finger protein family')
  enzyme, ec_class  EC number present; first EC digit (oxidoreductase ... translocase)
  n_tm, transmembrane, signal_peptide
  loc_*             subcellular location from the curated SUBCELLULAR LOCATION comment (multi-label)
  kw_*              molecular-function / process keywords (UniProt keywords)
  disease           'Involvement in disease' present
  n_modres, n_phospho  curated modified residues (all; phosphorylation) per protein
  n_interactors     curated binary interactors ('Interacts with')
  has_pdb, disprot, condensate (CD-CODE entry), coiled_coil, n_pfam, n_domains
  mobidb_pct        % of residues in UniProt 'Disordered' regions (MobiDB-lite), an independent disorder predictor

Outcomes
  pct_idr, length, n_iso, idr_range           canonical %IDR (MetaPredict), length, isoforms, %IDR range (idr_landscape)
  median_iso_pct_idr, dis_core, dis_var        reference-free (analysis/reference_free)
  gene_age, age_bin, is_tf                     GenOrigin age (Ma), TF label
"""
import os, re
import numpy as np, pandas as pd

ROOT = os.environ.get('TFIDR_ROOT', os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
HERE = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(HERE, '..'); OUT = os.path.join(HERE, 'outputs'); os.makedirs(OUT, exist_ok=True)
COLS = ['Entry', 'Length', 'Protein families', 'EC number', 'Subcellular location [CC]', 'Transmembrane', 'Signal peptide',
        'Involvement in disease', 'Modified residue', 'Interacts with', 'PDB', 'DisProt', 'CD-CODE', 'Coiled coil', 'Pfam',
        'Keywords', 'Domain [FT]', 'Region']
u = pd.read_csv(os.path.join(ROOT, 'data/raw/uniprot/uniprotkb_annotated'), sep='\t', usecols=COLS, low_memory=False).set_index('Entry')
g = pd.read_csv(os.path.join(A, 'idr_landscape/outputs/gene_disorder_class.csv'), index_col=0)
u = u.reindex(g.index)
s = lambda c: u[c].fillna('').astype(str)
D = pd.DataFrame(index=g.index)
D['family'] = s('Protein families').str.split(';').str[0].str.strip().replace('', np.nan)
D['family'] = D.family.str.split(',').str[0].str.strip()          # 'X family, Y subfamily' -> 'X family'
EC = {'1': 'oxidoreductase', '2': 'transferase', '3': 'hydrolase', '4': 'lyase', '5': 'isomerase', '6': 'ligase', '7': 'translocase'}
D['enzyme'] = s('EC number').ne('')
D['ec_class'] = s('EC number').str[0].map(EC)
D['n_tm'] = s('Transmembrane').str.count('TRANSMEM'); D['transmembrane'] = D.n_tm > 0
D['signal_peptide'] = s('Signal peptide').str.contains('SIGNAL')
loc = s('Subcellular location [CC]').str.lower().str.replace(r'\{[^}]*\}', '', regex=True)
LOCS = {'nucleus': r'nucle(us|olus|oplasm)|chromosome', 'cytoplasm': r'cytoplasm|cytosol', 'cell_membrane': r'cell membrane',
        'secreted': r'secreted|extracellular', 'mitochondrion': r'mitochondri', 'er': r'endoplasmic reticulum', 'golgi': r'golgi',
        'cytoskeleton': r'cytoskeleton|microtubule|centrosome|cilium', 'synapse': r'synapse|postsynaptic|presynaptic', 'vesicle': r'vesicle|endosome|lysosome'}
for k, pat in LOCS.items():
    D[f'loc_{k}'] = loc.str.contains(pat.replace('(', '(?:')).astype(float)
D.loc[loc.eq(''), [f'loc_{k}' for k in LOCS]] = np.nan
kw = s('Keywords')
KW = {'transcription': 'Transcription regulation', 'dna_binding': 'DNA-binding', 'rna_binding': 'RNA-binding', 'kinase': 'Kinase',
      'phosphatase': 'Protein phosphatase', 'receptor': 'Receptor', 'gpcr': 'G-protein coupled receptor', 'ion_channel': 'Ion channel',
      'transport': 'Transport', 'ubl_conjugation': 'Ubl conjugation pathway', 'chromatin_regulator': 'Chromatin regulator',
      'mrna_processing': 'mRNA processing', 'developmental': 'Developmental protein', 'cell_cycle': 'Cell cycle',
      'cytoskeleton_kw': 'Cytoskeleton', 'metal_binding': 'Metal-binding', 'zinc_finger': 'Zinc-finger', 'glycoprotein': 'Glycoprotein',
      'immunity': 'Immunity', 'signalling_adaptor': 'SH3 domain', 'neurogenesis': 'Neurogenesis'}
for k, w in KW.items():
    D[f'kw_{k}'] = kw.str.split(';').apply(lambda L, w=w: w in L)
D['disease'] = s('Involvement in disease').str.contains('DISEASE')
D['n_modres'] = s('Modified residue').str.count('MOD_RES')
D['n_phospho'] = s('Modified residue').str.count(r'/note="Phospho')
D['n_interactors'] = s('Interacts with').apply(lambda x: len([y for y in x.split(';') if y.strip()]))
D['has_pdb'] = s('PDB').ne(''); D['disprot'] = s('DisProt').ne(''); D['condensate'] = s('CD-CODE').ne('')
D['coiled_coil'] = s('Coiled coil').str.contains('COILED')
D['n_pfam'] = s('Pfam').apply(lambda x: len([y for y in x.split(';') if y.strip()]))
D['n_domains'] = s('Domain [FT]').str.count('DOMAIN ')
def mobidb(txt, L):
    tot = 0
    for m in re.finditer(r'REGION (\d+)\.\.(\d+); /note="Disordered"', txt):
        tot += int(m.group(2)) - int(m.group(1)) + 1
    return 100 * tot / L if L else np.nan
D['mobidb_pct'] = [mobidb(t, L) for t, L in zip(s('Region'), u.Length.fillna(0))]
# outcomes
D['pct_idr'] = g.pct_idr; D['length'] = g.Length; D['n_iso'] = g.n_iso; D['idr_range'] = g.range; D['is_tf'] = g.is_tf
D['dclass'] = g.dclass
rf = pd.read_csv(os.path.join(A, 'reference_free/outputs/gene_core_variable.csv')).set_index('gene')
for c in ['median_iso_pct_idr', 'dis_core', 'dis_var', 'n_core', 'n_var', 'core_pct_idr']:
    D[c] = rf[c]
age = pd.read_csv(os.path.join(ROOT, 'data/annotated/IDRisoforms_df_geneage.csv'), usecols=['base_accession', 'gene_age', 'age_bin', 'gene_name'],
                  low_memory=False).drop_duplicates('base_accession').set_index('base_accession')
D['gene_age'] = age.gene_age; D['age_bin'] = age.age_bin; D['gene_name'] = age.gene_name
D.index.name = 'gene'
D.to_csv(os.path.join(OUT, 'gene_annotations.csv'))
print(D.shape); print(D.drop(columns=['family', 'ec_class', 'dclass', 'age_bin', 'gene_name']).mean(numeric_only=True).round(3).to_string())
print('families:', D.family.nunique(), 'genes with family:', D.family.notna().sum())

# ---- broad protein class (mutually exclusive; first matching rule in this priority order wins)
kw = lambda c: D[c].fillna(False).astype(bool)
RULES = [('DNA-binding TFs', kw('is_tf')),
         ('Transcription co-regulators', kw('kw_transcription') | kw('kw_dna_binding')),
         ('Chromatin & RNA processing', kw('kw_chromatin_regulator') | kw('kw_mrna_processing') | kw('kw_rna_binding')),
         ('Kinases & phosphatases', kw('kw_kinase') | kw('kw_phosphatase')),
         ('Receptors & channels', kw('kw_receptor') | kw('kw_gpcr') | kw('kw_ion_channel')),
         ('Transporters', kw('kw_transport') & kw('transmembrane')),
         ('Metabolic enzymes', kw('enzyme')),
         ('Cytoskeletal & structural', kw('kw_cytoskeleton_kw') | D.loc_cytoskeleton.fillna(0).astype(bool)),
         ('Secreted', kw('signal_peptide') & ~kw('transmembrane')),
         ('Other membrane', kw('transmembrane'))]
pclass = pd.Series('Other / unclassified', index=D.index)
for name, m in reversed(RULES):
    pclass[m] = name
D['pclass'] = pclass
D.to_csv(os.path.join(OUT, 'gene_annotations.csv'))
print(D.groupby('pclass').agg(genes=('pct_idr', 'size'), median_pct_idr=('pct_idr', 'median'), mean_iso=('n_iso', 'mean'),
                              median_range=('idr_range', 'median'), frac_tf=('is_tf', 'mean')).sort_values('median_pct_idr').round(2).to_string())
