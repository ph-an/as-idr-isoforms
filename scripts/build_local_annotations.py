"""Build TF labels and gene ages for the whole reviewed human proteome from local files only.

Replaces the cluster-only Ensembl-gene -> UniProt table (data/external/hpc_snapshot/README.md) with UniProt's own
Ensembl cross-references (data/raw/uniprot/uniprot_ensembl_xref_raw.tsv, the file the genomic pipeline already uses).
Because the snapshot covered only genes with >=3 isoforms, any other cohort definition needs this.

Rules
  * ENSG -> UniProt: every unversioned ENSG listed in an entry's Ensembl cross-references (primary assembly).
  * TF labels: HumanTFs v1.01 rows with "Is TF? == Yes". When one ENSG maps to several entries, keep the entries whose
    UniProt primary gene name equals the HumanTFs HGNC symbol; Q13948 (CASP, a non-DNA-binding product of the CUX1
    locus) is excluded explicitly so CUX1 maps to P39880, as on the cluster.
  * Gene age: GenOrigin age per ENSG (">4290" -> 4300). An entry with several ENSGs keeps the age of the first ENSG in
    UniProt's cross-reference order that GenOrigin dates.

Check against the cluster snapshot (>=3-isoform cohort): TF labels identical (416 genes, same DBD); gene ages differ
for 15 of 5,108 accessions (0.3%), all entries with several ENSG IDs where the two tables pick a different one.

Outputs (tracked): data/external/local_mapping/tf_labels_uniprot_xref.csv, gene_age_uniprot_xref.csv
"""
import re
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
RAW = REPO / "data" / "raw"
OUT = REPO / "data" / "external" / "local_mapping"
OUT.mkdir(parents=True, exist_ok=True)
EXCLUDE_TF_ENTRIES = {"Q13948"}

x = pd.read_csv(RAW / "uniprot" / "uniprot_ensembl_xref_raw.tsv", sep="\t")
pairs = [(g, e) for e, s in zip(x["Entry"], x["Ensembl"].fillna("")) for g in dict.fromkeys(re.findall(r"(ENSG\d+)", s))]
m = pd.DataFrame(pairs, columns=["ENSG", "base_accession"])
names = pd.read_csv(RAW / "uniprot" / "uniprotkb_annotated", sep="\t", usecols=["Entry", "Gene Names (primary)"],
                    low_memory=False).set_index("Entry")["Gene Names (primary)"]

tf = pd.read_excel(RAW / "humantfs" / "DatabaseExtract_v_1.01.xlsx")
tf = tf[tf["Is TF?"] == "Yes"][["Ensembl ID", "HGNC symbol", "DBD"]]
t = m.merge(tf, left_on="ENSG", right_on="Ensembl ID")
multi = t.groupby("ENSG")["base_accession"].transform("nunique") > 1
t = t[~multi | (t["base_accession"].map(names) == t["HGNC symbol"])]
t = t[~t["base_accession"].isin(EXCLUDE_TF_ENTRIES)]
tf_labels = t[["base_accession", "DBD"]].drop_duplicates("base_accession").sort_values("base_accession")
tf_labels.to_csv(OUT / "tf_labels_uniprot_xref.csv", index=False)

go = pd.read_csv(RAW / "genorigin" / "Homo_sapiens.csv")
a = m.merge(go, left_on="ENSG", right_on="ensembl_gene_id")
a.loc[a["gene_age"] == ">4290", "gene_age"] = "4300"
a["gene_age"] = a["gene_age"].astype(int)
a = a.drop_duplicates("base_accession")  # first ENSG in UniProt cross-reference order
a["age_bin"] = pd.cut(a["gene_age"], [0, 100, 500, 1000, np.inf], labels=["<100", "100-500", "500-1000", ">1000"], right=False)
age = a.rename(columns={"base_accession": "upkb_accession"})[["upkb_accession", "gene_age", "age_bin"]].sort_values("upkb_accession")
age.to_csv(OUT / "gene_age_uniprot_xref.csv", index=False)

snap = REPO / "data" / "external" / "hpc_snapshot"
s_tf = pd.read_csv(snap / "tf_labels_from_hpc_mapping.csv")
s_age = pd.read_csv(snap / "gene_age_from_hpc_mapping.csv")
cmp_tf = s_tf.merge(tf_labels, on="base_accession", how="left", suffixes=("_hpc", "_local"))
j = s_age.merge(age, on="upkb_accession", how="left", suffixes=("_hpc", "_local"))
print(f"TF labels: {len(tf_labels)} proteome-wide; snapshot genes recovered {cmp_tf.DBD_local.notna().sum()}/{len(s_tf)}, "
      f"DBD identical {(cmp_tf.DBD_hpc == cmp_tf.DBD_local).sum()}")
print(f"Gene ages: {len(age)} proteome-wide; snapshot accessions differing {(j.gene_age_hpc != j.gene_age_local).sum()}/{len(s_age)}")
