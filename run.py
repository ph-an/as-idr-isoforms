"""Rebuild every table and figure from the raw inputs, in dependency order.

    python run.py                          # everything (~40 min on a laptop)
    python run.py --steps analysis paper   # manuscript analyses and figures only (needs data/ built)
    python run.py --list                   # show the steps

Steps
  inputs    check the SHA-256 of every input listed in data/manifest.tsv
  proteome  run notebooks/01-03 top to bottom in fresh kernels -> data/annotated/, figures/proteome/
  genomic   pipeline/run_pipeline.py (UniProt->Ensembl mapping, GENCODE coordinates, IDR projection,
            exon tables, splicing analyses) -> data/genomic/, figures/genomic/
  analysis  the manuscript analyses in analysis/ (exon usage, disorder landscape, IDR architecture,
            composition, genomic origin, disorder classes) -> analysis/*/outputs/ (every number in the paper)
  paper     manuscript figures -> manuscript/figures/{main,extended_data,supplementary}/

Logs go to results/logs/ and executed notebooks to results/executed_notebooks/ (neither tracked by git).
A complete run writes results/run_manifest.json (environment, parameters, timings, output checksums);
a partial run writes the same record to results/logs/partial_run_manifest.json.
PDFs carry a fixed creation date (SOURCE_DATE_EPOCH) so an unchanged figure is byte-identical after a rerun.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import importlib.metadata as md
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent
RESULTS = REPO / "results"
LOGS = RESULTS / "logs"
MANIFEST = REPO / "data" / "manifest.tsv"

NOTEBOOKS = [
    "notebooks/01_cohort_and_idrs",               # cohort filters, IDR segments, TF labels -> isoform master tables
    "notebooks/02_proteome_analysis",             # gene age; figures/proteome (fig1-fig7)
    "notebooks/03_appris_annotation",             # APPRIS gene annotation (read by the genomic pipeline)
]
# manuscript analyses, in dependency order (script path, extra arguments)
ANALYSIS_SCRIPTS = [
    ("analysis/reanalysis_2026-10-26/reanalysis.py", []),            # exon usage x disorder, gene-aware (R1-R6)
    ("analysis/reanalysis_2026-10-26/cds_level_usage.py", []),       # usage on coding coordinates (sensitivity)
    ("analysis/idr_landscape/idr_landscape.py", []),                 # disorder landscape, disorder classes
    ("analysis/idr_architecture/idr_architecture.py", []),           # canonical-alternative alignment, IDR fates
    ("analysis/idr_architecture/composition_from_chunks.py", []),    # composition + non-zinc-finger sensitivity
    ("analysis/idr_architecture/spliced_sequence_origin.py", []),    # genomic origin of alternative-only disorder
    ("analysis/idr_landscape/splicing_by_disorder_class.py", []),    # coupling by disorder class
    ("analysis/reference_free/reference_free.py", []),              # no reference isoform: core vs variable, all-pairs fates
    ("analysis/reference_free/codon_vs_exon.py", []),               # exon labels vs codon labels (Fig. 2b)
    ("analysis/reference_free/fate_method_check.py", []),           # alignment vs genomic-codon IDR fates
    ("analysis/protein_classes/gene_annotations.py", []),           # UniProt families, classes, features, gene age per gene
    ("analysis/protein_classes/protein_classes.py", []),            # disorder and splicing by protein class and family
    ("analysis/evolution/evolution.py", []),                        # gene age: disorder, splicing, IDR fates, composition, inheritance
    ("analysis/exon_age/exon_age.py", []),                          # exon age from splice-site conservation (405 mammals)
    ("analysis/exon_age/exon_age_robustness.py", []),               # exon age: alignability controls, other definitions, matched pairs
    ("analysis/reference_free/fate_robustness.py", []),             # IDR fates: alignment methods, junction codons, weighting, disorder kept
    ("analysis/sensitivity/predict_scores.py", []),                 # per-residue MetaPredict v3 scores (slow; restartable)
    ("analysis/sensitivity/predict_v1_subset.py", []),              # MetaPredict V1 (not trained on pLDDT), young-exon genes + matched genes
    ("analysis/sensitivity/predict_v1_subset.py", ["all"]),         # ... and the rest of the cohort
    ("analysis/sensitivity/sensitivity.py", []),                    # main results under other disorder calls
]
PAPER_SCRIPTS = [
    # manuscript/figures/main/
    ("analysis/manuscript_figures/cohort_counts.py", []),
    ("analysis/manuscript_figures/supp2_exons_by_disorder_class.py", []),   # also writes the exon statistics used by Fig. 4e
    ("analysis/manuscript_figures/main_fig1_graphical_abstract.py", []),
    ("analysis/manuscript_figures/main_fig2_overview.py", []),
    ("analysis/manuscript_figures/main_fig3_families.py", []),
    ("analysis/manuscript_figures/main_fig4_exons.py", []),
    ("analysis/manuscript_figures/main_fig5_isoforms.py", []),
    ("analysis/manuscript_figures/main_fig6_evolution.py", []),
    # manuscript/figures/extended_data/
    ("analysis/manuscript_figures/ed1_methods.py", []),
    ("analysis/manuscript_figures/ed6_composition.py", []),
    ("analysis/manuscript_figures/ed7_family_age.py", []),
    ("analysis/manuscript_figures/ed9_exon_age.py", []),
    ("analysis/manuscript_figures/ed10_exon_age_families.py", []),
    ("analysis/manuscript_figures/ed5_examples.py", ["candidates"]),
    ("analysis/manuscript_figures/ed5_examples.py", ["plot", "PDLIM3", "MDM4"]),
    ("analysis/manuscript_figures/ed4_fate_validation.py", []),
    ("analysis/manuscript_figures/ed2_robustness.py", []),
    ("analysis/manuscript_figures/ed3_disease.py", []),
    ("analysis/evolution/isoform_disorder_by_age.py", []),
    ("analysis/manuscript_figures/main_fig7_splice_sites.py", []),
    ("analysis/manuscript_figures/ed8_family_isoform_age.py", []),
]


OUTPUT_DIRS = ["data/annotated", "data/genomic/analysis", "data/genomic/processed", "data/genomic/isoform_mapping",
               "data/genomic/interim", "data/genomic/qc", "figures", "manuscript/figures"]
PACKAGES = ["numpy", "pandas", "scipy", "matplotlib", "seaborn", "statsmodels", "pyarrow", "openpyxl",
            "biopython", "tqdm", "requests", "nbclient", "ipykernel"]
ENV = {"PYTHONUTF8": "1", "PYTHONHASHSEED": "0", "MPLBACKEND": "Agg", "SOURCE_DATE_EPOCH": "0"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()                  # chunked, so Python 3.10 works too
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def installed_version(package: str) -> str:
    try:
        return md.version(package)
    except md.PackageNotFoundError:
        return "not installed"


def run_cmd(name: str, cmd: list[str], cwd: Path) -> None:
    print(f"  -> {name}", flush=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    with (LOGS / f"{name}.log").open("w", encoding="utf-8") as log:
        rc = subprocess.run(cmd, cwd=cwd, env=os.environ, stdout=log, stderr=subprocess.STDOUT).returncode
    if rc:
        raise SystemExit(f"{name} failed (exit {rc}); see {LOGS / (name + '.log')}")


def step_inputs() -> dict:
    rows = list(csv.DictReader(MANIFEST.open(encoding="utf-8"), delimiter="\t"))
    bad = []
    for r in rows:
        p = REPO / r["path"]
        if not p.exists():
            bad.append(f"MISSING  {r['path']}")
        elif sha256(p) != r["sha256"]:
            bad.append(f"CHANGED  {r['path']}")
    if bad:
        raise SystemExit("Input check failed (see data/README.md for how to obtain the inputs):\n  " + "\n  ".join(bad))
    print(f"  {len(rows)} inputs verified", flush=True)
    return {"verified_inputs": len(rows)}


def step_proteome() -> dict:
    import nbformat
    from nbclient import NotebookClient

    out_dir = RESULTS / "executed_notebooks"
    out_dir.mkdir(parents=True, exist_ok=True)
    timings = {}
    for rel in NOTEBOOKS:
        path = REPO / f"{rel}.ipynb"
        print(f"  -> {rel}.ipynb", flush=True)
        nb = nbformat.read(path, as_version=4)
        start = time.time()
        client = NotebookClient(nb, timeout=None, kernel_name="python3", resources={"metadata": {"path": str(path.parent)}})
        try:
            client.execute()
        finally:
            nbformat.write(nb, out_dir / f"{path.stem}.ipynb")
        timings[path.stem] = round(time.time() - start, 1)
    return {"notebook_seconds": timings}


def step_genomic() -> dict:
    start = time.time()
    run_cmd("genomic_pipeline", [sys.executable, "run_pipeline.py"], REPO / "pipeline")
    return {"genomic_seconds": round(time.time() - start, 1)}


def _run_scripts(prefix: str, scripts: list) -> dict:
    timings = {}
    for rel, extra in scripts:
        start = time.time(); path = REPO / rel
        name = f"{prefix}_{path.parent.name}_{path.stem}" + ("_" + "_".join(extra) if extra else "")
        run_cmd(name, [sys.executable, str(path), *extra], path.parent)
        timings[name] = round(time.time() - start, 1)
    return timings


def step_analysis() -> dict:
    return {"analysis_seconds": _run_scripts("analysis", ANALYSIS_SCRIPTS)}


def step_paper() -> dict:
    return {"paper_seconds": _run_scripts("paper", PAPER_SCRIPTS)}


STEPS = {"inputs": step_inputs, "proteome": step_proteome, "genomic": step_genomic, "analysis": step_analysis, "paper": step_paper}


def output_checksums() -> dict:
    sums = {}
    for d in OUTPUT_DIRS:
        for p in sorted((REPO / d).rglob("*")):
            if p.is_file():
                sums[p.relative_to(REPO).as_posix()] = sha256(p)
    return sums


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--steps", nargs="+", choices=list(STEPS), default=list(STEPS))
    parser.add_argument("--list", action="store_true", help="print the step descriptions and exit")
    args = parser.parse_args()
    if args.list:
        print(__doc__)
        return

    os.environ.update(ENV)   # inherited by notebook kernels and subprocesses
    min_idr = os.environ.get("TFIDR_MIN_IDR_LEN_AA", "20")
    record = {
        "started": dt.datetime.now().isoformat(timespec="seconds"),
        "steps": args.steps,
        "parameters": {"TFIDR_MIN_IDR_LEN_AA": min_idr, "TFIDR_MIN_ISOFORMS": os.environ.get("TFIDR_MIN_ISOFORMS", "2"),
                       "TFIDR_ANNOTATIONS": os.environ.get("TFIDR_ANNOTATIONS", "local")},
        "python": sys.version,
        "platform": platform.platform(),
        "packages": {p: installed_version(p) for p in PACKAGES},
    }
    if min_idr != "20":
        print(f"WARNING: TFIDR_MIN_IDR_LEN_AA={min_idr} (not the paper's IDR definition)", flush=True)
    for name in args.steps:
        print(f"[{name}]", flush=True)
        record[name] = STEPS[name]()
    record["finished"] = dt.datetime.now().isoformat(timespec="seconds")
    record["outputs_sha256"] = output_checksums()
    # the tracked manifest describes a complete rebuild; partial runs are logged next to their step logs
    target = RESULTS / "run_manifest.json" if args.steps == list(STEPS) else LOGS / "partial_run_manifest.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
    print(f"done -> {target}")


if __name__ == "__main__":
    main()
