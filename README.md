# Neuro Data Analysis

A computational neuroscience portfolio built during my research internship
at IISER Tirupati (Dr. Devanathan's Lab, 2026).

This repository contains Python-based analysis pipelines for real
neurological disease datasets, built from scratch as part of developing
skills at the intersection of biology and data science.

---

## Notebooks & Scripts

### 1. Gene Expression Analysis (`gene_expression_analysis.py`)
Simulated qRT-PCR data analysis across three conditions — Control,
RA-Differentiated, and Treated — for key neurological genes (CHAT, TH, BDNF).

**Skills demonstrated:**
- Biological replicate simulation with NumPy
- Statistical testing (t-test, p-values, significance stars)
- Publication-quality bar charts and heatmaps with seaborn/matplotlib
- Fold change visualization with baseline reference

---

### 2. Alzheimer's Disease — Real GEO Dataset (`alzheimers_analysis.py`)
Analysis of GSE5281 — a published microarray dataset from human
Alzheimer's brain tissue (161 samples across multiple brain regions).

**What this analysis finds:**
- 2,284 genes downregulated and 781 upregulated in Alzheimer's vs control
- Top hits include SNAP25, SYT1, VAMP2, NRXN1, NRXN3, GRIN2B — the core
  synaptic release machinery, consistent with the synaptic failure
  hypothesis of Alzheimer's disease
- Volcano plot showing genome-wide significance vs fold change

**Skills demonstrated:**
- Loading and parsing real NCBI GEO series matrix files
- Vectorized t-tests across 54,000+ probes
- Volcano plot generation
- Probe-to-gene symbol mapping via mygene.info API

---

### 3. Comparative Disease Analysis (`comparative_analysis.py`)
A master pipeline that loads any two GEO neurological disease datasets
and compares their gene expression signatures side by side.

Currently comparing:
- **Alzheimer's Disease** — GSE5281 (human brain, 161 samples)
- **Parkinson's Disease** — GSE7621 (substantia nigra, 25 samples)

**Output generated:**
- Volcano plot per disease with gene names annotated on top hits
- Bar chart of top 20 upregulated and downregulated genes per disease
- Combined volcano showing both disease signatures on one plot
- Shared gene analysis — genes dysregulated in both diseases
- Opposite direction analysis — genes up in one disease, down in the other
- Full results exported to CSV

**Key biological question this addresses:**
Which genes fail across multiple neurodegenerative diseases?
Genes consistently dysregulated in both AD and PD are strong candidates
for shared neurodegeneration mechanisms and potential therapeutic targets.

**Skills demonstrated:**
- Multi-dataset comparative transcriptomics
- GEO platform annotation file parsing (GPL NCBI download)
- Probe collapsing by variance (one gene = one row)
- FDR correction via Benjamini-Hochberg method
- Disk caching of API results for reproducibility
- Modular pipeline design — works with any two GEO .txt.gz files

---

## Tech Stack

Python | Pandas | NumPy | SciPy | Statsmodels
Matplotlib | Seaborn | Requests | mygene.info API | NCBI GEO

---

## Background

BSc Life Sciences student at UM-DAE CEBS, Mumbai.
Currently a research intern at IISER Tirupati in Dr. Devanathan's
neurobiology lab, which focuses on synaptic signaling, cholinergic and
dopaminergic systems, and neurodegeneration.

This portfolio was built independently alongside the internship to develop
computational skills applicable to molecular neuroscience research.
