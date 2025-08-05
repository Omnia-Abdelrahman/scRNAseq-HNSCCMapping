# HNSCC Single-Cell Atlas: Decoding Tumor Microenvironments One Cell at a Time

Welcome to the **HNSCC Single-Cell Mapping Project** project, which uses **single-cell RNA sequencing (scRNA-seq)** data to dive deep into the cellular world of **head and neck squamous cell carcinoma (HNSCC)** and to understand what makes tumor and non-tumor cells different.

This isn’t just a casual look at data. It's a step-by-step journey toward identifying targets that could eventually help develop better treatment/preventative decisions. Everything here is done openly and simply so others can learn, build, or collaborate.

---

## Project Overview

This project uses real scRNA-seq data to:

- Separate cancer cells from non-cancer cells
- Compare their gene activity
- Visualize how they cluster and differ
- Explore possible genes that could be targeted in therapy
- Recreate and build on the open dataset shared by **Puram et al. (2017)**

---

## Data Source and Credit

This work is possible thanks to the open data from:

**Puram et al., Cell (2017):**

> "Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor Ecosystems in Head and Neck Cancer"

You can find the dataset here: [**GEO Accession GSE103322**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322)

---

## Project Folder Structure

```
HNSCC-scRNAseq-Mapping/
├── scRNAseq-HNSCCMapping.py       # The main code for running all steps (so far, from Phase 1)
├── data/                          # Contains input files (ex. matrices)
├── figures/                       # For saving plots like UMAPs and violins
├── outputs/                       # Processed AnnData (.h5ad) and DE results
└── README.md                      # You’re reading it
```

---

## Step-by-Step Plan

### Phase 1: Setup and Visualization (COMPLETED)

- Load the scRNA-seq data file (annotated matrix)
- Separate expression from metadata
- Add labels to say which cells are cancer and which are not
- Normalize the gene counts
- Log-transform and reduce data with UMAP
- Plot gene expression (e.g. EPCAM, SPP1) across cells

### Phase 2: Find Key Genes

- Run statistical tests to find genes that differ between cancer and non-cancer cells
- Plot top genes using dot plots and heatmaps
- Create gene lists for further research

### Phase 3: Link to Biology

- Search for genes related to immunity or stemness
- Use public databases to learn more about those genes
- Suggest possible vaccine targets based on results

### Phase 4: Share the Project (Future)

- Build a basic app with Streamlit or Dash
- Let users search and visualize by gene or cell type
- Publish the results for everyone to explore

---

## Why This Project?

This project started from a simple goal: **learn and build something real with single-cell data**. The HNSCC tumor microenvironment is complex. By breaking it down one cell at a time, we aim to:

- Understand the tumor better
- Practice reproducible bioinformatics
- Open the door to hypothesis generation and future lab work

---

## Environment and Requirements

### Python Version

```
Python 3.10 
```

### Key Libraries

Install with pip:

```bash
pip install scanpy pandas numpy matplotlib seaborn
```

Or with conda:

```bash
conda install -c bioconda scanpy
```

---

## Example Outputs

- UMAPs colored by cancer status
- Violin plots showing marker expression
- Lists of genes differentially expressed in cancer
- Clean `.h5ad` file to reuse later

---

## How to Follow Along

This is an open, evolving project. New results and plots will be shared often. You can:

- Watch the repo for updates
- Fork and try it yourself
- Suggest new genes or analysis directions
- Help annotate or interpret findings

---

## Legal Note

This is a research and learning project built on **public data** from Puram et al. (2017). No patient-identifiable information is used. This work is legal under data reuse terms and is intended for academic and educational purposes only. No license included yet.

---

## Citation

If this helps your work or learning, please cite:

```
Puram SV, Tirosh I, Parikh AS, et al. Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor Ecosystems in Head and Neck Cancer. Cell. 2017.
```

---

## Let’s Build Together

If you’re into science, coding, or bioinformatics—and especially if you're learning, you’re welcome here.

Questions? Ideas? Feedback? Always appreciated.


