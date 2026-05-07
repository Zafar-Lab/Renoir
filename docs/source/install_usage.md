# Installation and Basic usage

## Requirements
All requirements are provided in the ```renoir.yml``` file. It is recommended to utilize the same versions as provided in ```renoir.yml``` file.

## Quick Start

This example walks through a complete Renoir workflow on a 10x Genomics Visium spatial transcriptomics dataset.

### Installation

```bash
conda env create -f renoir.yml
pip install .
```

> **Note:** `cell2location` must be installed separately. See [cell2location installation](https://cell2location.readthedocs.io/en/latest/).

---

### Step 1 — Imports

```python
import Renoir
import scanpy as sc
import pandas as pd
import numpy as np
```

---

### Step 2 — Required inputs

Renoir requires the following input files:

| Input | Description |
|---|---|
| `ST_path` | Spatial transcriptomics AnnData (`.h5ad`) |
| `SC_path` | Matched scRNA-seq reference AnnData with annotated cell types (`.h5ad`) |
| `pairs_path` | CSV of ligand–target pairs to score |
| `ligand_receptor_path` | Curated ligand–receptor pair table (e.g. from NATMI / OmniPath) |
| `celltype_proportions_path` | CSV of per-spot cell-type proportions (e.g. from cell2location) |
| `expins_path` | Cell-type-specific mRNA abundance pickle (see format note below) |

> **mRNA abundance file format**
>
> `mRNA_abundance.pkl` must be a dictionary where each key is a gene name and each value is the cell-type-specific mRNA abundance for that gene across all spots. Two formats are supported:
>
> **Format 1 — dictionary of DataFrames:**
> ```python
> {
>     'GENE1': pd.DataFrame(  # shape: (n_spots, n_celltypes)
>         index=['spot_1', 'spot_2', ...],       # spot barcodes
>         columns=['CellType_A', 'CellType_B', ...],
>     ),
>     'GENE2': pd.DataFrame(...),
>     ...
> }
> ```
> Each cell in the DataFrame holds the expression of that gene attributed to that cell type in that spot.
>
> **Format 2 — dictionary of CSR matrices with shared index keys:**
> ```python
> {
>     'GENE1': scipy.sparse.csr_matrix(...),  # shape: (n_spots, n_celltypes)
>     'GENE2': scipy.sparse.csr_matrix(...),
>     ...
>     'cells':     ['spot_1', 'spot_2', ...],  # row indices, shared across all matrices
>     'celltypes': ['CellType_A', 'CellType_B', ...],  # column indices, shared across all matrices
> }
> ```
> In Format 2, the `'cells'` and `'celltypes'` keys provide the row and column labels shared across all gene matrices.
>
> Cell-type-specific mRNA abundance values can be computed from [cell2location](https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html#Estimate-cell-type-specific-expression-of-every-gene-in-the-spatial-data-(needed-for-NCEM)). Currently Renoir also provides a naïve approach to computing cell-type-specific mRNA abundance by setting `expins_path=None` in the `compute_neighborhood_scores` function.

---

### Step 3 — Compute neighborhood scores

```python
neighborhood_scores = Renoir.compute_neighborhood_scores(
    SC_path='path/to/scRNA.h5ad',
    ST_path='path/to/ST.h5ad',
    pairs_path='path/to/lt_pairs.csv',
    ligand_receptor_path='path/to/All_human_lrpairs.csv',
    celltype_proportions_path='path/to/celltype_proportions.csv',
    expins_path='path/to/mRNA_abundance.pkl',
)
```

---

## References

If you use Renoir, please cite our paper:

Rao, N., Kumar, T., Kazemi, D. et al. Charting spatial ligand-target activity using Renoir. *Nature Communications* 17, 3983 (2026). [https://doi.org/10.1038/s41467-026-72388-7](https://doi.org/10.1038/s41467-026-72388-7)
