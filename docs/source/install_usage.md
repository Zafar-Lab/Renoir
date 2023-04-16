# Installation and Basic usage

## Requirements
All requirements are provided in the renoir.yml file. It is recommended to utilize the same versions as provided in renoir.yml.

## Installation
**Note**: Installation does not include cell2location. Please install cell2location separately.

- Download the repository and install the conda environment using `conda env create -f renoir.yml`
- Install Renoir by cd to this directory and running `pip install .`

## Usage

```
import Renoir
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
import pickle
```

**Read in ST and SC data with annotated celltype**

```
ST = sc.read_visium(path/to/data_folder, count_file=filename.h5, library_id = 'custom')
SC = anndata.read_h5ad(path/to/scRNA.h5ad)
```

**Read in ligand and targets you would like to work with**

```
pairs = pd.read_csv(path/to/ltpairs)
ligands = pairs['ligand']
targets = pairs['target']
```

**Get list of celltypes and estimated celltype proportions**

```
celltype_proportions = pd.read_csv(path/to/celltype_proportion.csv)
celltypes = list(celltype_proportions.columns)
```

**Read in mRNA abundance values generated from cell2location**

```
expins = pickle.load(cell2location/mRNA_abundance/gene_celltype_spot,'rb'))
genes = list(expins.keys())
expins_new = []
for gene in expins.keys():
    expins_new.append(expins[gene].to_numpy())
expins = np.array(expins_new)
```

> **NOTE**: cell type proportions and gene-cell type specific mRNA abundance have been generated via cell2location (v 0.06). To generate the gene-cell type specific mRNA abundance values, you can use 
> 
> `Renoir.compute_mRNA_abundance(model, genes)` 
> 
> where model is the cell2location model used to estimate cell type abundance values and genes are the list of genes you wish to calculate the mRNA abundance values for. This should only be considered if cell2location v0.6/v0.5 was utilized to estimate cell type abundance values.

**Compute neighborhood from ST data**

```
graph = Renoir.neighborhood(ST.obs['array_row'].tolist(), ST.obs['array_col'].tolist(), technology='visium')
```

**Get list of unique ligands, targets and ligand-target pairs**

```
ligand_target_index, ligand_target_pairs, ST_nonzero = Renoir.get_ligand_target(ligands, targets, ST, SC, genes)
```

**Get neighborhood scores**

```
neighborhood_scores = Renoir.compute_neighborhood_scores(SC, ST, celltypes, celltype_proportions, graph, ligand_target_index, ligand_target_pairs, ST_nonzero, expins, genes)
```

## References

References
