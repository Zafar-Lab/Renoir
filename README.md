# Renoir

**Renoir** is an information-theory-based scoring metric for quantifying the activity between a ligand and its target gene given a specific spatial context. Renoir can also infer spatial communication domains that harbor similar ligand-target activities. Renoir also identifies spatially enriched ligand-target gene sets (pathway activity) and characterizes domain-specific activities between ligands and targets.

<p align="center">
<img src="https://github.com/Zafar-Lab/Renoir/blob/main/docs/source/images/Overview.png" alt="Overview of RENOIR" width="600"/>
</p>

### Requirements
All requirements are provided in the ```renoir.yml``` file. It is recommended to utilize the same versions as provided in ```renoir.yml``` file.

### Installation

- Download the repository and install the conda environment using `conda env create -f renoir.yml`
- Install Renoir by cd to this directory and running `pip install`.

**Note**: Installation does not include **cell2location**. Please install cell2location separately.

### Usage

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

> **NOTE**: cell type proportions and cell type-specific mRNA abundances for genes have been inferred via cell2location (v 0.06). To quantify the cell type-specific mRNA abundance values for genes, you can use 
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

### Documentation

Documentation for Renoir is available at: https://renoir.readthedocs.io/en/latest/

### References

References
