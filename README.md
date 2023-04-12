# Renoir

Renoir is an information-theory-based scoring metric that delineates spatial communication domains that are cell-type, ligand, and target specific. It also identifies spatially enriched ligand-target gene sets and characterizes domain-specific activities between ligands and targets.

<p align="center">
<img src="https://github.com/narein97/Renoir/blob/main/images/Overview.png" alt="Overview of RENOIR" width="600"/>
</p>

### Requirements
**Note**: All requirements are provided in the renoir.yml file. Specific libraries required for the neighborhood score computation and downstream tasks are as follows,

- scanpy (1.9.1)
- numpy (1.21.6)
- numba (0.56.4)
- pandas (1.2.5)
- scipy (1.7.1)
- hdbscan (0.8.28)
- dynamictreecut (0.1.1)
- sklearn (0.24.2)
- seaborn (0.11.1)
- matplotlib (0.11.7)

### Installation
**Note**: Installation does not include cell2location. Please install cell2location separately.

- Download the repository and install the conda environment using `conda env create -f renoir.yml`
- Install Renoir by cd to this directory and running `pip install .`

### Usage

```
import Renoir
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
import pickle
```

**Load in spatial transcriptomics and corresponding single cell data**

```
#Load in spatial transcriptomics data
ST = sc.read_visium(path_to_visium_data, count_file=count_file.h5, library_id = 'custom')

#Load in celltype annotated single cell data
SC = anndata.read_h5ad(path_to_scRNA-seq_h5ad.h5ad)
```

**Load in list of ligand and target pairs to generate neighbourhood scores for**

```
#Get ligand-target pairs you would like to work with
ligand_target = pd.read_csv(path_to_ligand_target_list)
ligands = ligand_target['ligand'].tolist()
targets = ligand_target['target'].tolist()

#Get list of unique, indexed ligands and targets and ligand-target pairs
ligand_target_index, ligand_target_pairs = get_ligand_target(ligands, targets, ST, SC)
ligand_target_list = list(ligand_target_index.keys())
```

**Load in estimated cell type proportions and gene-cell type specific mRNA abundance**

```
#Get list of celltypes and estimated celltype proportions and mRNA abundance values
celltype_proportions = pd.read_csv(path_to_celltype_abundance)
celltype_proportions['sum']=celltype_proportions.sum(axis=1)
celltype_proportions=celltype_proportions[celltype_proportions.columns].div(celltype_proportions['sum'], axis=0)
celltype_proportions=celltype_proportions.drop('sum',axis=1)
celltypes = list(celltype_proportions.columns)
celltype_proportions = celltype_proportions.to_numpy()

#Read in mRNA abundance values generated from cell2location
expins = pickle.load(open(path_to_computed_mRNA_abundance,'rb'))
expins_new = []
for gene in ligand_target_list:
    expins_new.append(expins[gene].to_numpy())
expins = np.array(expins_new)
```

> **NOTE**: cell type proportions and gene-cell type specific mRNA abundance have been generated via cell2location (v 0.06). To generate the gene-cell type specific mRNA abundance values, you can use 
> 
> `Renoir.compute_mRNA_abundance(model, genes)` 
> 
> where model is the cell2location model used to estimate cell type abundance values and genes are the list of genes you wish to calculate the mRNA abundance values for.

**Compute neighbourhood scores for the selected set of ligand-target pairs**

```
#Compute graph for ST data
graph = create_graph(ST.obs['array_row'].tolist(), ST.obs['array_col'].tolist(), technology='visium')

#Compute ECS and ISM values
ISM_result, ECS_result = ISM_PEM(SC, expins, ligand_target_list, ligand_target_pairs, celltypes, celltype_proportions)
ECS_result = ECS_result[gene_indices]

#Get neighborhood scores
neighborhood_scores = get_neighborhood_score(ligand_target_pairs, graph, ECS_result, ISM_result)
```

### Documentation

Documentation is available at: link

### References

References
