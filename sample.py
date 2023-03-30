import Renoir
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
import pickle

#Read in ST and SC data
ST = sc.read_visium('/home/narein/Breast_Cancer_corrected/ST', count_file='ST.h5', library_id = 'custom')
SC = anndata.read_h5ad('/home/narein/Breast_Cancer_corrected/scRNA.h5ad')

#Get ligand and targets you would like to work with
temp = sc.read_visium('/home/narein/Breast_Cancer_corrected/ST', count_file='BC_neighborhood_scores_top_10_corrected.h5', library_id = 'custom')
sc.pp.filter_genes(temp, min_cells=1)
sc.pp.filter_cells(temp, min_genes=1)
ligands = []
targets = []
for pair in temp.var_names:
    ligands.append(pair.split(':')[0])
    targets.append(pair.split(':')[1])


#Get list of celltypes and estimated celltype proportions
celltype_proportions = pd.read_csv('/home/narein/Breast_Cancer_corrected/cell2loc_result/rep1/W_cell_density_q05.csv',index_col=0)
celltype_proportions.columns = celltype_proportions.columns.str.replace('q05_spot_factors', '')
celltype_proportions['sum']=celltype_proportions.sum(axis=1)
celltype_proportions=celltype_proportions[celltype_proportions.columns].div(celltype_proportions['sum'], axis=0)
celltype_proportions=celltype_proportions.drop('sum',axis=1)
celltypes = list(celltype_proportions.columns)
celltype_proportions = celltype_proportions.to_numpy()


#Compute graph for ST data
graph = create_graph(ST.obs['array_row'].tolist(), ST.obs['array_col'].tolist(), technology='visium')
#Get list of unique ligands and targets and ligand-target pairs
ligand_target_index, ligand_target_pairs = get_ligand_target(ligands, targets, ST, SC)

#Read in mRNA abundance values generated from cell2location
expins = pickle.load(open('/home/narein/Breast_Cancer_corrected/cell2loc_result/rep1/mRNA_subset_100.pkl','rb'))
expins_new = []
for gene in ligand_target_index.keys():
    expins_new.append(expins[gene].to_numpy())
expins = np.array(expins_new)

#Compute ECS and ISM values
ISM_result, ECS_result = ISM_PEM(SC, expins, list(ligand_target_index.keys()), ligand_target_pairs, celltypes, celltype_proportions)
#Get neighborhood scores
neighborhood_scores = get_neighborhood_score(ligand_target_pairs, graph, ECS_result, ISM_result)