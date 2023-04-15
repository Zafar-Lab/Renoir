import anndata
import math
import numpy as np
from numba import njit, prange
import pandas as pd


#Get unique set of ligands and targets and indexed ligand target pairs
#Returns: dict: {ligands and targets:index}; array: ligand x target; nonzero ST positions
def get_ligand_target(ligands, targets, ST, SC, expins_genes):
    SC_copy = SC.copy()
    ST_copy = ST.copy()
    if len(ligands) != len(targets):
        raise Exception("ERROR: No. of ligands and targets should be equal")
    
    ligands_and_targets = []
    ligands_subset = []
    targets_subset = []

    for index in range(len(ligands)):
        if (ligands[index] in SC_copy.var_names and ligands[index] in ST_copy.var_names and ligands[index] in expins_genes) and (targets[index] in SC_copy.var_names and targets[index] in ST_copy.var_names and targets[index] in expins_genes):
            ligands_subset.append(ligands[index])
            targets_subset.append(targets[index])
    
    for gene in ligands_subset+targets_subset:
        if gene not in ligands_and_targets:
            ligands_and_targets.append(gene)
    ligands_and_targets = dict(zip(ligands_and_targets, range(len(ligands_and_targets))))
    ligand_target_pairs = []
    for index in range(len(ligands_subset)):
        ligand_index = ligands_and_targets[ligands_subset[index]]
        target_index = ligands_and_targets[targets_subset[index]]
        if(ligand_index > len(ligand_target_pairs)-1):
            ligand_target_pairs.append([])
        if target_index not in ligand_target_pairs[ligand_index]:
            ligand_target_pairs[ligand_index].append(target_index)
    
    ST_nonzero = np.array(np.transpose(ST[:,list(ligands_and_targets.keys())].X.toarray()), dtype=bool)
    pad_ligand_target = len(max(ligand_target_pairs, key=len))
    return ligands_and_targets, np.asarray([ligand + [-1]*(pad_ligand_target - len(ligand)) for ligand in ligand_target_pairs]), ST_nonzero


#Generate neighbors for each spot given (X,Y) for each spot and spatial technology used
#Returns: graph: row index = spot index; graph[row index] = list of spot neighbors 
@njit(parallel=True, fastmath=True)
def create_graph(X_coord, Y_coord, technology, radius=0):
    graph = [ [-1] for _ in range(len(X_coord)) ]
    if technology == 'visium' or technology == 'ST':
        for spot1 in prange(len(X_coord)):
            for spot2 in prange(len(X_coord)):
                if (spot1 == spot2) or (abs(X_coord[spot1] - X_coord[spot2]) == 1 and abs(Y_coord[spot1] - Y_coord[spot2]) == 1) or (X_coord[spot1] == X_coord[spot2] and abs(Y_coord[spot1] - Y_coord[spot2]) == 2):
                    graph[spot1].append(spot2)
    
    elif technology == 'slideseq':
        for spot1 in prange(len(X_coord)):
            for spot2 in prange(len(X_coord)):
                if math.sqrt((X_coord[spot1] - X_coord[spot2])**2 + (Y_coord[spot1] - Y_coord[spot2])**2) <= radius:
                    graph[spot1].append(spot2)
    
    for spot in range(len(X_coord)):
        del graph[spot][0]

    pad = max([len(spot) for spot in graph])
    for spot in range(len(graph)):
        for _ in range(pad - len(graph[spot])):
            graph[spot].append(-1)

    return np.array(graph)

#Wrapper function for create_graph
#Returns: graph: row index = spot index; graph[row index] = list of spot neighbors
def neighborhood(X_coord, Y_coord, technology, radius):
    technologies=['visium','slideseq']
    if technology in technologies:
        graph = create_graph(X_coord, Y_coord, technology, radius)
        return graph
    else:
        raise Exception("technology currently not supported.")


#Calculate miller-madow entropy
@njit
def entropy(values):
    values = values.flatten()
    values = values[values != 0]
    n = np.sum(values)
    return -1*np.sum((values/n)*np.log(values/n)) + ((len(values)-1)/(2*n))


#Custom histogram2d
@njit 
def hist2d(ligand_data, target_data, ligand_bins, target_bins):
    ligand_bins += 1e-15
    ligand_bins[0] -= 1e-15
    target_bins += 1e-15
    target_bins[0] -= 1e-15
    ligand_hist = np.digitize(ligand_data, ligand_bins, right=False)
    ligand_hist[ligand_hist > ligand_bins.shape[0]-1] = ligand_bins.shape[0]-1
    target_hist = np.digitize(target_data, target_bins, right=False)
    target_hist[target_hist > target_bins.shape[0]-1] = target_bins.shape[0]-1
    table = np.zeros((ligand_bins.shape[0]-1,target_bins.shape[0]-1))
    for index in range(ligand_data.shape[0]):
        table[ligand_hist[index]-1][target_hist[index]-1]+=1
    return table



#@njit(parallel=True)
def ISM(ligand_target_pairs, H, n_celltypes, celltype_start_index, celltype_proportions, scdata):
    ISM_result = np.full((ligand_target_pairs.shape[0], ligand_target_pairs.shape[1], celltype_proportions.shape[0]), -1, dtype=np.float64)
    for ligand in range(ligand_target_pairs.shape[0]):
        for target_index in range(ligand_target_pairs[ligand].shape[0]):
            target = ligand_target_pairs[ligand][target_index]
            if target != -1:
                MI_ligand = np.zeros(n_celltypes, dtype=np.float64)
                for celltype in range(n_celltypes):
                    ligand_data = scdata[ligand][celltype_start_index[celltype]:celltype_start_index[celltype+1]]
                    target_data = scdata[target][celltype_start_index[celltype]:celltype_start_index[celltype+1]]
                    ligand_bin = np.linspace(np.min(scdata[ligand][celltype_start_index[celltype]:celltype_start_index[celltype+1]]),np.max(scdata[ligand][celltype_start_index[celltype]:celltype_start_index[celltype+1]]),math.ceil(math.sqrt(celltype_start_index[celltype+1] - celltype_start_index[celltype]))+1)
                    target_bin = np.linspace(np.min(scdata[target][celltype_start_index[celltype]:celltype_start_index[celltype+1]]),np.max(scdata[target][celltype_start_index[celltype]:celltype_start_index[celltype+1]]),math.ceil(math.sqrt(celltype_start_index[celltype+1] - celltype_start_index[celltype]))+1)
                    table = hist2d(ligand_data, target_data, ligand_bin, target_bin)
                    Hl = H[ligand][celltype]
                    Ht = H[target][celltype]
                    if min(Hl,Ht)==0:
                        ISM_temp = 1
                    else:
                        ISM_temp = (Hl + Ht - entropy(table))/min(Hl, Ht)
                        if ISM_temp <= 0 or np.isnan(ISM_temp) or np.isinf(ISM_temp):
                            ISM_temp=1
                    MI_ligand[celltype] = -1*np.log10(ISM_temp)
                ISM_result[ligand][target_index] = (celltype_proportions*MI_ligand).sum(axis=1)
    return ISM_result

'''
@njit
def MI(ligand_target_pairs, H, n_celltypes, celltype_start_index, celltype_proportions, scdata):
    ISM_result = np.full((ligand_target_pairs.shape[0], ligand_target_pairs.shape[1], celltype_proportions.shape[1]), -1, dtype=np.float64)
    for ligand in range(ligand_target_pairs.shape[0]):
        for target_index in range(ligand_target_pairs[ligand].shape[0]):
            target = ligand_target_pairs[ligand][target_index]
            if target != -1:
                MI_ligand = np.zeros(n_celltypes, dtype=np.float64)
                for celltype in range(n_celltypes):
                    ligand_data = scdata[ligand][celltype_start_index[celltype]:celltype_start_index[celltype+1]]
                    target_data = scdata[target][celltype_start_index[celltype]:celltype_start_index[celltype+1]]
                    ligand_bin = np.linspace(np.min(scdata[ligand][celltype_start_index[celltype]:celltype_start_index[celltype+1]]),np.max(scdata[ligand][celltype_start_index[celltype]:celltype_start_index[celltype+1]]),math.ceil(math.sqrt(celltype_start_index[celltype+1] - celltype_start_index[celltype]))+1)
                    target_bin = np.linspace(np.min(scdata[target][celltype_start_index[celltype]:celltype_start_index[celltype+1]]),np.max(scdata[target][celltype_start_index[celltype]:celltype_start_index[celltype+1]]),math.ceil(math.sqrt(celltype_start_index[celltype+1] - celltype_start_index[celltype]))+1)
                    table = hist2d(ligand_data, target_data, ligand_bin, target_bin)
                    Hl = H[ligand][celltype]
                    Ht = H[target][celltype]
                    if min(Hl,Ht)!=0:
                        ISM_temp = (Hl + Ht - entropy(table))/min(Hl, Ht)
                    else:
                        ISM_temp = 1
                    if ISM_temp <= 0 or np.isnan(ISM_temp) or np.isinf(ISM_temp):
                        ISM_temp = 1
                    MI_ligand[celltype] = -1*np.log10(ISM_temp)
                ISM_result[ligand][target_index] = MI_ligand
    
    return ISM_result
'''


#Fast histogram function
@njit
def fasthist(n_ligand_target, n_celltypes, celltype_start_index, scdata):
    H = np.zeros((n_ligand_target, n_celltypes))
    for gene in range(n_ligand_target):
        for celltype in range(n_celltypes):
            start = np.min(scdata[gene][celltype_start_index[celltype]:celltype_start_index[celltype+1]])
            stop = np.max(scdata[gene][celltype_start_index[celltype]:celltype_start_index[celltype+1]])
            num = math.ceil(np.sqrt(celltype_start_index[celltype+1] - celltype_start_index[celltype]))+1
            bin_temp = np.linspace(start,stop,num)
            bin_temp += 1e-15
            bin_temp[0] -= 1e-15
            hist, _ = np.histogram(scdata[gene][celltype_start_index[celltype]:celltype_start_index[celltype+1]], bins=bin_temp)
            H[gene][celltype] = entropy(hist)
    return H



#Calculate ISM and PEM for ligand targets considered
#Returns: array, ISM (l x t x s) and PEM (g x s)
def ISM_PEM(SC, expins, ligand_target_list, ligand_target_pairs, celltypes, celltype_proportions):
    scdata = np.empty((len(ligand_target_list), 0))
    celltype_start_index = np.zeros(len(celltypes)+1, dtype=np.int32)
    for celltype in range(len(celltypes)):
        celltype_start_index[celltype] = scdata.shape[1]
        scdata = np.concatenate((scdata, np.transpose(SC[SC.obs.index[SC.obs['celltype'] == celltypes[celltype]], ligand_target_list].X.toarray())), axis=1)
    celltype_start_index[len(celltypes)] = scdata.shape[1]
    for gene in range(len(ligand_target_list)):
        for celltype in range(len(celltypes)):
            if np.sum(scdata[gene][celltype_start_index[celltype]:celltype_start_index[celltype+1]]) == 0:
                index = np.random.randint(scdata[gene][celltype_start_index[celltype]:celltype_start_index[celltype+1]].shape[0])
                scdata[gene][celltype_start_index[celltype] + index] = 1e-20
    H = fasthist(len(ligand_target_list), len(celltypes), celltype_start_index, scdata)
    ISM_result = ISM(ligand_target_pairs, H, len(celltypes), celltype_start_index, celltype_proportions, scdata)
    ECS_result = PEM(expins, celltype_proportions)
    return ISM_result, ECS_result



#Calculate PEM
#Returns PEM of dimension g x s
@njit
def PEM(expins, celltype_proportions):
    S = expins.sum(axis=0)
    S_total = S.sum(axis=1)
    S = S/S_total.reshape((-1, 1))
    shape = S.shape
    S = S.ravel()
    S[np.isnan(S)] = 0
    S = S.reshape(shape)
    exp_total = expins.sum(axis=2).reshape((expins.shape[0], expins.shape[1], 1))
    E = S*exp_total
    PEM = np.log10(expins/E)
    shape = PEM.shape
    PEM = PEM.ravel()
    PEM[np.isnan(PEM)] = 0
    PEM[np.isinf(PEM)] = 0
    PEM = PEM.reshape(shape)
    PEM = (PEM*celltype_proportions).sum(axis=2)
    return PEM



#Find min/max PEM/ISM
@njit
def min_max(PEM, ISM, graph, ligand_target_pairs, ST_nonzero):
    min_PEM = 100
    min_ISM = 100
    max_PEM = -100
    max_ISM = -100
    for lig_index in range(ligand_target_pairs.shape[0]):
        for tar_index in range(ligand_target_pairs.shape[1]):
            target = ligand_target_pairs[lig_index][tar_index]
            if target == -1:
                break
            for spot in range(graph.shape[0]):
                for neighbor in graph[spot]:
                    if neighbor == -1:
                        break
                    if ST_nonzero[lig_index][spot] and ST_nonzero[target][neighbor]:
                        PEM_sn = PEM[lig_index][spot]+PEM[target][neighbor]
                        ISM_sn = ISM[lig_index][tar_index][spot] + ISM[lig_index][tar_index][neighbor]
                        if PEM_sn < min_PEM:
                            min_PEM = PEM_sn
                        if PEM_sn > max_PEM:
                            max_PEM = PEM_sn
                        if ISM_sn < min_ISM:
                            min_ISM = ISM_sn
                        if ISM_sn > max_ISM:
                            max_ISM = ISM_sn
    return min_PEM/2, max_PEM/2, min_ISM/2, max_ISM/2


#Compute the neighborhood scores given
#Return array of size ligand_target_pair vs spots
@njit(parallel=True)
def get_neighborhood_score(ligand_target_pairs, graph, PEM, ISM, ST_nonzero):
    neighborhood_scores = np.zeros((ligand_target_pairs.shape[0], ligand_target_pairs.shape[1], graph.shape[0]), dtype=np.float32)
    min_PEM, max_PEM, min_ISM, max_ISM = min_max(PEM,ISM,graph,ligand_target_pairs,ST_nonzero)
    for ligand in prange(ligand_target_pairs.shape[0]):
        for target_index in prange(ligand_target_pairs.shape[1]):
            target = ligand_target_pairs[ligand][target_index]
            if target!=-1:
                for spot in range(graph.shape[0]):
                    cumsum = 0
                    if ST_nonzero[ligand][spot] and ST_nonzero[target][spot]:
                        PEM_lt = (PEM[ligand][spot]+PEM[target][spot])/2
                        ISM_lt = (ISM[ligand][target_index][spot]+ISM[ligand][target_index][spot])/2
                        cumsum = ((PEM_lt - min_PEM)/(max_PEM - min_PEM)) * ((ISM_lt - min_ISM)/(max_ISM - min_ISM))
                    for neighbor in graph[spot]:
                        if neighbor!=-1 and neighbor!=spot:
                            if ST_nonzero[ligand][spot] and ST_nonzero[target][neighbor]:
                                PEM_lt1 = (PEM[ligand][spot]+PEM[target][neighbor])/2
                                ISM_lt1 = (ISM[ligand][target_index][spot]+ISM[ligand][target_index][neighbor])/2
                                cumsum += ((PEM_lt1 - min_PEM)/(max_PEM - min_PEM)) * ((ISM_lt1 - min_ISM)/(max_ISM - min_ISM))
                            if ST_nonzero[ligand][neighbor] and ST_nonzero[target][spot]:
                                PEM_lt2 = (PEM[ligand][neighbor]+PEM[target][spot])/2
                                ISM_lt2 = (ISM[ligand][target_index][neighbor]+ISM[ligand][target_index][spot])/2
                                cumsum += ((PEM_lt2 - min_PEM)/(max_PEM - min_PEM)) * ((ISM_lt2 - min_ISM)/(max_ISM - min_ISM))
                    neighborhood_scores[ligand][target_index][spot] = cumsum
    return neighborhood_scores

#Wrapper function for get neighborhood scores
#Return neighborhood scores: array of size ligand_target_pair vs spots OR anndata object
def compute_neighborhood_scores(SC, ST, celltypes, celltype_proportions, graph, ligand_target_index, ligand_target_pairs, ST_nonzero, expins, genes, return_adata=True):
    #Get gene indices
    gene_indices = []
    ligand_target_list = list(ligand_target_index.keys())
    for x in ligand_target_list:
        gene_indices.append(genes.index(x))
    #Compute ECS and ISM values
    ISM_result, ECS_result = ISM_PEM(SC, expins, list(ligand_target_index.keys()), ligand_target_pairs, celltypes, celltype_proportions)
    ECS_result = ECS_result[gene_indices]
    neighborhood_scores = get_neighborhood_score(ligand_target_pairs, graph, ECS_result, ISM_result, ST_nonzero)
    if return_adata:
        adata = toadata(neighborhood_scores, ST, ligand_target_pairs, ligand_target_list)
        return adata
    else:
        return neighborhood_scores


#Convert numpy neighborhood scores to anndata
#Return anndata object with pairs vs spots
def toadata(neighborhood_scores, ST, ligand_target_pairs, ligand_target_list):
    dataframe = []
    pairs = []
    for lig_index in range(ligand_target_pairs.shape[0]):
        for tar_index in range(ligand_target_pairs.shape[1]):
            if neighborhood_scores[lig_index][tar_index].sum() > 0:
                target = ligand_target_pairs[lig_index][tar_index]
                pairs.append(ligand_target_list[lig_index]+':'+ligand_target_list[target])
                dataframe.append(list(neighborhood_scores[lig_index][tar_index]))
    dataframe = pd.DataFrame(dataframe).transpose()
    dataframe.columns = pairs
    dataframe.index = ST.obs_names
    adata = anndata.AnnData(dataframe)
    adata.obs = ST.obs
    adata.var['gene_ids'] = pairs
    adata.uns = ST.uns
    adata.obsm = ST.obsm
    return adata

