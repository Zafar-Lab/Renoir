import scanpy as sc
import numpy as np
import pandas as pd
import math
import multiprocessing as mp
from numba import njit, prange
import time


#Get unique set of ligands and targets and indexed ligand target pairs
#Returns: dict: {ligands and targets:index}; list: ligand x target
def get_ligand_target(ligands, targets, ST, SC):
    sc.pp.filter_cells(ST, min_genes=1)
    sc.pp.filter_genes(ST, min_cells=1)
    sc.pp.filter_cells(SC, min_genes=1)
    sc.pp.filter_genes(SC, min_cells=1)
    if len(ligands) != len(targets):
        raise Exception("ERROR: No. of ligands and targets should be equal")
    ligands_and_targets = []
    for gene in ligands+targets:
        if gene not in ligands_and_targets and gene in SC.var_names and gene in ST.var_names:
            ligands_and_targets.append(gene)
    ligands_and_targets = dict(zip(ligands_and_targets, range(len(ligands_and_targets))))
    ligand_target_pairs = []
    for index in range(len(ligands)):
        ligand_index = ligands_and_targets[ligands[index]]
        target_index = ligands_and_targets[targets[index]]
        if(ligand_index > len(ligand_target_pairs)-1):
            ligand_target_pairs.append([])
        if target_index not in ligand_target_pairs[ligand_index]:
            ligand_target_pairs[ligand_index].append(target_index)
    
    pad = len(max(ligand_target_pairs, key=len))
    return np.asarray([ligand + [-1]*(pad - len(ligand)) for ligand in ligand_target_pairs])
    return ligands_and_targets, ligand_target_pairs



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
    ligand_bins[0]=min(ligand_bins)-1
    target_bins[0]=min(target_bins)-1
    ligand_hist = np.digitize(ligand_data, ligand_bins,right=True)
    target_hist = np.digitize(target_data, target_bins,right=True)
    table = np.zeros((max(ligand_hist),max(target_hist)))
    for index in range(ligand_hist.shape[0]):
        table[ligand_hist[index]-1][target_hist[index]-1]+=1
    return table


'''
#Compute ISM using fast histogram2d
def ISM(ligand, target, n_celltypes, H, celltype_proportions, scdata, bins):
    MI_ligand = []
    for celltype in range(n_celltypes):
        table, _, _ = np.histogram2d(scdata[celltype][ligand], scdata[celltype][target], bins=(bins[ligand][celltype], bins[target][celltype]))
        Hl = H[ligand][celltype]
        Ht = H[target][celltype]
        ISM_temp = (Hl + Ht - entropy(table))/min(Hl, Ht)
        if ISM_temp <= 0:
            ISM_temp=1
        MI_ligand.append(-1*np.log10(ISM_temp))
    return (celltype_proportions*np.array(MI_ligand)).sum(axis=1)

#Parallelize ISM calculation
def ISM_wrapper(ligand, ltpairs, H, n_celltypes, celltype_proportions, scdata, bins, pool):
    result = pool.starmap(ISM, [(ligand, target, n_celltypes, H, celltype_proportions, scdata, bins) for target in range(len(ltpairs[ligand]))])
    return result
'''

@njit(parallel=True)
def ISM(ligand_target_pairs, H, n_celltypes, celltype_start_index, celltype_proportions, scdata):
    ISM_result = np.full((ligand_target_pairs.shape[0], ligand_target_pairs.shape[1], celltype_proportions.shape[0]), -1, dtype=np.float32)
    for ligand in prange(ligand_target_pairs.shape[0]):
        for target_index in prange(ligand_target_pairs[ligand].shape[0]):
            target = ligand_target_pairs[ligand][target_index]
            if target != -1:
                MI_ligand = np.zeros(n_celltypes, dtype=np.float32)
                for celltype in range(n_celltypes):
                    ligand_data = scdata[ligand][celltype_start_index[celltype]:celltype_start_index[celltype+1]]
                    target_data = scdata[target][celltype_start_index[celltype]:celltype_start_index[celltype+1]]
                    ligand_bin = np.linspace(np.min(scdata[ligand][celltype_start_index[celltype]:celltype_start_index[celltype+1]]),np.max(scdata[ligand][celltype_start_index[celltype]:celltype_start_index[celltype+1]]),math.ceil(math.sqrt(celltype_start_index[celltype+1] - celltype_start_index[celltype]))+1)
                    target_bin = np.linspace(np.min(scdata[target][celltype_start_index[celltype]:celltype_start_index[celltype+1]]),np.max(scdata[target][celltype_start_index[celltype]:celltype_start_index[celltype+1]]),math.ceil(math.sqrt(celltype_start_index[celltype+1] - celltype_start_index[celltype]))+1)
                    table = hist2d(ligand_data, target_data, ligand_bin, target_bin)
                    Hl = H[ligand][celltype]
                    Ht = H[target][celltype]
                    ISM_temp = (Hl + Ht - entropy(table))/min(Hl, Ht)
                    if ISM_temp <= 0 or np.isnan(ISM_temp) or np.isinf(ISM_temp):
                        ISM_temp=1
                    MI_ligand[celltype] = -1*np.log10(ISM_temp)
                ISM_result[ligand][target_index] = (celltype_proportions*MI_ligand).sum(axis=1)
    return ISM_result


#Fast histogram function
@njit
def fasthist(n_ligand_target, n_celltypes, celltype_start_index, scdata):
    H = np.zeros((n_ligand_target, n_celltypes))
    for gene in range(n_ligand_target):
        for celltype in range(n_celltypes):
            start = np.min(scdata[gene][celltype_start_index[celltype]:celltype_start_index[celltype+1]])
            stop = np.max(scdata[gene][celltype_start_index[celltype]:celltype_start_index[celltype+1]])
            num = math.ceil(math.sqrt(celltype_start_index[celltype+1] - celltype_start_index[celltype]))+1
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
    celltype_start_index = np.zeros(len(celltypes)+1, dtype=np.int8)
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
    '''
    with mp.Manager() as manager:
        with manager.Pool(n_proc) as pool_l:
            with manager.Pool(n_proc) as pool_t:
                start = time.time()
                ISM = pool_l.starmap(ISM_wrapper, [(ligand, ligand_target_pairs, H, len(celltypes), celltype_proportions, scdata, bins, pool_t) for ligand in range(len(ligand_target_pairs))])
                end = time.time()
                print('Time taken: ',(end-start))
    '''
    PEM_result = PEM(expins, celltype_proportions)
    return ISM_result, PEM_result



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
def min_max(PEM, ISM, graph, ligand_target_pairs):
    min_PEM = 100
    min_ISM = 100
    max_PEM = -100
    max_ISM = -100
    for lig_index in range(ligand_target_pairs.shape[0]):
        for tar_index in ligand_target_pairs[lig_index]:
            if tar_index == -1:
                break
            for spot in range(graph.shape[0]):
                for neighbor in graph[spot]:
                    if neighbor == -1:
                        break
                    PEM_sn = PEM[lig_index][spot]+PEM[tar_index][neighbor]
                    ISM_sn = ISM[lig_index][tar_index][spot] + ISM[lig_index][tar_index][neighbor]
                    if PEM_sn < min_PEM:
                        min_PEM = PEM_sn
                    if PEM_sn > max_PEM:
                        max_PEM = PEM_sn
                    if ISM_sn < min_ISM:
                        min_ISM = PEM_sn
                    if ISM_sn > max_ISM:
                        max_ISM = ISM_sn
    return min_PEM, max_PEM, min_ISM, max_ISM



#Compute the neighborhood scores given
#Return array of size ligand_target_pair vs spots
@njit(parallel=True)
def get_neighborhood_score(ligand_target_pairs, graph, PEM, ISM):
    neighborhood_scores = np.zeros((ligand_target_pairs.shape[0], ligand_target_pairs.shape[1], graph.shape[0]), dtype=np.float32)
    min_PEM, max_PEM, min_ISM, max_ISM = min_max(PEM,ISM,graph,ligand_target_pairs)
    for ligand in prange(ligand_target_pairs.shape[0]):
        for target_index in prange(ligand_target_pairs.shape[1]):
            target = ligand_target_pairs[ligand][target_index]
            if target!=-1:
                for spot in range(graph.shape[0]):
                    PEM_lt = (PEM[ligand][spot]+PEM[target][spot])/2
                    ISM_lt = (ISM[ligand][target][spot]+ISM[ligand][target][spot])/2
                    cumsum = ((PEM_lt - min_PEM)/(max_PEM - min_PEM)) * ((ISM_lt - min_ISM)/(max_ISM - min_ISM))
                    for neighbor in graph[spot]:
                        if neighbor!=-1:
                            PEM_lt1 = (PEM[ligand][spot]+PEM[target][neighbor])/2
                            PEM_lt2 = (PEM[ligand][neighbor]+PEM[target][spot])/2
                            ISM_lt1 = (ISM[ligand][target][spot]+ISM[ligand][target][neighbor])/2
                            ISM_lt2 = (ISM[ligand][target][neighbor]+ISM[ligand][target][spot])/2
                            cumsum += ((PEM_lt1 - min_PEM)/(max_PEM - min_PEM)) * ((ISM_lt1 - min_ISM)/(max_ISM - min_ISM))
                            cumsum += ((PEM_lt2 - min_PEM)/(max_PEM - min_PEM)) * ((ISM_lt2 - min_ISM)/(max_ISM - min_ISM))
                    neighborhood_scores[ligand][target][spot] = cumsum
    return neighborhood_scores

