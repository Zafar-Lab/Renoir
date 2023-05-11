### hdbscan clusters  ++
### Using neighborhood score clusters (different heatmaps for WP_H and geneCluster)  ++
### pcscore vs neighborhood scores ++
### Spot vs Spot ++
### celltype vs lt pairs average
### Sankey plot ++

# from turtle import color
import pandas as pd
import scanpy as sc
import numpy as np
import random
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage
# from scipy import sparse
import hdbscan
from dynamicTreeCut import cutreeHybrid
from sklearn.decomposition import TruncatedSVD
import matplotlib.backends.backend_pdf
from sklearn.metrics.pairwise import cosine_similarity
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import seaborn as sns
import distinctipy
from matplotlib.patches import Patch

# create hdbscan/dynamic hclust gene clusters
def get_msig(species, path=None):
    if species == 'human':
        return pd.read_csv('./msigdb/msig_human_WP_H_KEGG')
    elif species == 'mouse':
        return pd.read_csv('./msigdb/msig_mouse_WP_H_KEGG')
    elif path is not None:
        return pd.read_csv(path)
    else:
        raise 'get_msig :: Invalid arg: provide a valid species (mouse/human)/path to the msigdb dataframe'


def create_cluster(neighbscore, msig, pathway_path=None, method='dhc', minpts=3, use_pathway=False, pathway_thresh=4, ltclust_thresh=6, restrict_to_KHW=False):
    """Creation of ligand-target pair sets using provided pathways; or by De novo clustering of ligand-target pairs.

    :param neighbscore: Neighborhood scores generated
    :type neighbscore: AnnData
    :param msig: Dataframe with columns 'gs_name' and 'gs_symbol'
    :type msig: DataFrame
    :param pathway_path: Subset of pathways to use within msig, defaults to None
    :type pathway_path: str, optional
    :param method: Clustering method for De novo clusters (available: ['hdbscan','dhc']), defaults to 'dhc'
    :type method: str, optional
    :param minpts: Minimum cluster size for hdbscan option, defaults to 3
    :type minpts: int, optional
    :param use_pathway: Use pathway_path to subset msig, defaults to False
    :type use_pathway: bool, optional
    :param pathway_thresh: Minimum number of ligand-target pairs in pathway ligand-target pair set, defaults to 4
    :type pathway_thresh: int, optional
    :param ltclust_thresh: Minimum number of ligand-target pairs in De novo ligand-target pair set, defaults to 6
    :type ltclust_thresh: int, optional
    :param restrict_to_KHW: Restrict pathways to KEGG, HALLMARK and WIKIPATHWAYS (Note: gs_name should have prefixes 'KEGG\_', 'HALLMARK\_', 'WP\_' in order to be considered), defaults to False
    :type restrict_to_KHW: bool, optional
    :return: Dictionary of pathways/De novo clusters and their corresponding ligand-target pair sets
    :rtype: dict
    """
    neighbscore_copy = neighbscore.copy()
    if use_pathway:
        pathwayset = pd.read_csv(pathway_path)
    sc.pp.filter_genes(neighbscore_copy, min_cells=1)
    neighbscore_copy = neighbscore_copy.to_df()
    #neighbscore_copy = neighbscore_copy.apply(min_max, axis=0)
    genes = []
    for ltpair in neighbscore_copy.columns:
        genes += ltpair.split(':')
    msig_subset = msig.loc[msig['gene_symbol'].isin(genes),]
    pathways = {}
    pathways = {k: g["gene_symbol"].tolist() for k, g in msig_subset[['gs_name', 'gene_symbol']].groupby("gs_name")}
    if use_pathway:
        pathways = {k: v for k, v in pathways.items() if k in pathwayset['pathways'].tolist()}
    for key, val in pathways.items():
        temp = []
        for genea in val:
            for geneb in val:
                if genea + ':' + geneb in neighbscore_copy.columns:
                    temp.append(genea + ':' + geneb)
        pathways[key] = list(set(temp))

    # Create distance matrix using pearson correlation
    if method == 'hdbscan':
        dist = pdist(abs(neighbscore_copy.corr() - 1))
        dist = squareform(dist)
        clusterer = hdbscan.HDBSCAN(min_cluster_size=minpts)
        clusterer.fit(dist)
        for i in range(len(clusterer.labels_)):
            if not ('cluster_' + str(clusterer.labels_[i]) in pathways.keys()) and clusterer.labels_[i] != -1:
                pathways['cluster_' + str(clusterer.labels_[i])] = []
                pathways['cluster_' + str(clusterer.labels_[i])].append(neighbscore_copy.columns[i])
            else:
                if clusterer.labels_[i] != -1:
                    pathways['cluster_' + str(clusterer.labels_[i])].append(neighbscore_copy.columns[i])
    elif method == 'dhc':
        dist = pdist(abs(neighbscore_copy.corr() - 1))
        dist = squareform(dist)
        link = linkage(dist, "average")
        clusters = cutreeHybrid(link, dist)
        for i in range(len(clusters['labels'])):
            if not ('cluster_' + str(clusters['labels'][i]) in pathways.keys()):
                pathways['cluster_' + str(clusters['labels'][i])] = []
                pathways['cluster_' + str(clusters['labels'][i])].append(neighbscore_copy.columns[i])
            else:
                pathways['cluster_' + str(clusters['labels'][i])].append(neighbscore_copy.columns[i])
    elif method is not None:
        raise 'create_cluster :: Invalid arg: provide a valid method (dhc/hdbscan)'

    # Consider those pathways / ltpair clusters with atleast threshold no. of pairs
    if restrict_to_KHW:
        pathways = {k: v for k, v in pathways.items() if ((k.startswith('KEGG') or k.startswith('HALLMARK') or k.startswith('WP_')) and len(v) > pathway_thresh) or (k.startswith('cluster_') and len(v) > ltclust_thresh)}
    else:
        pathways = {k: v for k, v in pathways.items() if (len(v) > pathway_thresh) or (k.startswith('cluster_') and len(v) > ltclust_thresh)}
    return pathways


def min_max(array):
    return (array - min(array)) / (max(array) - min(array))


def sum_norm(array):
    return array / sum(array)


def get_top_n_clust_pairs(neighbscore, neighbscore_df, n=20):
    top_n = []
    for cluster in sorted(neighbscore.obs['leiden'].unique()):
        top_spot = {}
        top_counts = {}
        cluster_spots = neighbscore.obs.loc[neighbscore.obs['leiden'] == str(cluster),].index
        for spot in cluster_spots:
            top_spot[spot] = list(neighbscore_df.loc[neighbscore_df.loc[:,spot] > 0, spot].sort_values(ascending=False).head(5).index)
        for spot in top_spot.keys():
            for index in range(len(top_spot[spot])):
                if (top_spot[spot][index] not in top_counts.keys()):
                    top_counts[top_spot[spot][index]] = 1
                else:
                    top_counts[top_spot[spot][index]] += 1
        top_n += sorted(top_counts, key=top_counts.get, reverse=True)[:n]
    return list(dict.fromkeys(top_n))
    # return list(set(top_n))


def downstream_analysis(neighbscore, ltpair_clusters=None, resolution=0.8, n_markers=20, n_top=20, pdf_path=None, return_cluster = False, return_pcs = False):
    """Performs leiden clustering over neighborhood scores, generates DE ligand-target pairs for each domain, finds top n highly expressed pairs across each domain and converts pathway/De novo cluster scores into AnnData objects. (Note: This function executes default DE analysis as is mentioned in scanpy docs. For more flexibility, it is recommended for users to perform analysis over generated neighborhood scores manually with the help of scanpy)

    :param neighbscore: adata object with pre computed neighborhood scores
    :type neighbscore: AnnData
    :param ltpair_clusters: Dictionary of pathway/gene set vs ligand-target pairs. Can be generated using create_cluster, defaults to None
    :type ltpair_clusters: dict, optional
    :param resolution: resoluton parameter for leiden clustering, defaults to 0.8
    :type resolution: float, optional
    :param n_markers: number of markers to display for each cluster, defaults to 20
    :type n_markers: int, optional
    :param n_top: number of highly expressed ligand-target pairs to display, defaults to 20
    :type n_top: int, optional
    :param pdf_path: path where the results are to be stored as a pdf, defaults to None
    :type pdf_path: str, optional
    :param return_cluster: return adata object with computed leiden clusters, defaults to False
    :type return_cluster: bool, optional
    :param return_pcs: return adata object with pathway/geneset and their corresponding scores, defaults to False
    :type return_pcs: bool, optional
    :return: return adata neighborhood score adata object with leiden clusters and/or adata object with pathway/geneset and their corresponding scores
    :rtype: AnnData
    """
    neighbscore_copy = neighbscore.copy()
    sc.pp.filter_genes(neighbscore_copy, min_cells=1)
    neighbscore_df = neighbscore_copy.to_df().T
    #neighbscore_df = neighbscore_df.apply(min_max, axis=1)
    if ltpair_clusters is None:
        raise 'ERROR: ltpair_clusters arg is None. Provide cluster object using create_cluster'
    pcs = {}
    for cluster, ltpairs in ltpair_clusters.items():
        pca = TruncatedSVD()
        pcs[cluster] = list(pca.fit_transform(neighbscore_df.loc[ltpairs,].T)[:, 0])
    pcs = pd.DataFrame.from_dict(pcs, orient='index')
    pcs.columns = neighbscore.obs_names
    pcs = pcs.applymap(abs)
    #pcs = pcs.apply(min_max, axis=1)
    pcs.fillna(0, inplace=True)
    neighbscore_copy = sc.AnnData(neighbscore_df.T).copy()
    neighbscore_copy.raw = neighbscore_copy
    neighbscore_df = neighbscore_copy.raw.to_adata().to_df().T
    pcs = sc.AnnData(pcs.T)
    pcs.raw = pcs
    neighbscore_copy.uns = neighbscore.uns
    neighbscore_copy.obs = neighbscore.obs
    neighbscore_copy.var['gene_ids'] = neighbscore_copy.var.index
    neighbscore_copy.obsm = neighbscore.obsm
    pcs.uns = neighbscore.uns
    pcs.obsm = neighbscore.obsm
    pcs.var['gene_ids'] = pcs.var.index
    pcs.obs = neighbscore.obs
    # Plot heatmaps for pc scores gained 
    # sc.pl.clustermap(sc.AnnData(pcs.to_df()[[x for x in pcs.var_names if not x.startswith('cluster')]].T),  show=True)
    # sc.pl.clustermap(sc.AnnData(pcs.to_df()[[x for x in pcs.var_names if x.startswith('cluster')]].T),  show=True)
    if pdf_path is not None:
        pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_path)
        if len([x for x in pcs.var_names if not x.startswith('cluster')])>2:
            sc.pl.clustermap(
                sc.AnnData(pcs.to_df()[[x for x in pcs.var_names if not x.startswith('cluster')]].T), show=False)
            temp_fig = plt.gcf()
            temp_fig.set_size_inches(25, 25)
            pdf.savefig(temp_fig)
            plt.close()
        if len([x for x in pcs.var_names if x.startswith('cluster')])>2:
            sc.pl.clustermap(
                sc.AnnData(pcs.to_df()[[x for x in pcs.var_names if x.startswith('cluster')]].T), show=False)
            temp_fig = plt.gcf()
            temp_fig.set_size_inches(25, 25)
            pdf.savefig(temp_fig)
            plt.close()
        # Perform clustering over PCs and neighborhood scores
    sc.pp.highly_variable_genes(pcs, n_top_genes=2000)
    sc.pp.highly_variable_genes(neighbscore_copy, n_top_genes=2000)
    sc.pp.scale(pcs)
    sc.pp.scale(neighbscore_copy)
    sc.tl.pca(pcs)
    sc.tl.pca(neighbscore_copy)
    sc.pp.neighbors(pcs, n_neighbors=20, n_pcs=pcs.varm['PCs'].shape[1])
    sc.pp.neighbors(neighbscore_copy, n_neighbors=20, n_pcs=neighbscore_copy.varm['PCs'].shape[1])
    # sc.tl.leiden(pcs, resolution = 0.3)
    sc.tl.leiden(neighbscore_copy, resolution=resolution)
    sc.tl.umap(pcs)
    sc.tl.umap(neighbscore_copy)
    # Set cluster labels of neighborhood scores to PCs
    pcs.uns['leiden'] = neighbscore_copy.uns['leiden']
    # Setup pdf to store images
    if pdf_path is not None:
        #raise "ERROR: pdf_path arg is None. Path to save pdf missing."
        # Plot UMAP / Spatial plot
        temp_fig = sc.pl.umap(pcs, components='all', color='leiden', return_fig=True, show=False)
        temp_fig.set_size_inches(25, 25)
        pdf.savefig(temp_fig)
        plt.close()
        temp_fig = sc.pl.umap(neighbscore_copy, components='all', color='leiden', return_fig=True, show=False)
        temp_fig.set_size_inches(25, 25)
        pdf.savefig(temp_fig)
        plt.close()
        sc.pl.spatial(neighbscore_copy, img_key="hires", color=["leiden"], size=1.4, show=False)
        temp_fig = plt.gcf()
        temp_fig.set_size_inches(25, 25)
        pdf.savefig(temp_fig)
        plt.close()
        # Compute top N markers per cluster
        sc.tl.rank_genes_groups(neighbscore_copy, "leiden", method="wilcoxon")
        sc.pl.rank_genes_groups_heatmap(neighbscore_copy, n_genes=n_markers, groupby="leiden", show_gene_labels=True, dendrogram = False,
                                        show=False, swap_axes=True, standard_scale='var')
        temp_fig = plt.gcf()
        temp_fig.set_size_inches(25, 25)
        pdf.savefig(temp_fig)
        plt.close()
        # Plot top N pairs per cluster
        top_n = get_top_n_clust_pairs(neighbscore_copy, neighbscore_df, n=n_top)
        sc.pl.heatmap(neighbscore_copy, var_names=top_n, groupby='leiden', show=False, swap_axes=True, show_gene_labels=True, standard_scale='var')
        temp_fig = plt.gcf()
        temp_fig.set_size_inches(25, 25)
        pdf.savefig(temp_fig)
        plt.close()
        pdf.close()
    if return_cluster and return_pcs:
        return neighbscore_copy, pcs
    else:
        if return_cluster:
            return neighbscore_copy
        if return_pcs:
            return pcs


def spot_v_spot(neighbscore, celltype, resolution=0.8, ltpair_clusters=None, n_markers=20, n_top=20, pdf_path=None):
    neighbscore_copy = neighbscore.copy()
    sc.pp.filter_genes(neighbscore_copy, min_cells=1)
    neighbscore_df = pd.DataFrame(neighbscore_copy.X.toarray())
    neighbscore_df.index = neighbscore_copy.obs_names
    neighbscore_df.columns = neighbscore_copy.var_names
    neighbscore_df = neighbscore_df.T
    neighbscore_df = neighbscore_df.apply(min_max, axis=1)
    if ltpair_clusters is None:
        raise 'ERROR: ltpair_clusters arg is None. Provide cluster object using create_cluster'
    pcs = {}
    for cluster, ltpairs in ltpair_clusters.items():
        pca = TruncatedSVD()
        pcs[cluster] = list(pca.fit_transform(neighbscore_df.loc[ltpairs,].T)[:, 0])
    pcs = pd.DataFrame.from_dict(pcs, orient='index')
    pcs.columns = neighbscore.obs_names
    pcs = pcs.applymap(abs)
    pcs = pcs.apply(min_max, axis=1)
    pcs.fillna(0, inplace=True)
    neighbscore_copy = sc.AnnData(neighbscore_df.T)
    neighbscore_copy.raw = neighbscore_copy
    neighbscore_df = neighbscore_copy.raw.to_adata().to_df().T
    pcs = sc.AnnData(pcs.T)
    neighbscore_copy.uns = neighbscore.uns
    neighbscore_copy.obs = neighbscore.obs
    neighbscore_copy.var['gene_ids'] = neighbscore_copy.var.index
    neighbscore_copy.obsm = neighbscore.obsm
    pcs.uns = neighbscore.uns
    pcs.obsm = neighbscore.obsm
    pcs.var['gene_ids'] = pcs.var.index
    pcs.obs = neighbscore.obs
    # Plot heatmaps for pc scores gained 
    # sc.pl.clustermap(sc.AnnData(pcs.to_df()[[x for x in pcs.var_names if not x.startswith('cluster')]].T),  show=True)
    # sc.pl.clustermap(sc.AnnData(pcs.to_df()[[x for x in pcs.var_names if x.startswith('cluster')]].T),  show=True)
    # pcs_heatmap_ltclust = sc.pl.clustermap(sc.AnnData(pcs.to_df()[[x for x in pcs.var_names if not x.startswith('cluster')]].T),  show=True)
    # pcs_heatmap_pathways = sc.pl.clustermap(sc.AnnData(pcs.to_df()[[x for x in pcs.var_names if x.startswith('cluster')]].T),  show=True)
    # Perform clustering over PCs and neighborhood scores
    # sc.pp.highly_variable_genes(pcs)
    sc.pp.highly_variable_genes(neighbscore_copy, n_top_genes=2000)
    # sc.pp.scale(pcs)
    sc.pp.scale(neighbscore_copy)
    # sc.tl.pca(pcs)
    sc.tl.pca(neighbscore_copy)
    # sc.pp.neighbors(pcs, n_neighbors=20, n_pcs=40)
    sc.pp.neighbors(neighbscore_copy, n_neighbors=20, n_pcs=neighbscore_copy.varm['PCs'].shape[1])
    # sc.tl.leiden(pcs, resolution = 0.3)
    sc.tl.leiden(neighbscore_copy, resolution=resolution)
    # sc.tl.umap(pcs)
    sc.tl.umap(neighbscore_copy)
    # Set cluster labels of neighborhood scores to PCs
    pcs.uns['leiden'] = neighbscore_copy.uns['leiden']
    pcs.obs['leiden'] = neighbscore_copy.obs['leiden']
    # Plot UMAP / Spatial plot
    pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_path)
    temp_fig = sc.pl.umap(neighbscore_copy, components='all', color='leiden', return_fig=True, show=False)
    temp_fig.set_size_inches(16, 16)
    pdf.savefig(temp_fig)
    plt.close()
    sc.pl.spatial(neighbscore_copy, img_key="hires", color=["leiden"], return_fig=True, size = 1.4, show=False, alpha_img=0)
    temp_fig = plt.gcf()
    temp_fig.set_size_inches(16, 16)
    pdf.savefig(temp_fig)
    plt.close()
    # Compute top N markers per cluster
    sc.tl.rank_genes_groups(neighbscore_copy, "leiden", method="wilcoxon")
    sc.pl.rank_genes_groups_heatmap(neighbscore_copy, n_genes=n_markers, groupby="leiden", show_gene_labels=True, dendrogram = False,
                                    show=False, swap_axes=True, standard_scale='var', cmap='viridis')
    temp_fig = plt.gcf()
    temp_fig.set_size_inches(20, 30)
    pdf.savefig(temp_fig)
    plt.close()
    # Plot top N pairs per cluster
    top_n = get_top_n_clust_pairs(neighbscore_copy, neighbscore_df, n=n_top)
    de_ltpairs = [pair for tup in neighbscore_copy.uns['rank_genes_groups']['names'][0:n_markers] for pair in tup]
    sc.pl.heatmap(neighbscore_copy, var_names=top_n, groupby='leiden', show=False, swap_axes=True, show_gene_labels=True, standard_scale='var', cmap='viridis')
    temp_fig = plt.gcf()
    temp_fig.set_size_inches(20, 30)
    pdf.savefig(temp_fig)
    plt.close()
    celltype.uns['leiden'] = neighbscore_copy.uns['leiden']
    celltype.obs['leiden'] = neighbscore_copy.obs['leiden']
    celltype_df = celltype.to_df().T.copy()
    celltype_df = celltype_df.apply(sum_norm, axis=0)
    celltype_df = celltype_df.fillna(0)
    celltype_df = celltype_df.apply(min_max, axis=1)
    celltype_df = celltype_df.fillna(0)
    neighbscore_df = neighbscore.to_df().T
    neighbscore_df = neighbscore_df.apply(sum_norm, axis=0)
    neighbscore_df = neighbscore_df.fillna(0)
    neighbscore_df = neighbscore_df.apply(min_max, axis=1)
    neighbscore_df = neighbscore_df.fillna(0)
    # Find cosine similarity amongst spots using celltype abundance and neighbscore of top
    # n pairs and DE ltpairs
    celltype_cosine = pd.DataFrame(data=cosine_similarity(celltype_df.T), columns=list(celltype_df.columns),
                                   index=list(celltype_df.columns))
    neighbscore_cosine = pd.DataFrame(data=cosine_similarity(neighbscore_df.T[list(set(top_n + de_ltpairs))]),
                                      columns=list(neighbscore_df.columns), index=list(neighbscore_df.columns))
    cosine = (celltype_cosine + neighbscore_cosine) / 2
    celltype_cosine = sc.AnnData(celltype_cosine)
    celltype_cosine.obs['leiden'] = neighbscore_copy.obs['leiden']
    celltype_cosine.obs['index'] = celltype_cosine.obs.index
    celltype_cosine.uns = neighbscore_copy.uns
    neighbscore_cosine = sc.AnnData(neighbscore_cosine)
    neighbscore_cosine.obs['leiden'] = neighbscore_copy.obs['leiden']
    neighbscore_cosine.obs['index'] = neighbscore_cosine.obs.index
    neighbscore_cosine.uns = neighbscore_copy.uns
    cosine = sc.AnnData(cosine)
    cosine.obs['leiden'] = neighbscore_copy.obs['leiden']
    cosine.obs['index'] = cosine.obs.index
    cosine.uns = neighbscore_copy.uns
    sc.pl.heatmap(celltype_cosine, celltype_cosine.obs.groupby('leiden').index.apply(list).to_dict(), groupby='leiden',
                  show=False, cmap='YlOrBr')
    temp_fig = plt.gcf()
    temp_fig.set_size_inches(16, 16)
    pdf.savefig(temp_fig)
    plt.close()
    sc.pl.heatmap(neighbscore_cosine, neighbscore_cosine.obs.groupby('leiden').index.apply(list).to_dict(),
                  groupby='leiden', show=False, cmap='YlOrBr')
    temp_fig = plt.gcf()
    temp_fig.set_size_inches(16, 16)
    pdf.savefig(temp_fig)
    plt.close()
    sc.pl.heatmap(cosine, cosine.obs.groupby('leiden').index.apply(list).to_dict(), groupby='leiden', show=False, cmap='YlOrBr')
    temp_fig = plt.gcf()
    temp_fig.set_size_inches(16, 16)
    pdf.savefig(temp_fig)
    plt.close()
    celltype_df = celltype.to_df().T.copy()
    #celltype_df = celltype_df.apply(sum_norm, axis=0)
    #celltype_df = celltype_df.fillna(0)
    celltype_df = celltype_df.apply(min_max, axis=1)
    celltype_df = celltype_df.fillna(0)
    celltype_copy = sc.AnnData(celltype_df.T)
    celltype_copy.obs['leiden'] = celltype.obs['leiden']
    celltype_copy.uns = celltype.uns
    # hierarchically cluster each group
    clusterobj = sc.pl.clustermap(celltype_copy, show=False)
    plt.close()
    reordered_var = list(celltype_copy.var_names[clusterobj.dendrogram_col.reordered_ind])
    celltype_copy = celltype_copy[:,reordered_var]
    clustered_spot_order = []
    for cluster in celltype_copy.obs.leiden.unique():
        clusterobj = sc.pl.clustermap(celltype_copy[celltype_copy.obs.leiden[celltype_copy.obs.leiden==cluster].index,:], show=False, col_cluster=False)
        plt.close()
        clustered_spot_order+=list(celltype_copy[celltype_copy.obs.leiden[celltype_copy.obs.leiden==cluster].index,:].obs_names[clusterobj.dendrogram_row.reordered_ind])
    celltype_copy= celltype_copy[clustered_spot_order,:]
    sc.pl.heatmap(celltype_copy, var_names =celltype_copy.var_names, groupby='leiden', show=False, swap_axes=True, standard_scale='var', show_gene_labels=True, cmap='bone_r')
    temp_fig = plt.gcf()
    temp_fig.set_size_inches(30, 25)
    pdf.savefig(temp_fig)
    plt.close()
    clustered_spot_order = []
    for cluster in pcs.obs.leiden.unique():
        clusterobj = sc.pl.clustermap(pcs[pcs.obs.leiden[pcs.obs.leiden==cluster].index,:], show=False, col_cluster=False)
        clustered_spot_order+=list(pcs[pcs.obs.leiden[pcs.obs.leiden==cluster].index,:].obs_names[clusterobj.dendrogram_row.reordered_ind])
    pcs= pcs[clustered_spot_order,:]
    if len([x for x in pcs.var_names if not x.startswith('cluster')])>2:
        clusterobj = sc.pl.clustermap(pcs[:,[x for x in pcs.var_names if not x.startswith('cluster')]], show=False)
        plt.close()
        reordered_var = list(pcs[:,[x for x in pcs.var_names if not x.startswith('cluster')]].var_names[clusterobj.dendrogram_col.reordered_ind])
        sc.pl.heatmap(pcs, var_names=reordered_var, groupby='leiden',
                      show=False, swap_axes=True, standard_scale='var', show_gene_labels=True, cmap='YlGnBu_r')
        temp_fig = plt.gcf()
        temp_fig.set_size_inches(40, 20)
        pdf.savefig(temp_fig)
        plt.close()
    if len([x for x in pcs.var_names if x.startswith('cluster')])>2:
        clusterobj = sc.pl.clustermap(pcs[:,[x for x in pcs.var_names if x.startswith('cluster')]], show=False)
        plt.close()
        reordered_var = list(pcs[:,[x for x in pcs.var_names if x.startswith('cluster')]].var_names[clusterobj.dendrogram_col.reordered_ind])
        sc.pl.heatmap(pcs, var_names=reordered_var, groupby='leiden', show=False, swap_axes=True, standard_scale='var', show_gene_labels=True, cmap='YlGnBu_r')
        temp_fig = plt.gcf()
        temp_fig.set_size_inches(40, 20)
        pdf.savefig(temp_fig)
        plt.close()
    pdf.close()


def sankeyPlot(neighbscore, celltype, ltpairs, n_celltype=5, clusters='All', title=None, labelsize=2, labelcolor='#000000'):
    """Summarizes ligand-target activity across celltypes and domains using a Sankey plot for a given set of ligand-target pairs.

    :param neighbscore: adata object of neighborhood scores with leiden clusters (obs column 'leiden' of type categorical)
    :type neighbscore: AnnData
    :param celltype: adata object with celltype distributions (proportiond or abundance) across spots
    :type celltype: AnnData
    :param ltpairs: set of ligand-target pairs to display
    :type ltpairs: list
    :param n_celltype: top n highly expressed celltypes to be considered per domain, defaults to 5
    :type n_celltype: int, optional
    :param clusters: list of domains to consider, defaults to 'All'
    :type clusters: list, optional
    :param title: Title of the sankey plot, defaults to None
    :type title: str, optional
    :param labelsize: size of labels in the plot, defaults to 2
    :type labelsize: int, optional
    :param labelcolor: color of labels in the plot, defaults to '#000000'
    :type labelcolor: str, optional
    """
    celltype_df = celltype.to_df().T
    celltype_df = celltype_df.apply(sum_norm, axis=0)
    celltype_copy = sc.AnnData(celltype_df.T)
    celltype_copy.obs['leiden'] = neighbscore.obs['leiden']
    celltype_copy.uns = neighbscore.uns
    if clusters=='All':
        clusters = list(celltype_copy.obs.leiden.cat.categories)
    if not isinstance(clusters, list):
        raise "ERROR: clusters should either be 'All' or list of cluster ids of type string"
    # Get average celltype score and average neighbscore
    celltype_avg = pd.DataFrame(columns=celltype_copy.var_names, index=clusters)
    for clust in clusters:
        celltype_avg.loc[clust] = celltype_copy[celltype_copy.obs['leiden'].isin([clust]), :].X.mean(0)
    celltype_avg = celltype_avg.T.to_dict('dict')
    ltpair_avg = pd.DataFrame(columns=neighbscore.var_names, index=clusters)
    for clust in clusters:
        ltpair_avg.loc[clust] = neighbscore.raw[neighbscore.obs['leiden'].isin([clust]), :].X.mean(0)
    ltpair_avg = ltpair_avg[ltpairs].T.to_dict('dict')
    # Create a dictionary of links and nodes
    ligands = []
    targets = []
    links = []
    clust_v_ligand = {}
    clust_v_ligand_sum = {}
    target_v_clust = {}
    target_v_clust_sum = {}
    lt_sum = {}
    temp_sum = {}
    for pair in ltpairs:
        ligand = pair.split(':')[0]
        target = pair.split(':')[1]
        ligands += ligand
        targets += target
        temp_sum[pair] = 0
        for clust in ltpair_avg.keys():
            if clust in clusters:
                temp_sum[pair] += ltpair_avg[clust][pair]
                lig_clust = clust + ':' + ligand
                tar_clust = clust + ':' + target
                if target in lt_sum.keys():
                    lt_sum[target] += ltpair_avg[clust][pair]
                else:
                    lt_sum[target] = ltpair_avg[clust][pair]
                if ligand in clust_v_ligand_sum.keys():
                    clust_v_ligand_sum[ligand] += ltpair_avg[clust][pair]
                else:
                    clust_v_ligand_sum[ligand] = ltpair_avg[clust][pair]
                if lig_clust in clust_v_ligand.keys():
                    clust_v_ligand[lig_clust] += ltpair_avg[clust][pair]
                else:
                    clust_v_ligand[lig_clust] = ltpair_avg[clust][pair]
                if target in target_v_clust_sum.keys():
                    target_v_clust_sum[target] += ltpair_avg[clust][pair]
                else:
                    target_v_clust_sum[target] = ltpair_avg[clust][pair]
                if tar_clust in target_v_clust.keys():
                    target_v_clust[tar_clust] += ltpair_avg[clust][pair]
                else:
                    target_v_clust[tar_clust] = ltpair_avg[clust][pair]

    value_temp = []
    for lig_clust in clust_v_ligand.keys():
        value_temp.append(clust_v_ligand[lig_clust]) # / clust_v_ligand_sum[lig_clust.split(':')[1]]
        links.append({'source': lig_clust.split(':')[0], 'target': lig_clust.split(':')[1],
                      'value': clust_v_ligand[lig_clust]}) #  / clust_v_ligand_sum[lig_clust.split(':')[1]]

    # max_temp = float(max(value_temp))
    clust_v_ligand_avg = sum(value_temp) / len(value_temp)

    value_temp = []
    for tar_clust in target_v_clust.keys():
        value_temp.append(target_v_clust[tar_clust]) #  / target_v_clust_sum[tar_clust.split(':')[1]]
    
    target_v_clust_avg = sum(value_temp)/len(value_temp)
    
    for tar_clust in target_v_clust.keys():
        links.append({'source': tar_clust.split(':')[1], 'target': tar_clust.split(':')[0]+'_',
                      'value': (target_v_clust[tar_clust])}) #  / target_v_clust_sum[tar_clust.split(':')[1]])*(clust_v_ligand_avg/target_v_clust_avg

    value_temp = []
    for pair in ltpairs:
        value_temp.append(temp_sum[pair]) #/ len(ltpair_avg.keys()))

    ligand_target_avg = sum(value_temp) / len(value_temp)

    for pair in ltpairs:
        ligand = pair.split(':')[0]
        target = pair.split(':')[1]
        links.append({'source': ligand, 'target': target, 'value': (temp_sum[pair])})  # lt_sum[target])*max_temp})   / len(ltpair_avg.keys())) * (clust_v_ligand_avg / ligand_target_avg

    # Get top n celltypes per cluster
    top_celltypes = []
    value_temp = []
    for clust in celltype_avg.keys():
        celltypes_temp = sorted(celltype_avg[clust], key=celltype_avg[clust].get, reverse=True)[:n_celltype]
        for ct in celltypes_temp:
            value_temp.append(celltype_avg[clust][ct])
    
    ct_avg = sum(value_temp)/len(value_temp)
    for clust in celltype_avg.keys():
        celltypes_temp = sorted(celltype_avg[clust], key=celltype_avg[clust].get, reverse=True)[:n_celltype]
        top_celltypes += celltypes_temp
        for ct in celltypes_temp:
            links.append({'source': ct, 'target': clust, 'value': celltype_avg[clust][ct] * (
                    clust_v_ligand_avg / ct_avg)})
    top_celltypes = list(set(top_celltypes))
    links = pd.DataFrame(links)
    links['source'] = links['source'].astype(str)
    links['target'] = links['target'].astype(str)
    links['value'] = links['value'].astype(float)
    nodes = list(set(links['source'].tolist() + links['target'].tolist()))
    # get node ids and link colors
    colors = {}
    for node in nodes:
        colors[node] = "#" + ''.join([random.choice('0123456789ABCDE') for j in range(6)])
    links['color'] = links['source']
    links['color'] = links['color'].map(colors)
    node_id = {}
    id = 0
    for node in nodes:
        node_id[node] = id
        id += 1
    links['source'] = links['source'].map(node_id)
    links['target'] = links['target'].map(node_id)
    # Plot Sankey plot
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=nodes,
            color=[colors[node] for node in nodes]
        ),
        link=dict(
            source=links['source'].tolist(),
            target=links['target'].tolist(),
            value=links['value'].tolist()
            #color='rgba(140,140,140,0.5)'
        ),
        textfont=dict(
            color=labelcolor,
            size=labelsize
        )
    )])
    if title is not None:
        fig.update_layout(title_text=title, font_size=10)
    return fig

def pcs_v_neighbscore(neighbscore, ltpair_clusters=None, pdf_path=None, spatialfeatureplot=True, clustermap=True, size=1.4, colormap = sns.color_palette("Spectral",as_cmap=True)):
    """Generate heatmaps of pathway/geneset scores vs individual liagnd-target contribution

    :param neighbscore: adata object of neighborhood scores
    :type neighbscore: AnnData
    :param ltpair_clusters: dictionary of genesets/pathways with corresponding ligand-target pairs (can be generated using create_cluster), defaults to None
    :type ltpair_clusters: dict, optional
    :param pdf_path: path to save plots as a pdf file, defaults to None
    :type pdf_path: str, optional
    :param spatialfeatureplot: Include spatial feature plots of geneset/pathway activity, defaults to True
    :type spatialfeatureplot: bool, optional
    :param clustermap: Generate clustermap of each pathway/geneset, defaults to True
    :type clustermap: bool, optional
    :param size: spot size (required if spatialfeatureplot is set to True), defaults to 1.4
    :type size: float, optional
    :param colormap: colormap for the heatmap and spatial fetaure plots, defaults to sns.color_palette("Spectral",as_cmap=True)
    :type colormap: _type_, optional
    """
    neighbscore_copy = neighbscore
    sc.pp.filter_genes(neighbscore_copy, min_cells=1)
    neighbscore_df = neighbscore_copy.to_df().T
    neighbscore_df = neighbscore_df.apply(min_max, axis=1)
    if ltpair_clusters is None:
        raise 'ERROR: ltpair_clusters arg is None. Provide cluster object using create_cluster'
    pcs = {}
    for cluster, ltpairs in ltpair_clusters.items():
        pca = TruncatedSVD()
        pcs[cluster] = list(pca.fit_transform(neighbscore_df.loc[ltpairs,].T)[:, 0])
    pcs = pd.DataFrame.from_dict(pcs, orient='index')
    pcs.columns = neighbscore.obs_names
    pcs = pcs.applymap(abs)
    pcs = pcs.apply(min_max, axis=1)
    pcs.fillna(0, inplace=True)
    pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_path)
    for cluster, ltpairs in ltpair_clusters.items():
        print(cluster)
        pc_v_neighb = pd.DataFrame(data=pcs.loc[cluster,].tolist(), index=neighbscore.obs_names, columns=[cluster]).T
        pc_v_neighb = pd.concat([pc_v_neighb, neighbscore_df.loc[ltpairs,]])
        pc_v_neighb_df = pc_v_neighb
        pc_v_neighb = sc.AnnData(pc_v_neighb.T)
        pc_v_neighb.obsm = neighbscore.obsm
        pc_v_neighb.obs = neighbscore.obs
        pc_v_neighb.uns = neighbscore.uns
        if clustermap:
            clusterobj = sc.pl.clustermap(sc.AnnData(pc_v_neighb_df), show=False)
            plt.close()
            reordered_rows = list(pc_v_neighb_df.index[clusterobj.dendrogram_row.reordered_ind])
            reordered_rows.remove(cluster)
            reordered_rows = [cluster]+reordered_rows
            sc.pl.clustermap(sc.AnnData(pc_v_neighb_df.loc[reordered_rows,:]), show=False, row_cluster=False)
            temp_fig = plt.gcf()
            temp_fig.set_size_inches(16, 16)
            pdf.savefig(temp_fig)
            plt.close()
        if spatialfeatureplot:
            sc.pl.spatial(pc_v_neighb, color = cluster, show = False, size=size, alpha_img=0)
            temp_fig = plt.gcf()
            temp_fig.set_size_inches(16, 16)
            pdf.savefig(temp_fig)
            plt.close()
    pdf.close()    

def ligand_ranking(neighbscore, celltype, scrna, ligand_receptor_pairs, ligand_target_regulatory_potential, domain, receptor_exp=0.1, markers={'top':100}, domain_celltypes=['top',5], celltype_colors={'auto':True}):
    """_summary_

    :param neighbscore: _description_
    :type neighbscore: _type_
    :param celltype: _description_
    :type celltype: _type_
    :param scrna: _description_
    :type scrna: _type_
    :param ligand_receptor_pairs: _description_
    :type ligand_receptor_pairs: _type_
    :param ligand_target_regulatory_potential: _description_
    :type ligand_target_regulatory_potential: _type_
    :param domain: _description_
    :type domain: _type_
    :param receptor_exp: _description_, defaults to 0.1
    :type receptor_exp: float, optional
    :param markers: _description_, defaults to {'top':100}
    :type markers: dict, optional
    :param domain_celltypes: _description_, defaults to ['top',5]
    :type domain_celltypes: list, optional
    :param celltype_colors: _description_, defaults to {'auto':True}
    :type celltype_colors: dict, optional
    :return: _description_
    :rtype: _type_
    """
    sc.pp.filter_genes(neighbscore, min_cells=1)
    neighbscore_df = neighbscore.raw.to_adata().to_df()
    celltype_df = celltype.to_df()
    celltype_df = celltype_df.apply(sum_norm, axis=1).T
    #Get top n celltypes
    
    celltype_copy = sc.AnnData(celltype_df.T)
    celltype_copy.obs['leiden'] = neighbscore.obs['leiden']
    celltype_copy.uns = neighbscore.uns
    celltype_copy = celltype_copy[celltype_copy.obs.loc[celltype_copy.obs.leiden == domain].index].copy()
    celltype_df_copy = celltype_copy.to_df()
    # Get average celltype score and average neighbscore
    celltype_avg = pd.DataFrame(columns=celltype_copy.var_names, index=celltype_copy.obs['leiden'].cat.categories)
    celltype_avg.loc[domain] = celltype_copy[celltype_copy.obs['leiden'].isin([domain]), :].X.mean(0)
    celltype_avg = celltype_avg.T.to_dict('dict')
    if domain_celltypes[0]=='top':
        top_celltypes = sorted(celltype_avg[domain], key=celltype_avg[domain].get, reverse=True)[:domain_celltypes[1]]
        top_celltypes = list(set(top_celltypes))
    else:
        top_celltypes = domain_celltypes
    # Get celltype markers
    if 'top' in markers.keys():
        sc.pp.filter_genes(scrna, min_cells=1)
        scrna = scrna[scrna.obs.loc[scrna.obs.celltype.isin(top_celltypes)].index].copy()
        scrna.obs.celltype = scrna.obs.celltype.astype("category")
        scrna.raw = scrna
        sc.tl.rank_genes_groups(scrna, "celltype", method="wilcoxon")
        markers = {}
        for col in scrna.uns['rank_genes_groups']['names'][0].dtype.names:
            markers[col] = []

        for index in range(len(scrna.uns['rank_genes_groups']['names'])):
            for col in range(len(scrna.uns['rank_genes_groups']['names'][index])):
                if scrna.uns['rank_genes_groups']['logfoldchanges'][index][col] > 0 and scrna.uns['rank_genes_groups']['pvals_adj'][index][col] < 0.05:
                    markers[scrna.uns['rank_genes_groups']['names'][0].dtype.names[col]].append(scrna.uns['rank_genes_groups']['names'][index][col])
        markers = pd.DataFrame(list(markers.values()), index=markers.keys()).T
        markers = markers.iloc[0:200,:]
    else:
        markers = pd.DataFrame(list(markers.values()), index=markers.keys()).T
    markers = markers[top_celltypes]
    # Select ligands and targets that are markers and receptors that are expressed
    receptors = {}
    for receptor in ligand_receptor_pairs.receptor.unique():
        if receptor in scrna.var_names:
            receptors[receptor] = []
            for ct in markers.columns:
                temp = scrna[scrna.obs.loc[scrna.obs.celltype == ct].index, receptor].to_df()
                if  np.count_nonzero(temp)/len(temp) >= receptor_exp:
                    receptors[receptor].append(ct)
    
    ligands, targets = zip(*(s.split(":") for s in neighbscore.var_names))
    ct_spec_mark = {'ligand':{},'target':{}}
    for i in markers.columns:
        ct_spec_mark['target'][i] = []

    for index in range(len(ligands)):
        ligand = ligands[index]
        target = targets[index]
        ligand_celltypes = list(markers.columns[markers.isin([ligand]).any()])
        target_celltypes = list(markers.columns[markers.isin([target]).any()])
        if len(ligand_celltypes)>=1:
            ct_spec_mark['ligand'][ligand] = []
            for ctl in ligand_celltypes:
                ct_spec_mark['ligand'][ligand].append(ctl)
        if len(target_celltypes)>=1:
            for ctt in target_celltypes:
                if target not in ct_spec_mark['target'][ctt]:
                    ct_spec_mark['target'][ctt].append(target)
    
    # For each domain find score for each ligand and target and its corresponding celltype
    lrt_dict = {}
    for i,row in ligand_receptor_pairs.iterrows():
        if row['ligand'] in ct_spec_mark['ligand'].keys() and row['receptor'] in receptors.keys():
            if row['ligand'] not in lrt_dict.keys():
                lrt_dict[row['ligand']] = {}
                lrt_dict[row['ligand']]['celltype'] = ct_spec_mark['ligand'][row['ligand']]
            for ct in receptors[row['receptor']]:
                for target in ct_spec_mark['target'][ct]:
                    if row['ligand']+':'+target in neighbscore.var_names:
                        if target not in lrt_dict[row['ligand']].keys():
                            lrt_dict[row['ligand']][target] = []
                        if ct not in lrt_dict[row['ligand']][target]:
                            lrt_dict[row['ligand']][target].append(ct)
    
    ligand_score = {}
    for ligand in lrt_dict.keys():
        ligand_score[ligand] = {}
        targets = list(lrt_dict[ligand].keys())
        ligand_ct = lrt_dict[ligand]['celltype']
        targets.remove('celltype')
        for target in targets:
            target_ct = lrt_dict[ligand][target]
            temp = celltype_df_copy.loc[:,list(set(ligand_ct+target_ct))]
            spots = temp[(temp.select_dtypes(include=['number']) != 0).all(1)].index
            ligand_score[ligand][target] = neighbscore_df.loc[spots,ligand+':'+target].mean()
    
    targets = []
    for ligand in ligand_score.keys():
        for target in ligand_score[ligand].keys():
            if target not in targets:
                targets.append(target)
    
    # Get ligand ranking score
    ligand_ranking = {}
    for ligand in ligand_score.keys():
        if len(ligand_score[ligand].keys()) > 0:
            ligand_score_temp = []
            ligand_regulatory_potential = []
            for target in targets:
                if target in ligand_score[ligand].keys():
                    ligand_score_temp.append(ligand_score[ligand][target])
                else:
                    ligand_score_temp.append(0)
                if target in ligand_target_regulatory_potential[ligand].keys():
                    ligand_regulatory_potential.append(ligand_target_regulatory_potential[ligand][target])
                else:
                    ligand_regulatory_potential.append(0)
            ligand_ranking[ligand] = pd.Series(ligand_score_temp).corr(pd.Series(ligand_regulatory_potential))
    
    # Plot avg neighborhood scores at colocalized spots ranked by ligand ranking score
    ligand_score_df = pd.DataFrame(ligand_score).fillna(0).T
    ligand_score_df =(ligand_score_df-ligand_score_df.min())/(ligand_score_df.max()-ligand_score_df.min())
    ligand_score_df = ligand_score_df.fillna(0)
    ligand_score_df = ligand_score_df.loc[sorted(ligand_ranking, key=lambda k:ligand_ranking[k], reverse=True),:]
    ligand_score_df = ligand_score_df.loc[~(ligand_score_df==0).all(axis=1)]
    ligand_score_df = ligand_score_df.loc[:, (ligand_score_df != 0).any(axis=0)]
    # Sort columns
    col_order = []
    for index in ligand_score_df.index:
        for target in ligand_score[index].keys():
            if target not in col_order and target in ligand_score_df.columns:
                col_order.append(target)
    ligand_score_df = ligand_score_df[col_order]
    ligand_score_df.insert(loc=0, column='ligand score', value=[ligand_ranking[x] for x in ligand_score_df.index])

    if 'auto' in celltype_colors.keys():
        colors=dict(zip(top_celltypes,distinctipy.get_colors(len(top_celltypes))))
    else:
        colors = {k:celltype_colors[k] for k in top_celltypes}
    
    #Create dummy joint plot
    sns.set(rc={'axes.facecolor':'#ffffff', 'figure.facecolor':'#ffffff'}, font_scale=1.7)
    g = sns.jointplot(data=ligand_score_df, x=ligand_score_df.iloc[:,1], y=ligand_score_df.iloc[:,1], kind='hist', bins=(len(ligand_score_df.columns), len(ligand_score_df.index)))
    g.ax_marg_y.cla()
    g.ax_marg_x.cla()

    #Generate heatmap
    mask = ligand_score_df.apply(lambda x: x if not x.name.endswith('ligand score') else 0,result_type='broadcast',axis=0).eq(0)
    sns.heatmap(data=ligand_score_df, ax=g.ax_joint, mask=mask, cbar=True, cmap='BuPu', cbar_kws = dict(use_gridspec=False,location="right", shrink=0.2, pad=0.1))
    plt.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1)
    # get the current positions of the joint ax and the ax for the marginal x
    pos_joint_ax = g.ax_joint.get_position()
    pos_marg_x_ax = g.ax_marg_x.get_position()
    # reposition the joint ax so it has the same width as the marginal x ax
    g.ax_joint.set_position([pos_joint_ax.x0, pos_joint_ax.y0, pos_marg_x_ax.width, pos_joint_ax.height])
    # reposition the colorbar using new x positions and y positions of the joint ax
    g.fig.axes[-1].set_position([.83, pos_joint_ax.y0, .07, pos_joint_ax.height])
    mask = ligand_score_df.apply(lambda x: x if x.name.endswith('ligand score') else 0,result_type='broadcast',axis=0).eq(0)
    sns.heatmap(data=ligand_score_df, ax=g.ax_joint, mask=mask, cbar=True, cmap='RdYlBu_r', cbar_kws = dict(use_gridspec=False,location="left", shrink=0.15, pad=0.2))
    plt.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1)
    # get the current positions of the joint ax and the ax for the marginal x
    pos_joint_ax = g.ax_joint.get_position()
    pos_marg_x_ax = g.ax_marg_x.get_position()
    # reposition the joint ax so it has the same width as the marginal x ax
    g.ax_joint.set_position([pos_joint_ax.x0, pos_joint_ax.y0, pos_marg_x_ax.width, pos_joint_ax.height])
    # reposition the colorbar using new x positions and y positions of the joint ax
    g.fig.axes[-1].set_position([.83, pos_joint_ax.y0, .07, pos_joint_ax.height])

    #celltype bar plot for ligands
    ct_bar_plot = pd.DataFrame(0, columns = top_celltypes, index = ligand_score_df.index)
    for col in ct_bar_plot.columns:
        for row in ct_bar_plot.index:
            if col in ct_spec_mark['ligand'][row]:
                ct_bar_plot.loc[row,col] = celltype_avg[domain][col]
    ct_bar_plot = ct_bar_plot.apply(sum_norm, axis=1)
    ct_bar_plot = ct_bar_plot.cumsum(axis=1)
    cols = list(ct_bar_plot.columns)
    cols.reverse()
    for col in cols:
        g.ax_marg_y.barh(np.arange(0.5, len(ligand_score_df.index)), ct_bar_plot[col].tolist(), color=colors[col])

    #celltype bar plot for targets
    ct_bar_plot = pd.DataFrame(0, columns = top_celltypes, index = list(ligand_score_df.columns)[1:])
    for col in ct_bar_plot.columns:
        for row in ct_bar_plot.index:
            if row in ct_spec_mark['target'][col]:
                ct_bar_plot.loc[row,col] = celltype_avg[domain][col]
    ct_bar_plot = ct_bar_plot.apply(sum_norm, axis=1)
    ct_bar_plot = ct_bar_plot.cumsum(axis=1)
    cols = list(ct_bar_plot.columns)
    cols.reverse()
    for col in cols:
        g.ax_marg_x.bar(np.arange(1.5, len(ligand_score_df.columns)), ct_bar_plot[col].tolist(), color=colors[col])

    g.ax_joint.set_xticks(np.arange(0.5, len(ligand_score_df.columns)))
    g.ax_joint.set_xticklabels(list(ligand_score_df.columns))
    g.ax_joint.set_yticks(np.arange(0.5, len(ligand_score_df.index)))
    g.ax_joint.set_yticklabels(list(ligand_score_df.index))

    # remove ticks between heatmap and histograms
    g.ax_marg_x.tick_params(axis='x', bottom=False, labelbottom=False)
    g.ax_marg_y.tick_params(axis='y', left=False, labelleft=False)
    # remove ticks showing the heights of the histograms
    g.ax_marg_x.tick_params(axis='y', left=False, labelleft=False)
    g.ax_marg_y.tick_params(axis='x', bottom=False, labelbottom=False)

    #Create legend
    g.fig.suptitle("Cluster "+domain)
    handles = [Patch(facecolor=colors[name]) for name in colors]
    plt.legend(handles, colors, bbox_to_anchor = (1,1), bbox_transform = plt.gcf().transFigure, ncol=3, loc='upper right')
    fig = plt.gcf()
    plt.close()
    return fig
