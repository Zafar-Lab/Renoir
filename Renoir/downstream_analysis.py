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


def sankeyPlot(neighbscore, celltype, ltpairs, n_celltype=5, clusters='All', title=None, path=None, labelsize=2, labelcolor='#000000'):
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
    if path is None:
        fig.show()
    else:
        fig.write_html(path)

def pcs_v_neighbscore(neighbscore, ltpair_clusters=None, pdf_path=None, spatialfeatureplot=True, clustermap=True, size=1.4, colormap = sns.color_palette("Spectral",as_cmap=True)):
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
