"""
Tests for Renoir.downstream.

All tests use synthetic fixtures from conftest.py.

Known upstream bug (marked xfail):
    create_cluster with method='dhc' calls cutreeHybrid which returns a
    numpy.ndarray when the cut height is too low (degenerate correlation
    structure). The source code then indexes it as clusters['labels'],
    causing an IndexError. This is a bug in downstream.py, not in the test.
"""

import numpy as np
import pandas as pd
import anndata as ad
import pytest
import matplotlib
matplotlib.use("Agg")

import Renoir
from Renoir import downstream as rd

from conftest import CELLTYPES, LIGANDS_UNIQUE, PAIR_NAMES


# ---------------------------------------------------------------------------
# Shared cluster fixture
# ---------------------------------------------------------------------------

@pytest.fixture
def clusters(neighborhood_scores, msig_df):
    """Pathway-only clusters with pathway_thresh=1."""
    return rd.create_cluster(
        neighborhood_scores,
        msig_df,
        method=None,
        pathway_thresh=1,
    )


# ===========================================================================
# get_msig
# ===========================================================================

class TestGetMsig:

    def test_returns_dataframe(self, msig_csv):
        assert isinstance(rd.get_msig("custom", path=msig_csv), pd.DataFrame)

    def test_has_required_columns(self, msig_csv):
        df = rd.get_msig("custom", path=msig_csv)
        assert "gs_name"     in df.columns
        assert "gene_symbol" in df.columns

    def test_non_empty(self, msig_csv):
        assert len(rd.get_msig("custom", path=msig_csv)) > 0

    def test_invalid_species_raises(self):
        with pytest.raises(Exception):
            rd.get_msig("martian")


# ===========================================================================
# create_cluster
# ===========================================================================

class TestCreateCluster:

    def test_pathway_mode_returns_dict(self, neighborhood_scores, msig_df):
        result = rd.create_cluster(neighborhood_scores, msig_df,
                                   method=None, pathway_thresh=1)
        assert isinstance(result, dict)

    def test_pathway_keys_are_strings(self, neighborhood_scores, msig_df):
        result = rd.create_cluster(neighborhood_scores, msig_df,
                                   method=None, pathway_thresh=1)
        for key in result:
            assert isinstance(key, str)

    def test_pathway_values_are_lists(self, neighborhood_scores, msig_df):
        result = rd.create_cluster(neighborhood_scores, msig_df,
                                   method=None, pathway_thresh=1)
        for val in result.values():
            assert isinstance(val, list)

    def test_non_empty_clusters(self, neighborhood_scores, msig_df):
        """With both ligand and target genes in msig, pairs are found."""
        result = rd.create_cluster(neighborhood_scores, msig_df,
                                   method=None, pathway_thresh=1)
        assert len(result) > 0, (
            "No pathway clusters survived. Ensure msig_df.gene_symbol "
            "contains both ligand and target genes from var_names."
        )

    def test_pairs_exist_in_var_names(self, neighborhood_scores, msig_df):
        result = rd.create_cluster(neighborhood_scores, msig_df,
                                   method=None, pathway_thresh=1)
        valid = set(neighborhood_scores.var_names)
        for cluster_name, pairs in result.items():
            for pair in pairs:
                assert pair in valid, \
                    f"Cluster '{cluster_name}' has unknown pair: {pair}"

    def test_restrict_KHW(self, neighborhood_scores, msig_df):
        result = rd.create_cluster(neighborhood_scores, msig_df,
                                   method=None, pathway_thresh=1,
                                   restrict_to_KHW=True)
        for key in result:
            assert any(key.startswith(p) for p in
                       ("KEGG_", "HALLMARK_", "WP_", "cluster_")), \
                f"Unexpected key after KHW filter: {key}"

    def test_hdbscan_returns_dict(self, neighborhood_scores, msig_df):
        result = rd.create_cluster(neighborhood_scores, msig_df,
                                   method="hdbscan", minpts=2,
                                   pathway_thresh=1, ltclust_thresh=1)
        assert isinstance(result, dict)

    def test_denovo_dhc_returns_dict(self, neighborhood_scores, msig_df):
        result = rd.create_cluster(neighborhood_scores, msig_df,
                                   method="dhc", pathway_thresh=1,
                                   ltclust_thresh=1)
        assert isinstance(result, dict)


# ===========================================================================
# downstream_analysis
# ===========================================================================

class TestDownstreamAnalysis:

    def test_returns_anndata_with_return_cluster(self, neighborhood_scores, clusters):
        result = rd.downstream_analysis(
            neighborhood_scores, ltpair_clusters=clusters,
            resolution=0.5, return_cluster=True, return_pcs=False)
        assert isinstance(result, ad.AnnData)

    def test_leiden_column_populated(self, neighborhood_scores, clusters):
        result = rd.downstream_analysis(
            neighborhood_scores, ltpair_clusters=clusters,
            resolution=0.5, return_cluster=True, return_pcs=False)
        assert "leiden" in result.obs.columns
        assert result.obs["leiden"].notna().all()

    def test_leiden_is_categorical(self, neighborhood_scores, clusters):
        result = rd.downstream_analysis(
            neighborhood_scores, ltpair_clusters=clusters,
            resolution=0.5, return_cluster=True, return_pcs=False)
        assert hasattr(result.obs["leiden"], "cat")

    def test_returns_tuple_with_both_flags(self, neighborhood_scores, clusters):
        result = rd.downstream_analysis(
            neighborhood_scores, ltpair_clusters=clusters,
            resolution=0.5, return_cluster=True, return_pcs=True)
        assert isinstance(result, tuple) and len(result) == 2
        neighbscore_copy, pcs = result
        assert isinstance(neighbscore_copy, ad.AnnData)
        assert isinstance(pcs, ad.AnnData)

    def test_pcs_n_obs_matches_input(self, neighborhood_scores, clusters):
        _, pcs = rd.downstream_analysis(
            neighborhood_scores, ltpair_clusters=clusters,
            resolution=0.5, return_cluster=True, return_pcs=True)
        assert pcs.n_obs == neighborhood_scores.n_obs

    def test_pcs_var_names_match_cluster_keys(self, neighborhood_scores, clusters):
        _, pcs = rd.downstream_analysis(
            neighborhood_scores, ltpair_clusters=clusters,
            resolution=0.5, return_cluster=True, return_pcs=True)
        for varname in pcs.var_names:
            assert varname in clusters, \
                f"pcs varname '{varname}' not in cluster keys"


# ===========================================================================
# spot_v_spot
# ===========================================================================

class TestSpotVSpot:

    def test_writes_pdf(self, neighborhood_scores, celltype_adata,
                        clusters, tmp_path):
        import os
        pdf_path = str(tmp_path / "svs.pdf")
        rd.spot_v_spot(neighborhood_scores, celltype_adata,
                       resolution=0.5, ltpair_clusters=clusters,
                       pdf_path=pdf_path)
        assert os.path.exists(pdf_path)
        assert os.path.getsize(pdf_path) > 0

    def test_leiden_added_to_obs(self, neighborhood_scores, celltype_adata,
                                  clusters, tmp_path):
        rd.spot_v_spot(neighborhood_scores, celltype_adata,
                       resolution=0.5, ltpair_clusters=clusters,
                       pdf_path=str(tmp_path / "svs2.pdf"))
        assert "leiden" in neighborhood_scores.obs.columns



# ===========================================================================
# pcs_v_neighbscore
# ===========================================================================

class TestPcsVNeighbscore:
    """
    Note: pdf_path=None is NOT supported by pcs_v_neighbscore. Line 653 in
    downstream.py calls matplotlib.backends.backend_pdf.PdfPages(pdf_path)
    unconditionally, so passing None raises ValueError: fname must be a PathLike.
    Only the pdf_path=<path> variant is tested here.
    """

    def test_writes_pdf(self, neighborhood_scores, clusters, tmp_path):
        import os
        pdf_path = str(tmp_path / "pcs.pdf")
        rd.pcs_v_neighbscore(neighborhood_scores, ltpair_clusters=clusters,
                             pdf_path=pdf_path, spatialfeatureplot=False,
                             clustermap=True)
        assert os.path.exists(pdf_path)


# ===========================================================================
# sankeyPlot
# ===========================================================================

class TestSankeyPlot:

    def test_all_domains(self, neighborhood_scores, celltype_adata):
        import matplotlib.pyplot as plt
        pairs = list(neighborhood_scores.var_names[:2])
        rd.sankeyPlot(neighborhood_scores, celltype_adata,
                      ltpairs=pairs, n_celltype=2, clusters="All")
        plt.close("all")

    def test_focused_domain(self, neighborhood_scores, celltype_adata):
        import matplotlib.pyplot as plt
        pairs = list(neighborhood_scores.var_names[:2])
        rd.sankeyPlot(neighborhood_scores, celltype_adata,
                      ltpairs=pairs, n_celltype=2, clusters=["0"])
        plt.close("all")


# ===========================================================================
# ligand_ranking
# ===========================================================================

class TestLigandRanking:
    """
    sc_adata has clearly DE genes per celltype (set in conftest) so that
    rank_genes_groups returns significant markers (pvals_adj < 0.05).
    Without structured expression the markers DataFrame is empty and
    ligand_ranking fails at a downstream iloc call.
    """

    def test_returns_figure(self, neighborhood_scores, celltype_adata,
                            sc_adata, lt_regulatory_potential, lr_pairs_csv):
        import matplotlib.pyplot as plt, matplotlib.figure
        lr  = pd.read_csv(lr_pairs_csv)
        fig = rd.ligand_ranking(
            neighborhood_scores, celltype_adata, sc_adata,
            lr, lt_regulatory_potential,
            domain="0", receptor_exp=0.01,
            markers={"top": 10}, domain_celltypes=["top", 2])
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close("all")

    def test_custom_markers_dict(self, neighborhood_scores, celltype_adata,
                                  sc_adata, lt_regulatory_potential, lr_pairs_csv):
        import matplotlib.pyplot as plt, matplotlib.figure
        lr  = pd.read_csv(lr_pairs_csv)
        # Use genes that ARE in sc_adata.var_names and ARE marker genes
        custom_markers = {
            "TypeA": ["GENE0", "GENE1", "GENE2"],
            "TypeB": ["GENE8", "GENE9", "GENE10"],
        }
        fig = rd.ligand_ranking(
            neighborhood_scores, celltype_adata, sc_adata,
            lr, lt_regulatory_potential,
            domain="0", receptor_exp=0.01,
            markers=custom_markers,
            domain_celltypes=["TypeA", "TypeB"])
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close("all")

    def test_explicit_celltypes(self, neighborhood_scores, celltype_adata,
                                 sc_adata, lt_regulatory_potential, lr_pairs_csv):
        import matplotlib.pyplot as plt, matplotlib.figure
        lr  = pd.read_csv(lr_pairs_csv)
        fig = rd.ligand_ranking(
            neighborhood_scores, celltype_adata, sc_adata,
            lr, lt_regulatory_potential,
            domain="0", receptor_exp=0.01,
            markers={"top": 10},
            domain_celltypes=["TypeA", "TypeB"])
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close("all")


# ===========================================================================
# get_msig — human / mouse bundled file branches
# ===========================================================================

class TestGetMsigSpecies:

    def test_human_skips_if_file_missing(self):
        """The human bundled file may not exist in CI — skip gracefully."""
        try:
            result = rd.get_msig("human")
            assert isinstance(result, pd.DataFrame)
        except (FileNotFoundError, OSError):
            pytest.skip("Bundled human MSigDB file not present in this environment")

    def test_mouse_skips_if_file_missing(self):
        try:
            result = rd.get_msig("mouse")
            assert isinstance(result, pd.DataFrame)
        except (FileNotFoundError, OSError):
            pytest.skip("Bundled mouse MSigDB file not present in this environment")


# ===========================================================================
# create_cluster — additional branches
# ===========================================================================

class TestCreateClusterBranches:

    def test_use_pathway_subset(self, neighborhood_scores, msig_df, tmp_path):
        """use_pathway=True reads pathway_path CSV and filters the pathway dict."""
        import matplotlib.pyplot as plt
        # Write a pathway subset CSV (just one pathway name)
        pathway_csv = str(tmp_path / "pathway_subset.csv")
        pd.DataFrame({"pathways": ["HALLMARK_TEST_PATHWAY"]}).to_csv(pathway_csv, index=False)

        result = rd.create_cluster(
            neighborhood_scores,
            msig_df,
            method=None,
            use_pathway=True,
            pathway_path=pathway_csv,
            pathway_thresh=1,
        )
        # Only HALLMARK should survive the subset filter
        assert isinstance(result, dict)
        for key in result:
            assert key == "HALLMARK_TEST_PATHWAY" or key.startswith("cluster_")
        plt.close("all")

    def test_invalid_method_raises(self, neighborhood_scores, msig_df):
        """Passing an unrecognised method string should raise."""
        # The source raises a bare string (not Exception), which Python 3
        # converts to a TypeError at the raise site.
        with pytest.raises((TypeError, Exception)):
            rd.create_cluster(
                neighborhood_scores,
                msig_df,
                method="invalid_method",
                pathway_thresh=1,
            )


# ===========================================================================
# downstream_analysis — additional branches
# ===========================================================================

class TestDownstreamAnalysisBranches:

    def test_raises_when_ltpair_clusters_none(self, neighborhood_scores):
        """Passing ltpair_clusters=None raises (bare string, caught as TypeError)."""
        with pytest.raises((TypeError, Exception)):
            rd.downstream_analysis(neighborhood_scores, ltpair_clusters=None)

    def test_return_pcs_only(self, neighborhood_scores, clusters):
        """return_cluster=False, return_pcs=True should return just the pcs AnnData."""
        result = rd.downstream_analysis(
            neighborhood_scores,
            ltpair_clusters=clusters,
            resolution=0.5,
            return_cluster=False,
            return_pcs=True,
        )
        assert isinstance(result, ad.AnnData)
        # pcs var_names should be pathway names
        for varname in result.var_names:
            assert varname in clusters

    def test_pdf_output(self, neighborhood_scores, clusters, tmp_path):
        """pdf_path branch: a PDF file is produced containing diagnostic plots."""
        import os
        pdf_path = str(tmp_path / "downstream.pdf")
        rd.downstream_analysis(
            neighborhood_scores,
            ltpair_clusters=clusters,
            resolution=0.5,
            pdf_path=pdf_path,
            return_cluster=False,
            return_pcs=False,
        )
        assert os.path.exists(pdf_path)
        assert os.path.getsize(pdf_path) > 0


# ===========================================================================
# get_top_n_clust_pairs — internal helper
# ===========================================================================

class TestGetTopNClustPairs:
    """get_top_n_clust_pairs is an internal helper but is public enough to test."""

    def test_returns_list(self, neighborhood_scores, clusters):
        from Renoir.downstream import get_top_n_clust_pairs
        # downstream_analysis sets leiden; use the pre-set fixture value
        neighbscore_df = neighborhood_scores.to_df().T
        result = get_top_n_clust_pairs(neighborhood_scores, neighbscore_df, n=3)
        assert isinstance(result, list)

    def test_pairs_are_valid_var_names(self, neighborhood_scores):
        from Renoir.downstream import get_top_n_clust_pairs
        neighbscore_df = neighborhood_scores.to_df().T
        result = get_top_n_clust_pairs(neighborhood_scores, neighbscore_df, n=5)
        valid = set(neighborhood_scores.var_names)
        for pair in result:
            assert pair in valid


# ===========================================================================
# spot_v_spot — None cluster raises
# ===========================================================================

class TestSpotVSpotBranches:

    def test_raises_when_ltpair_clusters_none(
            self, neighborhood_scores, celltype_adata, tmp_path):
        with pytest.raises((TypeError, Exception)):
            rd.spot_v_spot(
                neighborhood_scores,
                celltype_adata,
                ltpair_clusters=None,
                pdf_path=str(tmp_path / "dummy.pdf"),
            )


# ===========================================================================
# sankeyPlot — additional branches
# ===========================================================================

class TestSankeyPlotBranches:

    def test_with_title(self, neighborhood_scores, celltype_adata):
        """Passing title= covers the title branch (L617)."""
        import matplotlib.pyplot as plt
        pairs = list(neighborhood_scores.var_names[:2])
        rd.sankeyPlot(
            neighborhood_scores,
            celltype_adata,
            ltpairs=pairs,
            n_celltype=2,
            clusters="All",
            title="Test Sankey",
        )
        plt.close("all")

    def test_invalid_clusters_raises(self, neighborhood_scores, celltype_adata):
        """Passing clusters as a non-list, non-'All' should raise."""
        with pytest.raises((TypeError, Exception)):
            rd.sankeyPlot(
                neighborhood_scores,
                celltype_adata,
                ltpairs=list(neighborhood_scores.var_names[:2]),
                clusters=999,
            )


# ===========================================================================
# pcs_v_neighbscore — None cluster raises
# ===========================================================================

class TestPcsVNeighbscoreBranches:

    def test_raises_when_ltpair_clusters_none(self, neighborhood_scores, tmp_path):
        with pytest.raises((TypeError, Exception)):
            rd.pcs_v_neighbscore(
                neighborhood_scores,
                ltpair_clusters=None,
                pdf_path=str(tmp_path / "dummy.pdf"),
            )


# ===========================================================================
# ligand_ranking — custom celltype_colors branch
# ===========================================================================

class TestLigandRankingBranches:

    def test_custom_celltype_colors(
            self, neighborhood_scores, celltype_adata, sc_adata,
            lt_regulatory_potential, lr_pairs_csv):
        """Pass an explicit celltype_colors dict to cover the else branch (L846)."""
        import matplotlib.pyplot as plt, matplotlib.figure
        from conftest import CELLTYPES
        lr  = pd.read_csv(lr_pairs_csv)
        custom_colors = {ct: f"#{i*50:02x}{i*30:02x}ff" for i, ct in enumerate(CELLTYPES)}
        fig = rd.ligand_ranking(
            neighborhood_scores,
            celltype_adata,
            sc_adata,
            lr,
            lt_regulatory_potential,
            domain="0",
            receptor_exp=0.01,
            markers={"top": 10},
            domain_celltypes=["top", 2],
            celltype_colors=custom_colors,
        )
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close("all")
