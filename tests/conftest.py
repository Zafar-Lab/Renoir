"""
Shared pytest fixtures for the Renoir test suite.

Design rules (learned from tracing the source):

1.  msig_df gene_symbol must include BOTH the ligand AND the target gene of
    every pair in neighborhood_scores.var_names.
    create_cluster builds pairs as: for genea in pathway_genes:
                                        for geneb in pathway_genes:
                                            if genea:geneb in var_names → keep
    So if GENE0:GENE5 is a var_name, BOTH 'GENE0' and 'GENE5' must appear in
    pathway gene_symbol.

2.  neighborhood_scores.X must be a scipy sparse matrix.
    spot_v_spot calls neighbscore_copy.X.toarray().

3.  sc_adata must have clearly differentially expressed genes per celltype.
    ligand_ranking runs rank_genes_groups on scrna and filters to
    pvals_adj < 0.05.  On random data no genes pass → markers_df is empty →
    downstream iloc crashes. Fix: give each celltype a set of exclusively
    highly-expressed genes.

4.  sc_adata must have N_CELLS_PER_CT cells per celltype (not just N_SPOTS
    total) so rank_genes_groups has enough cells in each group.

5.  neighborhood_scores.raw must be set before ligand_ranking is called
    (it calls neighbscore.raw.to_adata()).
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import pytest

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

N_SPOTS = 30
N_GENES = 40
N_CELLTYPES = 4
N_CELLS_PER_CT = 10  # cells per celltype in sc_adata (40 cells total)
N_LT_PAIRS = 12

CELLTYPES = ["TypeA", "TypeB", "TypeC", "TypeD"]
GENE_NAMES = [f"GENE{i}" for i in range(N_GENES)]

# Ligands: GENE0-GENE4 | Targets: GENE5-GENE9
# All are drawn from GENE_NAMES so get_ligand_target can find them
LIGANDS_LIST = [
    "GENE0",
    "GENE1",
    "GENE2",
    "GENE3",
    "GENE4",
    "GENE0",
    "GENE1",
    "GENE2",
    "GENE3",
    "GENE4",
    "GENE0",
    "GENE1",
]
TARGETS_LIST = [
    "GENE5",
    "GENE6",
    "GENE7",
    "GENE8",
    "GENE9",
    "GENE6",
    "GENE7",
    "GENE8",
    "GENE9",
    "GENE5",
    "GENE8",
    "GENE9",
]
PAIR_NAMES = [f"{lgt}:{tgt}" for lgt, tgt in zip(LIGANDS_LIST, TARGETS_LIST)]
LIGANDS_UNIQUE = list(dict.fromkeys(LIGANDS_LIST))  # GENE0..GENE4
TARGET_GENES = list(dict.fromkeys(TARGETS_LIST))  # GENE5..GENE9

# ALL genes that appear in any pair name (both sides)
ALL_PAIR_GENES = sorted(set(LIGANDS_LIST + TARGETS_LIST))

# Marker genes per celltype — designed so the full ligand→receptor→target chain
# in ligand_ranking always works regardless of which two celltypes are selected:
#   TypeA: GENE0-9  = all LIGANDS (GENE0-4) + all TARGETS (GENE5-9)
#   TypeB: GENE20-29 = all RECEPTORS (GENE20-24) used in lr_pairs_csv
#   TypeC: GENE10-19
#   TypeD: GENE30-39
# With top_celltypes=['TypeA','TypeB'], ligand_ranking finds:
#   - ligands (GENE0-4) as TypeA markers → ct_spec_mark['ligand'] populated
#   - receptors (GENE20-24) expressed in TypeB cells → receptors_out populated
#   - targets (GENE5-9) as TypeA markers → ct_spec_mark['target'] populated
#   → ligand_score_df gets target columns → iloc[:,1] succeeds
CT_MARKER_GENES = {
    "TypeA": [f"GENE{j}" for j in range(10)],  # GENE0-9: ligands + targets
    "TypeB": [f"GENE{j}" for j in range(20, 30)],  # GENE20-29: receptors
    "TypeC": [f"GENE{j}" for j in range(10, 20)],  # GENE10-19
    "TypeD": [f"GENE{j}" for j in range(30, 40)],  # GENE30-39
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _rng(seed=42):
    return np.random.default_rng(seed)


def _sparse(arr):
    return sp.csr_matrix(arr.astype(np.float32))


# ---------------------------------------------------------------------------
# Spatial transcriptomics AnnData (Visium-like)
# ---------------------------------------------------------------------------


@pytest.fixture
def st_adata():
    rng = _rng(0)
    spot_names = [f"spot_{i}" for i in range(N_SPOTS)]
    X = _sparse(rng.uniform(0.01, 5.0, size=(N_SPOTS, N_GENES)))
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=spot_names),
        var=pd.DataFrame(index=GENE_NAMES),
    )
    adata.obsm["spatial"] = rng.uniform(0, 500, size=(N_SPOTS, 2))
    adata.obs["celltype"] = pd.Categorical(np.resize(CELLTYPES, N_SPOTS))
    adata.obs["library_id"] = "library"
    # array_row / array_col are required by compute_neighborhood_scores (line 612
    # of renoir.py) which calls neighborhood(ST.obs['array_row'], ST.obs['array_col'])
    adata.obs["array_row"] = np.tile(np.arange(5), N_SPOTS // 5 + 1)[:N_SPOTS]
    adata.obs["array_col"] = np.tile(np.arange(6), N_SPOTS // 6 + 1)[:N_SPOTS]
    adata.uns["spatial"] = {
        "library": {
            "images": {"hires": rng.uniform(0, 1, (100, 100, 3)).astype(np.float32)},
            "scalefactors": {"tissue_hires_scalef": 1.0, "spot_diameter_fullres": 10},
        }
    }
    return adata


# ---------------------------------------------------------------------------
# scRNA-seq reference — each celltype has clearly DE marker genes
# so rank_genes_groups returns significant hits for ligand_ranking
# ---------------------------------------------------------------------------


@pytest.fixture
def sc_adata():
    n_cells = N_CELLS_PER_CT * N_CELLTYPES  # 40 cells
    rng = _rng(1)

    # Build expression: high for celltype-specific markers, low elsewhere
    X = np.zeros((n_cells, N_GENES), dtype=np.float32)
    celltype_labels = []
    for ct_idx, ct in enumerate(CELLTYPES):
        rows = range(ct_idx * N_CELLS_PER_CT, (ct_idx + 1) * N_CELLS_PER_CT)
        marker_cols = [int(g[4:]) for g in CT_MARKER_GENES[ct]]  # gene index
        for row in rows:
            X[row, marker_cols] = rng.uniform(8.0, 12.0, len(marker_cols))
            other_cols = [c for c in range(N_GENES) if c not in marker_cols]
            X[row, other_cols] = rng.uniform(0.0, 0.3, len(other_cols))
        celltype_labels += [ct] * N_CELLS_PER_CT

    cell_names = [f"cell_{i}" for i in range(n_cells)]
    obs = pd.DataFrame(index=cell_names)
    obs["celltype"] = pd.Categorical(celltype_labels)

    return ad.AnnData(
        X=_sparse(X),
        obs=obs,
        var=pd.DataFrame(index=GENE_NAMES),
    )


# ---------------------------------------------------------------------------
# Cell-type proportions
# ---------------------------------------------------------------------------


@pytest.fixture
def celltype_proportions():
    rng = _rng(2)
    spot_names = [f"spot_{i}" for i in range(N_SPOTS)]
    # Use a skewed Dirichlet so TypeA and TypeB are reliably the top two
    # celltypes. This guarantees ligand_ranking's top_2 celltypes always
    # include TypeA (which holds all ligand + target marker genes) and TypeB
    # (which holds all receptor marker genes), so the full signalling chain
    # ligand→receptor→target is always intact.
    raw = rng.dirichlet(alpha=[10.0, 10.0, 1.0, 1.0], size=N_SPOTS)
    return pd.DataFrame(raw, index=spot_names, columns=CELLTYPES)


# ---------------------------------------------------------------------------
# Celltype AnnData (for ligand_ranking / sankeyPlot)
# obs index MUST match neighborhood_scores obs index
# ---------------------------------------------------------------------------


@pytest.fixture
def celltype_adata(st_adata, celltype_proportions):
    spot_names = [f"spot_{i}" for i in range(N_SPOTS)]
    X = _sparse(celltype_proportions.values)
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=spot_names),
        var=pd.DataFrame(index=CELLTYPES),
    )
    adata.obsm["spatial"] = st_adata.obsm["spatial"]
    adata.uns["spatial"] = st_adata.uns["spatial"]
    adata.obs["library_id"] = "library"
    return adata


# ---------------------------------------------------------------------------
# Neighbourhood scores — mimics compute_neighborhood_scores output
# X is SPARSE (spot_v_spot calls .toarray())
# .raw is set (ligand_ranking calls neighbscore.raw.to_adata())
# ---------------------------------------------------------------------------


@pytest.fixture
def neighborhood_scores(st_adata):
    rng = _rng(4)
    spot_names = [f"spot_{i}" for i in range(N_SPOTS)]

    # Give the score matrix clear block structure so Leiden always finds
    # at least 2 clusters (random uniform scores collapse to 1 cluster).
    # Spots 0-14 have high scores for pairs 0-5; spots 15-29 for pairs 6-11.
    X = np.zeros((N_SPOTS, N_LT_PAIRS), dtype=np.float32)
    X[:15, :6] = rng.uniform(3.0, 5.0, (15, 6))
    X[:15, 6:] = rng.uniform(0.0, 0.3, (15, 6))
    X[15:, :6] = rng.uniform(0.0, 0.3, (15, 6))
    X[15:, 6:] = rng.uniform(3.0, 5.0, (15, 6))
    X = _sparse(X)

    obs = st_adata.obs.copy()
    obs.index = spot_names
    obs["leiden"] = pd.Categorical(
        ["0"] * 15 + ["1"] * 15  # pre-set 2 domains matching the block structure
    )

    adata = ad.AnnData(
        X=X,
        obs=obs,
        var=pd.DataFrame(index=PAIR_NAMES),
    )
    adata.obsm["spatial"] = st_adata.obsm["spatial"]
    adata.uns["spatial"] = st_adata.uns["spatial"]
    adata.obs["library_id"] = "library"
    adata.raw = adata.copy()
    return adata


# ---------------------------------------------------------------------------
# MSigDB pathway DataFrame
# CRITICAL: gene_symbol must include BOTH ligand AND target genes from PAIR_NAMES
# create_cluster matches pairs as genea:geneb where BOTH genea and geneb
# appear in the pathway's gene_symbol list.
# ---------------------------------------------------------------------------


@pytest.fixture
def msig_df():
    rows = []
    for gene in ALL_PAIR_GENES:  # includes GENE0..GENE4 (ligands) AND GENE5..GENE9 (targets)
        for gs in ["HALLMARK_TEST_PATHWAY", "KEGG_TEST_PATHWAY", "WP_TEST_PATHWAY"]:
            rows.append({"gs_name": gs, "gene_symbol": gene})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# MSigDB CSV on disk
# ---------------------------------------------------------------------------


@pytest.fixture
def msig_csv(tmp_path, msig_df):
    path = str(tmp_path / "msig.csv")
    msig_df.to_csv(path, index=False)
    return path


# ---------------------------------------------------------------------------
# Ligand-target regulatory potential (ligand_ranking)
# Indexed by LIGANDS_UNIQUE, columns = all GENE_NAMES
# ---------------------------------------------------------------------------


@pytest.fixture
def lt_regulatory_potential():
    rng = _rng(5)
    return pd.DataFrame(
        rng.uniform(0, 1, size=(len(LIGANDS_UNIQUE), N_GENES)),
        index=LIGANDS_UNIQUE,
        columns=GENE_NAMES,
    )


# ---------------------------------------------------------------------------
# Ligand-receptor CSV — receptors must be in sc_adata.var_names (GENE_NAMES)
# ---------------------------------------------------------------------------


@pytest.fixture
def lr_pairs_csv(tmp_path):
    df = pd.DataFrame(
        {
            "ligand": LIGANDS_UNIQUE,
            "receptor": [f"GENE{i + 20}" for i in range(len(LIGANDS_UNIQUE))],
        }
    )
    path = str(tmp_path / "lr_pairs.csv")
    df.to_csv(path, index=False)
    return path


# ---------------------------------------------------------------------------
# Ligand-target pairs CSV (for compute_neighborhood_scores integration tests)
# ---------------------------------------------------------------------------


@pytest.fixture
def lt_pairs_csv(tmp_path):
    df = pd.DataFrame({"ligand": LIGANDS_LIST, "target": TARGETS_LIST})
    path = str(tmp_path / "lt_pairs.csv")
    df.to_csv(path, index=False)
    return path


# ---------------------------------------------------------------------------
# Cell-type proportions CSV (for integration tests)
# ---------------------------------------------------------------------------


@pytest.fixture
def celltype_proportions_csv(tmp_path, celltype_proportions):
    path = str(tmp_path / "celltype_proportions.csv")
    celltype_proportions.to_csv(path)
    return path


# ---------------------------------------------------------------------------
# h5ad files on disk (for integration tests)
# ---------------------------------------------------------------------------


@pytest.fixture
def st_h5ad(tmp_path, st_adata):
    path = str(tmp_path / "ST.h5ad")
    st_adata.write_h5ad(path)
    return path


@pytest.fixture
def sc_h5ad(tmp_path, sc_adata):
    path = str(tmp_path / "SC.h5ad")
    sc_adata.write_h5ad(path)
    return path
