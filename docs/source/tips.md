# Tips for Running Renoir Effectively

This page collects practical guidance for getting the most out of Renoir across
different spatial transcriptomics technologies and experimental designs. Each tip
references the exact parameters you need to adjust.

> **NOTE:** We are currently adding more tips to this page to improve tool performance.

---

## 1. Choosing the Optimal Neighborhood Size

The neighborhood definition is set in `compute_neighborhood_scores` via the
`technology`, `use_radius`, and `radius` parameters. The right choice depends on
your platform's resolution and the biological scale of the signalling you expect.

### Visium (55 µm spots, ~6 spot diameter)

Visium spots are large and sparsely packed. The default hexagonal ring
neighborhood (the 6 immediately adjacent spots) is a good starting point and
requires no extra parameters:

```python
neighborhood_scores = Renoir.compute_neighborhood_scores(
    ...,
    single_cell=False,   # spot-level mode
    use_radius=False,    # use default hex-ring neighbor definition
)
```

If the signal looks too local (e.g., domains are only 1–2 spots wide), expand
to a wider neighborhood by switching to radius mode:

```python
neighborhood_scores = Renoir.compute_neighborhood_scores(
    ...,
    single_cell=False,
    use_radius=True,
    radius=200,         # 200 coordinate units — check your dataset's coordinate scale
)
```

### Visium HD (8 µm bins)

Visium HD produces much denser data. Because individual bins are small,
a larger radius is needed to capture biologically meaningful neighborhoods.
Start with a moderate value and increase if domains appear fragmented:

```python
neighborhood_scores = Renoir.compute_neighborhood_scores(
    ...,
    single_cell=False,
    use_radius=True,
    radius=50,          # in coordinate units — adjust based on your dataset's scale
)
```

### CosMx / Xenium / MERSCOPE (single-cell resolution)

For single-cell platforms, enable `single_cell=True` and set `radius` in
the same coordinate units as your AnnData's spatial coordinates:

```python
neighborhood_scores = Renoir.compute_neighborhood_scores(
    ...,
    single_cell=True,
    use_radius=True,
    radius=150,         # in coordinate units — check obsm['spatial'] for your dataset's scale
)
```

> **Important:** `radius` is always in the same units as the X, Y coordinates
> stored in your AnnData (typically `obsm['spatial']`). These units differ
> between datasets and platforms. Always inspect your coordinates before
> choosing a radius value:
>
> ```python
> import pandas as pd
> coords = pd.DataFrame(neighborhood_scores.obsm['spatial'], columns=['x', 'y'])
> print(coords.describe())  # check the range to understand the coordinate scale
> ```

**General rule of thumb:**

| Technology | Neighborhood mode | `single_cell` |
|---|---|---|
| Visium | `use_radius=False` (hex ring default) | `False` |
| Visium HD | `use_radius=True`, tune radius to your coordinate scale | `False` |
| CosMx | `use_radius=True`, tune radius to your coordinate scale | `True` |
| Xenium | `use_radius=True`, tune radius to your coordinate scale | `True` |
| MERSCOPE | `use_radius=True`, tune radius to your coordinate scale | `True` |

> **Tip:** If you are unsure, run `downstream_analysis` at two or three radii
> and compare the spatial coherence of the resulting domains. Domains that look
> biologically sensible (contiguous regions matching known tissue compartments)
> are a good sign that the neighborhood size is appropriate.

---

## 2. Creating Your Own Curated Ligand–Receptor and Ligand–Target Pair Lists

Renoir accepts custom pair lists via the `ligand_receptor_path` and `pairs_path`
arguments in `compute_neighborhood_scores`. Both are simply CSV files with
specific column names.

### Ligand–receptor pairs (`ligand_receptor_path`)

The file must have at least two columns named `ligand` and `receptor`:

```
ligand,receptor
TGFB1,TGFBR1
TGFB1,TGFBR2
IL6,IL6R
VEGFA,FLT1
VEGFA,KDR
```

**Sources to build from:**
- [CellChat DB](https://github.com/sqjin/CellChat) — curated, literature-backed
- [OmniPath](https://omnipathdb.org/) — large, multi-resource aggregate
- [NATMI](https://github.com/asrhou/NATMI) — what the tutorials use
- [CellPhoneDB](https://www.cellphonedb.org/) — includes multi-subunit complexes

To filter to your biology of interest (e.g., only cytokine signalling):

```python
import pandas as pd

lr = pd.read_csv('All_human_lrpairs.csv')

# Keep only cytokine-related ligands (example gene list)
cytokines = ['IL6', 'IL1B', 'TNF', 'IFNG', 'CXCL10', 'CCL2']
lr_cytokine = lr[lr['ligand'].isin(cytokines)]
lr_cytokine.to_csv('cytokine_lrpairs.csv', index=False)
```

### Ligand–target pairs (`pairs_path`)

The file must have columns named `ligand` and `target`:

```
ligand,target
IL6,STAT3
IL6,MYC
TGFB1,SNAI1
VEGFA,HIF1A
```

**How to generate a ranked list:**

The bundled top-N files (`top_10_target_opt_both_ordered.csv`,
`top_100_target_opt_both_ordered.csv`) are derived from NicheNet's regulatory
potential matrix. To build your own ranked list from NicheNet:

```python
import pandas as pd

# Load the full NicheNet regulatory potential matrix
# (available at https://zenodo.org/record/3260758)
reg = pd.read_csv('ligand_target_matrix.csv', index_col=0)

# Keep only your ligands of interest
ligands_of_interest = ['IL6', 'TGFB1', 'VEGFA']
reg_subset = reg.loc[ligands_of_interest]

# For each ligand, keep the top-N targets by regulatory potential score
top_n = 20
rows = []
for ligand in reg_subset.index:
    top_targets = reg_subset.loc[ligand].nlargest(top_n)
    for target, score in top_targets.items():
        rows.append({'ligand': ligand, 'target': target, 'score': score})

custom_pairs = pd.DataFrame(rows).sort_values('score', ascending=False)
custom_pairs.to_csv('my_ligand_target_pairs.csv', index=False)
```

> **Tip:** Start with a smaller, focused pair list (top 10–20 pairs per ligand)
> rather than top 100+. Downstream clustering is faster and domains are often
> more interpretable when the signal is not diluted by low-scoring pairs.

---

## 3. Generating Your Own Ligand–Target Regulatory Potential Scores

The `ligand_target_regulatory_potential` object used by `ligand_ranking` is a
precomputed dictionary (or DataFrame) mapping each ligand to regulatory potential
scores across its top target genes. There are three ways to produce this.

### Option A — Use the NicheNet matrix directly (recommended)

The NicheNet team publish a precomputed human and mouse matrix on Zenodo
(record 3260758). Load it and convert to the format Renoir expects:

```python
import pandas as pd, pickle

# Download from: https://zenodo.org/record/3260758
reg = pd.read_csv('ligand_target_matrix.csv', index_col=0)  # ligands × targets

# Keep only the top 500 targets per ligand (matching the bundled file)
top_500 = reg.apply(lambda row: row.nlargest(500), axis=1)

# Save as a pickle
with open('top_500_target_opt_both_scores.pkl', 'wb') as f:
    pickle.dump(top_500, f)
```

### Option B — Use an existing regulatory database

Rather than inferring a gene regulatory network from data (which is unreliable
in practice), it is better to derive regulatory potentials from a curated
transcription factor–target database. Two well-maintained options are:

- **[DoRothEA](https://saezlab.github.io/dorothea/)** — curated TF–target
  regulons with confidence levels (A–E); available for human and mouse.
- **[CollecTRI](https://github.com/saezlab/CollecTRI)** — a comprehensive
  signed TF–target network compiled from literature and ChIP-seq data.

Both can be accessed via the `decoupler` Python package and converted to the
format Renoir expects:

```python
import decoupler as dc
import pandas as pd, pickle

# Load DoRothEA regulons (human, confidence levels A and B only)
dorothea = dc.get_dorothea(organism='human', levels=['A', 'B'])
# dorothea has columns: 'source' (TF), 'target', 'weight', 'confidence'

# Map TFs to their upstream ligands using your ligand-receptor table
# (i.e., if a receptor activates a TF, then the ligand -> TF -> target chain
# gives you the ligand's regulatory potential over that target)
lr = pd.read_csv('All_human_lrpairs.csv')  # columns: ligand, receptor

# Build a ligand -> target score table via receptor -> TF -> target
rows = []
for _, row in lr.iterrows():
    ligand, receptor = row['ligand'], row['receptor']
    # Find targets regulated by TFs known to be downstream of this receptor
    # (requires a receptor -> TF mapping from e.g. OmniPath kinase-substrate)
    tf_targets = dorothea[dorothea['source'] == receptor]
    for _, t in tf_targets.iterrows():
        rows.append({'ligand': ligand, 'target': t['target'], 'score': t['weight']})

reg_potential = pd.DataFrame(rows)
reg_potential = reg_potential.groupby(['ligand', 'target'])['score'].max().unstack(fill_value=0)

with open('dorothea_reg_potential.pkl', 'wb') as f:
    pickle.dump(reg_potential, f)
```

> **Note:** A full receptor → TF mapping requires an additional signalling
> resource such as [OmniPath](https://omnipathdb.org/) (via `omnipath` Python
> package) or [SignaLink](https://signalink.org/). The NicheNet matrix (Option A)
> already integrates all of these layers internally, which is why it remains the
> recommended starting point.

### Option C — Restrict the bundled matrix to your ligands of interest

If you only care about a subset of ligands, slice the bundled pickle to reduce
memory and computation during `ligand_ranking`:

```python
import pickle, pandas as pd

with open('top_500_target_opt_both_scores.pkl', 'rb') as f:
    reg = pickle.load(f)

ligands_of_interest = ['IL6', 'TGFB1', 'VEGFA', 'CXCL10']
reg_subset = reg.loc[[l for l in ligands_of_interest if l in reg.index]]

with open('subset_reg_potential.pkl', 'wb') as f:
    pickle.dump(reg_subset, f)
```

---

## 4. Defining Communication Domains

Renoir offers two strategies for identifying communication domains. Which one
to use depends on whether you want data-driven discovery or hypothesis-driven
annotation.

### Strategy A — Data-driven Leiden clustering (`downstream_analysis`)

Use this when you have no prior knowledge of how many regions exist or where
they are. Renoir reduces the score matrix to pathway PCs and clusters with
Leiden:

```python
neighbscore_copy, pcs = Renoir.downstream_analysis(
    neighborhood_scores,
    ltpair_clusters=pathways,
    resolution=0.6,        # key parameter — see below
    return_cluster=True,
    return_pcs=True,
)
```

**Tuning `resolution`:**

The effect of resolution is highly data-dependent and no single range applies
universally. From experiments across datasets, the following has been observed
as a rough guide **for Visium data**:

| Resolution | Typical outcome (Visium) |
|---|---|
| 0.1 – 0.3 | 3–5 broad domains (tumor / stroma / immune) |
| 0.4 – 0.7 | 5–10 mid-grain domains — good starting point |
| 0.8 – 1.5 | 10+ fine-grain domains; risk of over-fragmentation |

For **high-resolution single-cell platforms** (CosMx, Xenium, MERSCOPE), the
optimal resolution is often much lower. Values below 0.1 have been used
successfully in practice — the much higher cell density means that even a very
low resolution produces a meaningful number of domains. Start at 0.05 and
increase gradually.

The best approach regardless of platform is to sweep a range and evaluate
visually:

```python
for res in [0.3, 0.5, 0.8, 1.0]:
    copy, _ = Renoir.downstream_analysis(
        neighborhood_scores,
        ltpair_clusters=pathways,
        resolution=res,
        return_cluster=True,
    )
    sc.pl.spatial(copy, color='leiden', title=f'resolution={res}', size=1.4)
```

### Strategy B — Cell-type-informed clustering (`spot_v_spot`)

Use `spot_v_spot` when cell-type co-localisation matters as much as signalling
similarity. It combines the L–T score matrix with pairwise cosine similarity of
cell-type abundances, so domains reflect both what signalling is happening and
which cell types are co-localised:

```python
rd.spot_v_spot(
    neighborhood_scores,
    celltype,
    resolution=0.8,
    ltpair_clusters=pathways,
    pdf_path='spot_v_spot_output.pdf',
)
```

### Strategy C — User-defined regions of interest

If you already know which spots belong to a region of interest (e.g., from
pathologist annotations, RCTD output, or manual selection in Loupe Browser),
assign the labels directly to `obs['leiden']` and skip clustering entirely:

```python
import pandas as pd

# Load your manual annotations: a CSV with columns 'barcode' and 'region'
annotations = pd.read_csv('manual_annotations.csv', index_col='barcode')

# Assign to the neighborhood scores object
neighborhood_scores.obs['leiden'] = annotations.loc[
    neighborhood_scores.obs_names, 'region'
].astype('category')

# Now run DE and ligand ranking as normal — Renoir treats the labels the same
# regardless of whether they came from clustering or manual annotation
sc.tl.rank_genes_groups(neighborhood_scores, 'leiden', method='wilcoxon')
```

> **Tip:** Manual annotations and Leiden clusters can be mixed. Annotate the
> regions you understand (tumor core, necrotic zone) and let Leiden cluster the
> rest — then merge the two `obs` columns before running `ligand_ranking`.

---

## 5. Choosing Cell Types and Providing Custom Markers for Ligand Ranking

`ligand_ranking` has two parameters that give you fine-grained control over
which cell types are analysed and which ligand–target pairs are used as the
ranking signal.

### Controlling which cell types are included (`domain_celltypes`)

**Option 1 — Top-N by abundance (default):**

```python
fig = Renoir.ligand_ranking(
    ...,
    domain_celltypes=['top', 5],   # top 5 most abundant cell types in the domain
)
```

**Option 2 — Explicit cell type list:**

Pass a list of cell-type names to restrict analysis to only those types,
regardless of their abundance in the domain. This is useful when you have a
biological hypothesis (e.g., "I only care about T cell – tumour interactions"):

```python
fig = Renoir.ligand_ranking(
    ...,
    domain_celltypes=['Cancer Basal SC', 'T cells CD8+', 'Macrophage'],
)
```

The cell-type names must match the column names in your `celltype` AnnData
exactly (check with `celltype.var_names` or `celltype.to_df().columns`).

### Controlling the ranking signal (`markers`)

The `markers` parameter defines which ligand–target pairs are used to score
each ligand's predicted activity in the domain.

**Option 1 — Top-N DE marker pairs (default):**

```python
fig = Renoir.ligand_ranking(
    ...,
    markers={'top': 100},   # use the top 100 DE pairs of the domain
)
```

Increase `top` if you want a broader ranking signal; decrease it to focus on
only the most domain-specific pairs.

**Option 2 — User-defined cell-type marker genes:**

Instead of using DE pairs from the domain, you can supply your own curated
marker genes per cell type as a dictionary. The keys are cell-type names and
the values are lists of marker genes. Renoir uses these to identify which
ligands best explain the activity of those markers within the domain:

```python
# Define known marker genes per cell type
custom_markers = {
    'Cancer Basal SC': ['KRT5', 'KRT14', 'TP63', 'CDH3'],
    'T cells CD8+':    ['CD8A', 'CD8B', 'GZMB', 'PRF1', 'IFNG'],
    'Macrophage':      ['CD68', 'CD163', 'MRC1', 'CSF1R'],
    'CAFs myCAF-like': ['ACTA2', 'FAP', 'POSTN', 'THY1'],
}

fig = Renoir.ligand_ranking(
    ...,
    markers=custom_markers,
)
```

This is particularly useful when:
- The domain of interest is small and `rank_genes_groups` returns few DE pairs.
- You have strong prior knowledge about which cell types are biologically
  relevant in the domain and want to anchor the ranking to known biology.
- You want results that are directly comparable across datasets or studies,
  where DE-derived markers would differ due to batch or composition differences.

The marker gene names must match the gene names in your scRNA-seq reference
AnnData (`SC.var_names`).

### Filtering ligands by receptor expression (`receptor_exp`)

`receptor_exp` sets the minimum fraction of spots in the domain that must
express a ligand's receptor. Raise it to be more stringent (only ligands the
domain is definitely listening to), or lower it to capture weak but potentially
important signals:

```python
# Strict: receptor must be expressed in ≥10% of domain spots (default)
fig = Renoir.ligand_ranking(..., receptor_exp=0.1)

# Permissive: receptor expressed in ≥1% of domain spots
fig = Renoir.ligand_ranking(..., receptor_exp=0.01)
```

> **Tip:** If `ligand_ranking` returns very few ligands, it is usually because
> `receptor_exp` is too high for your data. Try lowering it to `0.01` and check
> whether the additional ligands make biological sense before committing to the
> lower threshold.

---

## Quick Reference

| Goal | Parameter | Where |
|---|---|---|
| Wider neighborhood | `use_radius=True, radius=N` (N in your coordinate units) | `compute_neighborhood_scores` |
| Single-cell mode | `single_cell=True` | `compute_neighborhood_scores` |
| Fewer, broader domains (Visium) | `resolution=0.2` | `downstream_analysis` / `spot_v_spot` |
| More, finer domains (Visium) | `resolution=1.0` | `downstream_analysis` / `spot_v_spot` |
| High-res platforms (CosMx etc.) | `resolution=0.05` or lower | `downstream_analysis` / `spot_v_spot` |
| Cell-type-aware domains | use `spot_v_spot` instead of `downstream_analysis` | — |
| Manual region annotations | assign `obs['leiden']` directly | — |
| Focus on specific cell types | `domain_celltypes=['CellA', 'CellB']` | `ligand_ranking` |
| DE-driven ranking signal | `markers={'top': 100}` | `ligand_ranking` |
| Custom marker gene ranking | `markers={'CellType': ['gene1', 'gene2']}` | `ligand_ranking` |
| Strict receptor filter | `receptor_exp=0.1` | `ligand_ranking` |
| Permissive receptor filter | `receptor_exp=0.01` | `ligand_ranking` |
