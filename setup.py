from setuptools import setup, find_packages

VERSION = "1.0.0"
DESCRIPTION = (
    "Charting spatial ligand-target activity using Renoir " "(ligand-taRgEt iNteractions acrOss spatIal topogRaphy)"
)
LONG_DESCRIPTION = open("README.md", encoding="utf-8").read()

setup(
    name="renoir-spatial",
    version=VERSION,
    author="Narein Rao (Zafar-Lab)",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/Zafar-Lab/Renoir",
    project_urls={
        "Documentation": "https://renoir.readthedocs.io/en/latest/",
        "Bug Tracker": "https://github.com/Zafar-Lab/Renoir/issues",
        "Source Code": "https://github.com/Zafar-Lab/Renoir",
    },
    packages=find_packages(),
    python_requires=">=3.11",
    install_requires=[
        # ── Core scientific stack ─────────────────────────────────────────
        "numpy",
        "scipy",
        "pandas",
        "scikit-learn>=1.6.0",
        # ── Single-cell / spatial analysis ───────────────────────────────
        "scanpy>=1.10.4",
        "squidpy>=1.6.2",
        "anndata",
        # ── Clustering ───────────────────────────────────────────────────
        "leidenalg",
        "igraph",
        "hdbscan>=0.8.40",
        # ── Visualisation ─────────────────────────────────────────────────
        "plotly>=5.24.1",
        "distinctipy>=1.3.4",
        "matplotlib",
        "seaborn",
        # ── Pip-only helpers ─────────────────────────────────────────────
        "harmonypy>=0.0.10",
        "dynamictreecut>=0.1.1",
        # ── SpatialData ecosystem ─────────────────────────────────────────
        "spatialdata>=0.7.0",
        "spatialdata-io",
        "spatialdata-plot",
        # ── Utilities ─────────────────────────────────────────────────────
        "tqdm",
        "setuptools<82",
    ],
    extras_require={
        "test": [
            "pytest>=7.4",
            "pytest-cov",
        ],
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    keywords=[
        "spatial transcriptomics",
        "cell-cell communication",
        "ligand-target",
        "Visium",
        "CosMx",
        "scRNA-seq",
        "bioinformatics",
    ],
)
