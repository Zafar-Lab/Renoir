from setuptools import setup

VERSION = "1.0.0-beta"
DESCRIPTION = "Charting spatial ligand-target activity using Renoir"

# Setting up
setup(
    name="Renoir",
    version=VERSION,
    author="Narein Rao (Zafar-Lab)",
    description=DESCRIPTION,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
    ],
    extras_require={
        "test": [
            "pytest>=7.4",
            "pytest-cov",
            "anndata",
            "scanpy",
            "pandas",
            "numpy",
            "scipy",
            "matplotlib",
            "distinctipy",
            "hdbscan",
            "plotly",
            "dynamictreecut",
            "harmonypy",
            "spatialdata",
            "spatialdata-io",
            "spatialdata-plot",
            "igraph",
            "leidenalg",
        ]
    },
)
