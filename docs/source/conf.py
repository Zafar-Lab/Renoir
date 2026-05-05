# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys


sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------

project = 'Renoir'
copyright = '2023, Narein Rao'
author = 'Narein Rao'

# The full version, including alpha/beta/rc tags
release = '0.5.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc', 'nbsphinx', 'myst_parser',
]

autodoc_mock_imports = [
    "scanpy",
    "squidpy",
    "anndata",
    "numpy",
    "pandas",
    "scipy",
    "sklearn",
    "seaborn",
    "matplotlib",
    "hdbscan",
    "plotly",
    "harmonypy",
    "dynamictreecut",
    "distinctipy",
    "spatialdata",
    "spatialdata_io",
    "spatialdata_plot",
]

# Do not execute notebooks, render saved outputs only
nbsphinx_execute = 'never'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

napoleon_google_docstring = True
napoleon_numpy_docstring = True


# -- Options for HTML output -------------------------------------------------

html_theme = "pydata_sphinx_theme"

html_theme_options = {
    "github_url": "https://github.com/Zafar-Lab/Renoir",
    "navbar_end": ["navbar-icon-links"],
    "logo": {
        "text": "Renoir",
    },
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/Zafar-Lab/Renoir",
            "icon": "fa-brands fa-github",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/Renoir/",
            "icon": "fa-brands fa-python",
        },
    ],
    "show_version_warning_banner": True,
    "footer_start": ["copyright"],
    "footer_end": ["theme-version"],
}

# Optional: add a logo image if you have one
# html_logo = "_static/renoir_logo.png"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
