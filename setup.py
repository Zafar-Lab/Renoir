from setuptools import setup, find_packages
import codecs
import os

VERSION = '0.9.0'
DESCRIPTION = 'Charting spatial ligand-target activity using Renoir'

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
    ]
)
