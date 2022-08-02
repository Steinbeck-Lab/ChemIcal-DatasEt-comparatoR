# ChemIcal DatasEt comparatoR (CIDER)
[![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/MIT)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-blue.svg)](https://GitHub.com/hannbus/ChemIcal_DatasEt_compaRator/graphs/commit-activity)
[![GitHub issues](https://img.shields.io/github/issues/hannbus/ChemIcal_DatasEt_compaRator.svg)](https://GitHub.com/hannbus/ChemIcal_DatasEt_compaRator/issues/)
[![GitHub contributors](https://img.shields.io/github/contributors/hannbus/ChemIcal_DatasEt_compaRator.svg)](https://GitHub.com/hannbus/ChemIcal_DatasEt_compaRator/graphs/contributors/)
[![DOI](https://zenodo.org/badge/501949039.svg)](https://zenodo.org/badge/latestdoi/501949039)
[![Documentation Status](https://readthedocs.org/projects/ChemIcal_DatasEt_compaRator/badge/?version=latest)](https://ChemIcal_DatasEt_compaRator.readthedocs.io/en/latest/?badge=latest)
[![GitHub release](https://img.shields.io/github/release/hannbus/ChemIcal_DatasEt_compaRator.svg)](https://github.com/hannbus/ChemIcal_DatasEt_compaRator/releases/)
[![PyPI version fury.io](https://badge.fury.io/py/cider-chem.svg)](https://pypi.python.org/pypi/cider-chem/)

[![GitHub Logo](https://github.com/hannbus/ChemIcal_DatasEt_compaRator/blob/main/Workflow/Cider.png?raw=true)](https://pypi.python.org/pypi/cider-chem/)

- ChemIcal DatasEt comparatoR (CIDER) is a Python package and ready-to-use Jupyter Notebook workflow which primarily utilizes RDKit to compare two or more chemical structure datasets (SD files) in different aspects, e.g. size, overlap, molecular descriptor distributions, chemical space clustering, etc., most of which can be visually inspected in the notebook.

## Usage
-  To use CIDER, clone the repository to your local disk and make sure you install all the necessary requirements.

### We recommend to use CIDER inside a Conda environment to facilitate the installation of the dependencies.

- Conda can be downloaded as part of the [Anaconda](https://www.anaconda.com/) or the [Miniconda](https://conda.io/en/latest/miniconda.html) platforms (Python 3.8). We recommend to install miniconda3. Using Linux you can get it with:

```shell
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
```
### Installation

```shell
$ git clone https://github.com/hannbus/ChemIcal_DatasEt_compaRator.git
$ cd ChemIcal_DatasEt_compaRator
$ conda create --name ChemIcal_DatasEt_compaRator python=3.10
$ conda activate ChemIcal_DatasEt_compaRator
$ conda install pip
$ python -m pip install -U pip #Upgrade pip
$ pip install .
```
### Alternative
```shell
$ python -m pip install -U pip #Upgrade pip
$ pip install git+https://github.com/hannbus/ChemIcal_DatasEt_compaRator.git
```

### Install from PyPI
```shell
$ pip install cider-chem
```

### Basic usage: 
- Todo

### Cite us
- Busch, H., Schaub, J., Brinkhaus, H. O., Rajan, K., & Steinbeck, C. (2022). ChemIcal DatasEt comparatoR CIDER (Version 0.0.1-dev) [Computer software]. https://doi.org/10.5281/zenodo.6630494

## More information about our research group

[![GitHub Logo](https://github.com/Kohulan/DECIMER-Image-to-SMILES/blob/master/assets/CheminfGit.png?raw=true)](https://cheminf.uni-jena.de)
