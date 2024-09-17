<p align="center"><a href="https://chemical-dataset-comparator.readthedocs.io/en/latest/" target="_blank"><img src="Cider_white.png" width="400" alt="CMS Logo"></a></p>

[![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/MIT)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-blue.svg)](https://GitHub.com/Steinbeck-Lab/ChemIcal-DatasEt-comparatoR/graphs/commit-activity)
[![GitHub issues](https://img.shields.io/github/issues/Steinbeck-Lab/ChemIcal-DatasEt-comparatoR.svg)](https://GitHub.com/Steinbeck-Lab/ChemIcal-DatasEt-comparatoR/issues/)
[![GitHub contributors](https://img.shields.io/github/contributors/Steinbeck-Lab/ChemIcal-DatasEt-comparatoR.svg)](https://GitHub.com/Steinbeck-Lab/ChemIcal_DatasEt_compaRator/graphs/contributors/)
[![GitHub release](https://img.shields.io/github/release/Steinbeck-Lab/ChemIcal-DatasEt-comparatoR.svg)](https://github.com/Steinbeck-Lab/ChemIcal-DatasEt-comparatoR/releases/)
[![RDKit badge](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)
![Workflow](https://github.com/Steinbeck-Lab/ChemIcal-DatasEt-comparatoR/actions/workflows/main.yml/badge.svg)
[![DOI](https://zenodo.org/badge/501949039.svg)](https://zenodo.org/badge/latestdoi/501949039)
[![Documentation Status](https://readthedocs.org/projects/chemical-dataset-comparator/badge/?version=latest)](https://chemical-dataset-comparator.readthedocs.io/en/latest/?badge=latest)
[![PyPI version fury.io](https://badge.fury.io/py/cider-chem.svg)](https://pypi.python.org/pypi/cider-chem/)

# Overview of ChemIcal DatasEt comparatoR (CIDER) :globe_with_meridians:
- ChemIcal DatasEt comparatoR (CIDER) is a Python package and ready-to-use Jupyter Notebook workflow which primarily utilizes RDKit to compare two or more chemical structure datasets (SD files) in different aspects, e.g. size, overlap, molecular descriptor distributions, chemical space clustering, etc., most of which can be visually inspected in the notebook.

## Usage
-  To use CIDER, clone the repository to your local disk and make sure you install all the necessary requirements.

### We recommend to use CIDER inside a Conda environment to facilitate the installation of the dependencies.

- Conda can be downloaded as part of the [Anaconda](https://www.anaconda.com/) or the [Miniconda](https://conda.io/en/latest/miniconda.html) platforms (Python 3.10). We recommend to install miniconda3. Using Linux you can get it with:

```shell
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
```
### Installation

```shell
$ git clone https://github.com/Steinbeck-Lab/ChemIcal_DatasEt_compaRator.git
$ cd ChemIcal_DatasEt_compaRator
$ conda create --name cider_chem python=3.10
$ conda activate cider_chem
$ conda install pip
$ python -m pip install -U pip #Upgrade pip
$ pip install .
```
- Note: Make sure all installations are working correctly by running the tests. You can do this by running the pytest command in the repository root folder.

### Alternative
```shell
$ python -m pip install -U pip #Upgrade pip
$ pip install git+https://github.com/Steinbeck-Lab/ChemIcal_DatasEt_compaRator.git
```

### Install from PyPI
```shell
$ pip install cider-chem
```

### Basic usage:
```python
from CIDER import ChemicalDatasetComparator
cider = ChemicalDatasetComparator()

data_dir = './data/'  # dir with sd files containing molecules
testdict = cider.import_as_data_dict(data_dir)
cider.get_number_of_molecules(testdict)

```
### Documentation
- The documentation for the CIDER package can be found [here](https://chemical-dataset-comparator.readthedocs.io/en/latest/?badge=latest).

### Cite us
- Busch, H., Schaub, J., Brinkhaus, H. O., Rajan, K., & Steinbeck, C. (2022). ChemIcal DatasEt comparatoR CIDER (Version 0.0.1-dev) [Computer software]. https://doi.org/10.5281/zenodo.6630494

## Maintained by :wrench:
ChemIcal DatasEt comparatoR is developed and maintained by the [Steinbeck group](https://cheminf.uni-jena.de) at the [Friedrich Schiller University](https://www.uni-jena.de/en/) Jena, Germany.
The code for this web application is released under the [MIT license](https://opensource.org/licenses/MIT). Copyright © CC-BY-SA 2024



[![GitHub Logo](https://github.com/Kohulan/DECIMER-Image-to-SMILES/blob/master/assets/CheminfGit.png?raw=true)](https://cheminf.uni-jena.de)
