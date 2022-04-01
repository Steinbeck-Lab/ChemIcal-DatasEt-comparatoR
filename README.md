# Comparison of datasets

This repository contains a Jupyter Notebook that can be used to compare SDFiles with each other and create visuale output from the comparison. 

## Set-up

The notebook can be run in a conda environment with RDKit.

```shell
$ conda create -c rdkit -n myenv rdkit
$ activate new_env
$ conda install -c conda-forge rdkit
$ conda install jupyter
$ pip install matplotlib==3.5.1 seaborn==0.11.2 chemplot==1.2.0 matplotlib_venn==0.11.6
```
## How to use the notebook

### Import datasets

The first step after the set-up is to load the SDFiles for comparison into the notebook with the 'import_data_as_dict function'. The function will create a dictionary which content will updated with every subsequent function used. This main dictionary contains a sub-dictionary for every SDFile named after the SDFile name.

### Overview

To get an overview of the datasets the number of molecules in every dataset can be determined using the 'get_number_of_molecules' function. The results for every dataset are inclued in the dataframe. Additionally one can create a grit image of molecules from each dataset with 'draw_molecules' function. The number of displayed molecules from every datasets can be specified as well as the size of the images. The images will also be saved in an output folder.

### Identifier

For the subsequent comparison the molecules need a string representation and therefore one can use SMILES, InChI or InChIKey strings. The 'get_identifier_list_key' function gets the chosen identifier strings for every molecules and puts them in the dataframe.
If one uses datasets from a certain database the database's ID for the molecules can also be extracted from the SDFiles and stored in the dataframe ('get_database_id'). 

### Molecule comparision

With 'get_shared_molecules_key' one get the number and identifier string for those molecules which are present in all of the compared datasets. If there are no more than three datasets a Venn diagramm of the intersection of the datasets can be created using 'visualize_intersection'. The image will be saved in the output file.

### Descriptor and descriptor value distribution 



### Lipinski Rule of 5

### Chemical Space Visualization

