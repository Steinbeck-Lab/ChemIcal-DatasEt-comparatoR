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

The first step after the set-up is to load the SDFiles for comparison into the notebook with the 'import_data_as_dict' function. The function will create a dictionary which content will updated with every subsequent function used. This main dictionary contains a sub-dictionary for every SDFile named after the SDFile name.  
  
&nbsp;&nbsp;&nbsp;&nbsp; **function:**  
import_data_as_dict' 

### Overview

To get an overview of the datasets the number of molecules in every dataset can be determined using the 'get_number_of_molecules' function. The results for every dataset are inclued in the dataframe. Additionally one can create a grit image of molecules from each dataset with 'draw_molecules' function. The number of displayed molecules from every datasets can be specified (variable: number_of_mols) as well as the size of the images (variable: image_size). The images will also be saved in an output folder in a chosen format (variable: data_type).

### Identifier

For the subsequent comparison the molecules need a string representation and therefore one can use SMILES, InChI or InChIKey strings. The 'get_identifier_list_key' function gets the chosen identifier strings (variable: id_type) for every molecules and puts them in the dataframe.
If one uses datasets from a certain database and the database provides an own ID, this dataset ID can also be extracted from the SDFiles and stored in the dataframe. Therefore one can use the function 'get_database_id' and input the ID name (variable: id_name).

### Molecule comparision

With 'get_shared_molecules_key' one get the number and identifier string for those molecules which are present in all of the compared datasets. If there are no more than three datasets a Venn diagramm of the intersection of the datasets can be created using 'visualize_intersection'. The image will be saved in the output folder in a chosen format (variable: data_type).

### Descriptor and descriptor value distribution 

The function 'get_descriptor_list_key' utilize a callable function (variable: descriptor) form RDKit to get descriptor values for every molecule. For example the callable functions can be from the rdkit.Chem.rdMolDescriptors module the rdkit.Chem.Descriptors module. The values are saved in the dataframe under a chosen name (variable: descriptor_list_keyname).
To get and visualize the distribution of the descriptor values the function 'descriptor_counts_and_plot' can be used. The function distinguishes between continuous and discrete distributed descriptor vales and for the continuous values one can chose the binning size (variable: width_of_bins). The distribution of values is exported as csv-file and the visualization with a selectable format (variable: data_type) is also saved in the output folder.
With the database ID one can also search for the descriptor value of a specific molecule using the 'get_value_from_id' function. The function tells in which SDFile the molecule is found and the value of the descriptor.

### Lipinski Rule of 5

With 'get_lipinski_key' the number of broken Lipinki Rules for every molecule is calculated and a summary for every SDFile is created. Subsequently the 'lipinski_plot' function visualizes the number of broken rules. Again the results are exported as csv-file and the bar-plot is also saved with a selectable format (variable: data_type) in the output folder. 

### Chemical Space Visualization

For the visualization of the chemical space the chemplot module is used. The extended connectivity fingerprints can be specified with the fingerprint radius (variable: fp_radius) and the size (variable: fp_bits). For the dimension reduction PCA, t-SNE or UMAP can be chosen (variable: dimension_reduction). The chemical space plot is saved in the output folder except when choosing to create an interactive plot (variable: interactive). Then the plot will be displayed and can be manually saved.

