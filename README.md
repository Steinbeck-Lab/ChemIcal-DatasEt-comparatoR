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
  
| **function**  |  import_data_as_dict |   
|---|---|
|  **parameter**  | path_to_data : *str*  |
|   | &nbsp;&nbsp; Path to the SDFiles which are to be compared. |

### Overview

To get an overview of the datasets the number of molecules in every dataset can be determined using the 'get_number_of_molecules' function. The results for every dataset are inclued in the dataframe.  
  
| **function**  |  get_number_of_molecules |   
|---|---|
|  **parameter**  | all_dict : *dict*  |
|   | &nbsp;&nbsp; Name of the dictionary with the imported SDFiles. | 
  
Additionally one can create a grit image of molecules from each dataset with 'draw_molecules' function. The number of displayed molecules from every datasets can be specified (variable: number_of_mols) as well as the size of the images (variable: image_size). The images will also be saved in an output folder in a chosen format (variable: data_type).  
  
| **function**  |  draw_molecules |   
|---|---|
|  **parameter**  | all_dict : *dict*  |
|   | &nbsp;&nbsp; Name of the dictionary with the imported SDFiles. |
|   | number_of_mols : *int, default=12*  |
|   | &nbsp;&nbsp; Number of molecules from every dataset to be displayed. |
|   | image_size : *int, default=200*  |
|   | &nbsp;&nbsp; Image size for one molecule in the grid image. |
|   | data_type : *str, default='png'*  |
|   | &nbsp;&nbsp; Data type for the exported image. |

![mol_grit_set_phenole](https://user-images.githubusercontent.com/95417135/165088366-00867664-d01c-4748-97fe-1337cf2bffae.png)
|:--:| 
| *Example of the 'draw_molecules' function with the visualization of set_phenole.sdf as grit image (number_of_molecules = 6).* |


### Identifier

For the subsequent comparison the molecules need a string representation and therefore one can use SMILES, InChI or InChIKey strings. The 'get_identifier_list_key' function gets the chosen identifier strings (variable: id_type) for every molecules and puts them in the dataframe.  
  
| **function**  |  get_identifier_list_key |   
|---|---|
|  **parameter**  | all_dict : *dict*  |
|   | &nbsp;&nbsp; Name of the dictionary with the imported SDFiles. |
|   | id_type : *str, default='inchi'*  |
|   | &nbsp;&nbsp; Type of Identifier ("inchi", "inchikey" or "smiles"). |  
  
If one uses datasets from a certain database and the database provides an own ID, this dataset ID can also be extracted from the SDFiles and stored in the dataframe. Therefore one can use the function 'get_database_id' and input the ID name (variable: id_name).  
  
| **function**  |  get_database_id |   
|---|---|
|  **parameter**  | all_dict : *dict*  |
|   | &nbsp;&nbsp; Name of the dictionary with the imported SDFiles. |
|   | id_name : *str*  |
|   | &nbsp;&nbsp; Name of the database ID in the SDFiles. |

### Molecule comparision

With 'get_shared_molecules_key' one get the number and identifier string for those molecules which are present in all of the compared datasets.  
  
| **function**  |  get_shared_molecules_key |   
|---|---|
|  **parameter**  | all_dict : *dict*  |
|   | &nbsp;&nbsp; Name of the dictionary with the imported SDFiles. |  
  
If there are no more than three datasets a Venn diagramm of the intersection of the datasets can be created using 'visualize_intersection'. The image will be saved in the output folder in a chosen format (variable: data_type).  
  
| **function**  |  visualize_intersection |   
|---|---|
|  **parameter**  | all_dict : *dict*  |
|   | &nbsp;&nbsp; Name of the dictionary with the imported SDFiles. |
|   | data_type : *str, default='png'*  |
|   | &nbsp;&nbsp; Data type for the exported image. |  

![intersection](https://user-images.githubusercontent.com/95417135/165087450-b3336c13-a3bf-4da5-98e1-97e2fa91167f.png)
|:--:| 
| *Example of an intersection from three compared datasets created with the 'visualize_intersection' function.* |

### Descriptor and descriptor value distribution 

The function 'get_descriptor_list_key' utilize a callable function (variable: descriptor) form RDKit to get descriptor values for every molecule. For example the callable functions can be from the rdkit.Chem.rdMolDescriptors module the rdkit.Chem.Descriptors module. The values are saved in the dataframe under a chosen name (variable: descriptor_list_keyname).  
  
| **function**  |  get_descriptor_list_key |   
|---|---|
|  **parameter**  | all_dict : *dict*  |
|   | &nbsp;&nbsp; Name of the dictionary with the imported SDFiles. |
|   | descriptor : *callable*  |
|   | &nbsp;&nbsp; RDKit function returning a molecular descriptor value.  | 
|   | descriptor_list_keyname : *str*  |
|   | &nbsp;&nbsp; Name for refering to descriptor values.  | 
  
To get and visualize the distribution of the descriptor values the function 'descriptor_counts_and_plot' can be used. The function distinguishes between continuous and discrete distributed descriptor vales and for the continuous values one can chose the binning size (variable: width_of_bins). The distribution of values is exported as csv-file and the visualization with a selectable format (variable: data_type) is also saved in the output folder.  
  
| **function**  | descriptor_counts_and_plot  |   
|---|---|
|  **parameter**  | all_dict : *dict*  |
|   | &nbsp;&nbsp; Name of the dictionary with the imported SDFiles. |
|   | descriptor_list_keyname : *str*  |
|   | &nbsp;&nbsp; Name refering to descriptor values which are to be visualized.  | 
|   | width_of_bins : *int, default=10*  |
|   | &nbsp;&nbsp; Intervall size for binning of continuous descriptor values.  |
|   | data_type : *str, default='png'*  |
|   | &nbsp;&nbsp; Data type for the exported image. |  
|   | save_dataframe : *bool, default=True*  |
|   | &nbsp;&nbsp; Option to export the counts of the descriptor values as csv. |  

![distribution_of_LogP](https://user-images.githubusercontent.com/95417135/165087902-6788db96-7230-4829-83f3-8fd3b6f791fa.png)
|:--:| 
| *Example of the a descriptor value distribution. The LogP values of the three datasets are binned in intervalls of 5.* |

  
With the database ID one can also search for the descriptor value of a specific molecule using the 'get_value_from_id' function. The function tells in which SDFile the molecule is found and the value of the descriptor.  
  
| **function**  | get_value_from_id  |   
|---|---|
|  **parameter**  | all_dict : *dict*  |
|   | &nbsp;&nbsp; Name of the dictionary with the imported SDFiles. |
|   | wanted_id : *str*  |
|   | &nbsp;&nbsp; Database ID from the molecule in question.  | 
|   | descriptor_list_keyname : *str*  |
|   | &nbsp;&nbsp; Name refering to descriptor values which are to be visualized.  |   

### Lipinski Rule of 5

With 'get_lipinski_key' the number of broken Lipinki Rules for every molecule is calculated and a summary for every SDFile is created.  
  
| **function**  | get_lipinski_key  |   
|---|---|
|  **parameter**  | all_dict : *dict*  |
|   | &nbsp;&nbsp; Name of the dictionary with the imported SDFiles. |  
  
Subsequently the 'lipinski_plot' function visualizes the number of broken rules. Again the results are exported as csv-file and the bar-plot is also saved with a selectable format (variable: data_type) in the output folder.  
  
| **function**  | lipinski_plot  |   
|---|---|
|  **parameter**  | all_dict : *dict*  |
|   | &nbsp;&nbsp; Name of the dictionary with the imported SDFiles. |
|   | data_type : *str, default='png'*  |
|   | &nbsp;&nbsp; Data type for the exported image. |  
|   | save_dataframe : *bool, default=True*  |
|   | &nbsp;&nbsp; Option to export the counts of the descriptor values as csv. | 

![lipinski_rules_plot](https://user-images.githubusercontent.com/95417135/165089273-82852c1d-9e50-41d3-86f1-043c30d5ebf1.png)
| *Example of a Lipinski Plot with three datasets.* |
  
### Chemical Space Visualization

For the visualization of the chemical space the chemplot module is used. The extended connectivity fingerprints can be specified with the fingerprint radius (variable: fp_radius) and the size (variable: fp_bits). For the dimension reduction PCA, t-SNE or UMAP can be chosen (variable: dimension_reduction). The chemical space plot is saved in the output folder except when choosing to create an interactive plot (variable: interactive). Then the plot will be displayed and can be manually saved.  
  
| **function**  | lipinski_plot  |   
|---|---|
|  **parameter**  | all_dict : *dict*  |
|   | &nbsp;&nbsp; Name of the dictionary with the imported SDFiles. |
|   | fp_radius : *int, default=2*  |
|   | &nbsp;&nbsp; Radius of the Extended Connectivity Fingerprints. |  
|   | fp_bits : *int, default=2048*  |
|   | &nbsp;&nbsp; Size of the Extended Connetivity Fingerprints. |  
|   | dimension_reduction : *str, default='pca'*  |
|   | &nbsp;&nbsp; Method of dimension reduction ("pca", "umap" or "tsne"). |  
|   | interactive : *bool, default=True*  |
|   | &nbsp;&nbsp; Option to create an interactive plot. |  


![chemical_space](https://user-images.githubusercontent.com/95417135/165089065-dea082d6-2600-41dd-8c57-6e724a25f474.png)
|:--:| 
| *Example of the chemical space visualization with three datasets (interactive=False).* |
