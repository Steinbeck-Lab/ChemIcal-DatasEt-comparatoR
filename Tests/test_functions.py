from CIDER import ChemicalDatasetComparator
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import IPythonConsole
import os
import pytest

cider = ChemicalDatasetComparator()
testdict = cider.import_as_data_dict('unittest_data')

def test_import_as_data_dict():
   testdict = cider.import_as_data_dict('unittest_data')
   # Assert that the function generates the dictionary
   assert list(testdict.keys()) == ['set_A.sdf', 'set_B.sdf', 'set_D.sdf']

def test_get_number_of_molecules():
   cider.get_number_of_molecules(testdict)
   # Assert that the function generates new entries in the dictionary
   # and that the correct number of molecules are found in the datasets
   assert testdict['set_A.sdf']['number_of_molecules'] == 3
   assert testdict['set_B.sdf']['number_of_molecules'] == 4
   assert testdict['set_D.sdf']['number_of_molecules'] == 6

def test_draw_molecules():
    cider.draw_molecules(testdict, number_of_mols = 3)
    # Assert that the function gernerates a new entry in the dictionary
    assert any(key == 'molecule_picture' for key in list(testdict['set_A.sdf'].keys())) == True
    # Assert that the pictures are exported
    assert os.path.exists("output/mol_grit_set_A.png") 

def test_get_database_id():
    cider.get_database_id(testdict, 'coconut_id')
    # Assert that the function gernerates a new entry in the dictionary
    assert any(key == 'coconut_id_keyname' for key in list(testdict['set_A.sdf'].keys())) == True
    # Assert that the function gets the correct IDs 
    assert list(testdict['set_A.sdf']['coconut_id_keyname']) == ['CNP0206286', 'CNP0284887', 'CNP0080171']

def test_get_identifier_list_key():
    cider.get_identifier_list_key(testdict)
    # Assert that the function gernerates a new entry in the dictionary
    assert any(key == 'identifier_list' for key in list(testdict['set_A.sdf'].keys())) == True
    # Assert that at default configuration the correct InChI string is generated.
    assert testdict['set_A.sdf']['identifier_list'][0] == 'InChI=1S/C6H3Cl3/c7-4-1-2-5(8)6(9)3-4/h1-3H'
    # Assert that with id_type = 'smiles' the correct smiles string is generated.
    cider.get_identifier_list_key(testdict, 'smiles')
    assert testdict['set_A.sdf']['identifier_list'][0] == 'Clc1ccc(Cl)c(Cl)c1'
    # Assert that with id_type = 'inchikey' the correct InChIKey is generated.
    cider.get_identifier_list_key(testdict, 'inchikey')
    assert testdict['set_A.sdf']['identifier_list'][0] == 'PBKONEOXTCPAFI-UHFFFAOYSA-N'

def test_get_duplicate_key():
    cider.get_duplicate_key(testdict)
    # Assert that the function generates a new entry in the dictionary
    assert any(key == 'number_of_duplicates' for key in list(testdict['set_A.sdf'].keys())) == True
    # Assert that the function returns the right number of duplicates
    assert testdict['set_A.sdf']['number_of_duplicates'] == 0
    assert testdict['set_B.sdf']['number_of_duplicates'] == 0
    assert testdict['set_D.sdf']['number_of_duplicates'] == 1
   
def test_get_shared_molecules_key():
    cider.get_shared_molecules_key(testdict)
    # Assert that the function generates a new entry in the dictionary
    assert any(key == 'number_of_shared_molecules' for key in list(testdict['set_A.sdf'].keys())) == True
    assert any(key == 'shared_molecules' for key in list(testdict['set_A.sdf'].keys())) == True
    # Assert that the correct number and identity of the shared molecules are returned 
    assert testdict['set_A.sdf']['number_of_shared_molecules'] == 1
    assert testdict['set_A.sdf']['shared_molecules'] == 'ZPQOPVIELGIULI-UHFFFAOYSA-N' or 'InChI=1S/C6H4Cl2/c7-5-2-1-3-6(8)4-5/h1-4H' or 'Clc1cccc(Cl)c1'

def test_visualize_intersection():
    cider.visualize_intersection(testdict)
    assert os.path.exists("output/intersection.png") 

def test_get_descriptor_list_key():
    # Assertion for continuous values
    cider.get_descriptor_list_key(testdict,  Descriptors.MolWt, 'Molecular Weight')
    # Assert that the function generates a new entry in the dictionary
    assert any(key == 'Molecular Weight' for key in list(testdict['set_A.sdf'].keys())) == True
    # Assert that the correct values for the molecular weight are returned from rdkit.Chem Descriptors
    assert list(testdict['set_A.sdf']['Molecular Weight']) == [181.449, 147.00399999999996, 147.004]
    # Assert that the picture and the table are exported
    assert os.path.exists("output/distribution_of_Molecular Weight.png") 
    assert os.path.exists("output/table_Molecular Weight.csv")
    # Assertion for discrete values
    cider.get_descriptor_list_key(testdict, rdMolDescriptors.CalcMolFormula, 'Molecular Formula')
    # Assert that the function generates a new entries in the dictionary
    assert any(key == 'Molecular Formula' for key in list(testdict['set_A.sdf'].keys())) == True
    # Assert that the correct strings for the molecular formula are returned from rdkit.Chem rdMolDescriptors
    assert list(testdict['set_A.sdf']['Molecular Formula']) == ['C6H3Cl3', 'C6H4Cl2', 'C6H4Cl2']
    # Assert that the picture and the table are exported
    assert os.path.exists("output/distribution_of_Number of Rings.png") 
    assert os.path.exists("output/table_Number of Rings.csv")

def test_get_value_from_id(capfd):
    cider.get_value_from_id(testdict, 'CNP0206286', 'Molecular Weight')
    out, err = capfd.readouterr()
    assert out == 'Molecule found in set_A.sdf\nMolecular Weight value for ID CNP0206286: 181.449\nMolecule not found in set_B.sdf\nMolecule not found in set_D.sdf\n'

def test_descriptor_counts_and_plot():
    cider.descriptor_counts_and_plot(testdict, 'Molecular Weight', 5)
    # Assert that the function generates a new entry in the dictionary for binned continuous values
    assert any(key == 'binned Molecular Weight' for key in list(testdict['set_A.sdf'].keys())) == True
    cider.get_descriptor_list_key(testdict, Descriptors.RingCount, 'Number of Rings')
    cider.descriptor_counts_and_plot(testdict, 'Number of Rings')
    # Assert that the function generates a new entry in the dictionary for binned discrete values
    assert any(key == 'binned Number of Rings' for key in list(testdict['set_A.sdf'].keys())) == True

def test_get_lipinski_key():
    cider.get_lipinski_key(testdict)
    # Assert that the function generates new entries in the dictionary
    assert any(key == 'number_of_broken_Lipinski_Rules' for key in list(testdict['set_A.sdf'].keys())) == True
    assert any(key == 'Lipinski_Rule_of_5_summary' for key in list(testdict['set_A.sdf'].keys())) == True
    # Assert that the correct numbers are returned
    assert testdict['set_A.sdf']['number_of_broken_Lipinski_Rules'] == [0,0,0]
    assert testdict['set_A.sdf']['Lipinski_Rule_of_5_summary'] ==  {'lipinski_molecules': 3,
                                                                    '1_rule_broken': 0,
                                                                    '2_rules_broken': 0,
                                                                    '3_rules_broken': 0,
                                                                    '4_rules_broken': 0}

def test_lipinski_plot():
   cider.lipinski_plot(testdict)
   # Assert that the picture is exported
   assert os.path.exists("output/lipinski_rules_plot.png")

def test_chemical_space_visualization():
   cider.chemical_space_visualization(testdict, interactive = False)
   # Assert that the picture is exported
   assert os.path.exists("output/chemical_space.png")

def test_export_single_dict_values():
   cider.export_single_dict_values(testdict)
   # Assert that the picture is exported
   assert os.path.exists("output/descriptor_values_set_A.csv")
   assert os.path.exists("output/descriptor_values_set_B.csv")
   assert os.path.exists("output/descriptor_values_set_D.csv")

def test_export_all_picture_pdf():
   cider.export_all_picture_pdf()
   # Assert that the picture is exported
   assert os.path.exists("output/all_images.pdf")
