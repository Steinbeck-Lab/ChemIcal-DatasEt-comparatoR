from CIDER import ChemicalDatasetComparator
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
import os

cider = ChemicalDatasetComparator()
test_dict_path = os.path.join(os.path.split(__file__)[0],
                              "unittest_data")
testdict = cider.import_as_data_dict(test_dict_path)


def test_import_as_data_dict():
    # Assert that the function generates the dictionary
    assert set(testdict.keys()) == {"set_A.sdf", "set_B.sdf", "set_D.sdf", cider.figure_dict_keyname}

def test_check_invalid_mols_in_SDF():
    test_dict_path = os.path.join(os.path.split(__file__)[0], "unittest_data_invalid")
    invalid_testdict = cider.import_as_data_dict(test_dict_path)
    # Assert that faulty molecules are detected and optionally deleted
    cider._check_invalid_mols_in_SDF(invalid_testdict, delete=False)
    assert len(invalid_testdict["set_D_invalid.sdf"][cider.import_keyname]) == 7
    cider._check_invalid_mols_in_SDF(invalid_testdict, delete=True)
    assert len(invalid_testdict["set_D_invalid.sdf"][cider.import_keyname]) == 5

def test_get_number_of_molecules():
    cider.get_number_of_molecules(testdict)
    # Assert that the function generates new entries in the dictionary
    # and that the correct number of molecules are found in the datasets
    assert testdict["set_A.sdf"][cider.dataset_length_keyname] == 4
    assert testdict["set_B.sdf"][cider.dataset_length_keyname] == 4
    assert testdict["set_D.sdf"][cider.dataset_length_keyname] == 7

def test_draw_molecules():
    cider.draw_molecules(testdict)
    # Assert that the picture is exported
    test_path = os.path.join(os.path.split(__file__)[0],
                             "output", "mol_grid.png")
    assert os.path.exists(test_path)
    # Assert that a new entry in self.figure_dict_keyname is generated
    assert (
        any(key == "mol_grid" for key in list(testdict[cider.figure_dict_keyname].keys()))
    )

def test_get_database_id():
    cider.get_database_id(testdict, "coconut_id")
    # Assert that the function generates a new entry in the dictionary
    assert (
        any(key == "coconut_id" for key in list(testdict["set_A.sdf"].keys()))
    )
    # Assert that the function gets the correct IDs
    assert list(testdict["set_A.sdf"][cider.database_id_keyname]) == [
        "CNP0206286",
        "CNP0284887",
        "CNP0080171",
        "CNP0100364",
    ]


def test_get_identifier_list():
    testset = testdict["set_A.sdf"][cider.import_keyname]
    # Assert that the correct InChI strings are generated.
    expected_inchi = (
        [
            "InChI=1S/C6H3Cl3/c7-4-1-2-5(8)6(9)3-4/h1-3H",
            "InChI=1S/C6H4Cl2/c7-5-3-1-2-4-6(5)8/h1-4H",
            "InChI=1S/C6H4Cl2/c7-5-2-1-3-6(8)4-5/h1-4H",
            "InChI=1S/C9H4BrClO2/c10-7-4-13-8-2-1-5(11)3-6(8)9(7)12/h1-4H",
        ],
        0,
    )
    inchi = cider._get_identifier_list(testset)
    assert expected_inchi == inchi
    # Assert that the correct InChIKeys are generated.
    expected_inchikey = (
        [
            "PBKONEOXTCPAFI-UHFFFAOYSA-N",
            "RFFLAFLAYFXFSW-UHFFFAOYSA-N",
            "ZPQOPVIELGIULI-UHFFFAOYSA-N",
            "ZILSBPAUJJBEFF-UHFFFAOYSA-N",
        ],
        0,
    )
    inchikey = cider._get_identifier_list(testset, "inchikey")
    assert expected_inchikey == inchikey
    # Assert that the correct SMILES strings are generated.
    expected_smiles = (
        [
            "Clc1ccc(Cl)c(Cl)c1",
            "Clc1ccccc1Cl",
            "Clc1cccc(Cl)c1",
            "O=c1c(Br)coc2ccc(Cl)cc12",
        ],
        0,
    )
    smiles = cider._get_identifier_list(testset, "smiles")
    assert expected_smiles == smiles

def test_get_identifier_list_key():
    cider.get_identifier_list_key(testdict)
    # Assert that the function generates a new entry in the dictionary
    assert (
        any(
            key == cider.identifier_keyname
            for key in list(testdict["set_A.sdf"].keys())
        )
    )
    # Assert that at default configuration the correct InChI string is generated.
    assert (
        testdict["set_A.sdf"][cider.identifier_keyname][0]
        == "InChI=1S/C6H3Cl3/c7-4-1-2-5(8)6(9)3-4/h1-3H"
    )
    # Assert that with id_type = 'smiles' the correct smiles string is generated.
    cider.get_identifier_list_key(testdict, "smiles")
    assert testdict["set_A.sdf"][cider.identifier_keyname][0] == "Clc1ccc(Cl)c(Cl)c1"
    # Assert that with id_type = 'inchikey' the correct InChIKey is generated.
    cider.get_identifier_list_key(testdict, "inchikey")
    assert (
        testdict["set_A.sdf"][cider.identifier_keyname][0]
        == "PBKONEOXTCPAFI-UHFFFAOYSA-N"
    )


def test_get_duplicate_key():
    cider.get_duplicate_key(testdict)
    # Assert that the function generates a new entry in the dictionary
    assert (
        any(
            key == cider.duplicates_keyname
            for key in list(testdict["set_A.sdf"].keys())
        )
    )
    # Assert that the function returns the right number of duplicates
    assert testdict["set_A.sdf"][cider.duplicates_keyname] == 0
    assert testdict["set_B.sdf"][cider.duplicates_keyname] == 0
    assert testdict["set_D.sdf"][cider.duplicates_keyname] == 1

def test_get_shared_molecules_key():
    cider.get_shared_molecules_key(testdict)
    # Assert that the function generates a new entry in the dictionary
    assert (
        any(
            key == cider.shared_mols_keyname
            for key in list(testdict["set_A.sdf"].keys())
        )
    )
    assert (
        any(
            key == cider.shared_mols_id_keyname
            for key in list(testdict["set_A.sdf"].keys())
        )
    )
    # Assert that the correct number and identity of the shared molecules are returned
    assert testdict["set_A.sdf"][cider.shared_mols_keyname] == 1
    assert (
        testdict["set_A.sdf"][cider.shared_mols_id_keyname]
        == "ZPQOPVIELGIULI-UHFFFAOYSA-N"
        or "InChI=1S/C6H4Cl2/c7-5-2-1-3-6(8)4-5/h1-4H"
        or "Clc1cccc(Cl)c1"
    )

def test_visualize_intersection():
    cider.visualize_intersection(testdict)
    # Assert that the figure is exported
    test_path = os.path.join(os.path.split(__file__)[0],
                             "output", "intersection.png")
    assert os.path.exists(test_path)
    # Assert that a new entry in self.figure_dict_keyname is generated
    assert (
        any(key == "intersection" for key in list(testdict[cider.figure_dict_keyname].keys()))
    )


def test_get_descriptor_list():
    testset = testdict["set_A.sdf"][cider.import_keyname]
    # Assert that the correct LogP values are returned
    expected_LogP = [3.6468000000000007, 2.9934000000000003, 2.993400000000001, 3.2089000000000008]
    LogP = cider._get_descriptor_list(testset, Descriptors.MolLogP)
    assert expected_LogP == LogP

def test_get_descriptor_list_key_1():
    # Assertion for continuous values
    cider.get_descriptor_list_key(testdict, Descriptors.MolWt, "Molecular Weight")
    # Assert that the function generates a new entry in the dictionary
    assert (
        any(key == "Molecular Weight" for key in list(testdict["set_A.sdf"].keys()))
    )
    # Assert that the correct values for the molecular weight are returned from rdkit.Chem Descriptors
    assert list(testdict["set_A.sdf"]["Molecular Weight"]) == [
        181.449,
        147.00399999999996,
        147.004,
        259.48600000000005,
    ]

def test_get_descriptor_list_key_2():
    # Assertion for discrete values
    cider.get_descriptor_list_key(
        testdict, rdMolDescriptors.CalcMolFormula, "Molecular Formula"
    )
    # Assert that the function generates a new entry in the dictionary
    assert (
        any(key == "Molecular Formula" for key in list(testdict["set_A.sdf"].keys()))
    )
    # Assert that the correct strings for the molecular formula are returned from rdkit.Chem rdMolDescriptors
    assert list(testdict["set_A.sdf"]["Molecular Formula"]) == [
        "C6H3Cl3",
        "C6H4Cl2",
        "C6H4Cl2",
        "C9H4BrClO2",
    ]

def test_get_discrete_descriptor_counts():
    cider.get_descriptor_list_key(
        testdict, Descriptors.NumHDonors, "Number of H-Donors"
    )
    cider._get_discrete_descriptor_counts(testdict, "Number of H-Donors")
    # Assert that the function generates a new entry in the dictionary
    assert any(
        key == "binned_Number of H-Donors" for key in list(testdict["set_A.sdf"].keys())
    )
    # Assert the correct values for the bins
    assert list(testdict["set_A.sdf"]["binned_Number of H-Donors"]) == [4, 0]

def test_get_continuous_descriptor_counts():
    cider.get_descriptor_list_key(testdict, Descriptors.MolLogP, "LogP")
    cider._get_continuous_descriptor_counts(testdict, "LogP", 2)
    # Assert that the function generates a new entry in the dictionary
    assert any(key == "binned_LogP" for key in list(testdict["set_A.sdf"].keys()))
    # Assert the correct values for the bins
    assert list(testdict["set_A.sdf"]["binned_LogP"]) == [0, 4, 0, 0, 0]

def test_discrete_descriptor_plot():
    cider._discrete_descriptor_plot(testdict, "Number of H-Donors")
    # Assert that the image and the table are exported
    test_path_NumHDpng = os.path.join(os.path.split(__file__)[0],
                                      "output", "distribution_of_Number of H-Donors.png")
    test_path_NumHDcsv = os.path.join(os.path.split(__file__)[0],
                                      "output", "table_Number of H-Donors.csv")
    assert os.path.exists(test_path_NumHDpng)
    assert os.path.exists(test_path_NumHDcsv)
    # Assert that a new entry in self.figure_dict_keyname is generated
    assert (
        any(key == "distribution_of_Number of H-Donors" for key in list(testdict[cider.figure_dict_keyname].keys()))
    )

def test_continuous_descriptor_plot():
    cider._continuous_descriptor_plot(testdict, "LogP")
    # Assert that the image and the table are exported
    test_path_logppng = os.path.join(os.path.split(__file__)[0],
                                     "output", "distribution_of_LogP.png")
    test_path_logpcsv = os.path.join(os.path.split(__file__)[0],
                                     "output", "table_LogP.csv")
    assert os.path.exists(test_path_logppng)
    assert os.path.exists(test_path_logpcsv)
    # Assert that a new entry in self.figure_dict_keyname is generated
    assert (
        any(key == "distribution_of_LogP" for key in list(testdict[cider.figure_dict_keyname].keys()))
    )

def test_descriptor_counts_and_plot():
    # Assertion for continuous descriptor
    cider.descriptor_counts_and_plot(testdict, "Molecular Weight", 10)
    # Assert that the function generates a new entry in the dictionary for binned continuous values
    assert any(
        key == "binned_Molecular Weight" for key in list(testdict["set_A.sdf"].keys())
    )
    # Assert the correct values for the bins
    assert list(testdict["set_A.sdf"]["binned_Molecular Weight"]) == [
        2,
        0,
        0,
        0,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        1,
    ]
    # Assert that a new entry in self.figure_dict_keyname is generated
    assert (
        any(key == "distribution_of_Molecular Weight" for key in list(testdict[cider.figure_dict_keyname].keys()))
    )

    # Assertion for discrete descriptor
    cider.get_descriptor_list_key(testdict, Descriptors.RingCount, "Number of Rings")
    cider.descriptor_counts_and_plot(testdict, "Number of Rings", 2)
    # Assert that the function generates a new entry in the dictionary for binned discrete values
    assert any(
        key == "binned_Number of Rings" for key in list(testdict["set_A.sdf"].keys())
    )
    # Assert the correct values for the bins
    assert list(testdict["set_A.sdf"]["binned_Number of Rings"]) == [0, 3, 1]
    # Assert that a new entry in self.figure_dict_keyname is generated
    assert (
        any(key == "distribution_of_Number of Rings" for key in list(testdict[cider.figure_dict_keyname].keys()))
    )

def test_test_for_lipinski():
    testset = testdict["set_A.sdf"][cider.import_keyname]
    # Assert the correct number of broken rules are returned
    expected_num_of_break = [0, 0, 0, 0]
    num_of_break = cider._test_for_lipinski(testset)
    assert expected_num_of_break == num_of_break

def test_get_lipinski_key():
    cider.get_lipinski_key(testdict)
    # Assert that the function generates new entries in the dictionary
    assert (
        any(
            key == cider.lipinski_list_keyname
            for key in list(testdict["set_A.sdf"].keys())
        )
    )
    assert (
        any(
            key == cider.lipinski_summary_keyname
            for key in list(testdict["set_A.sdf"].keys())
        )
    )
    # Assert that the correct numbers are returned
    assert testdict["set_A.sdf"][cider.lipinski_list_keyname] == [0, 0, 0, 0]
    assert testdict["set_A.sdf"][cider.lipinski_summary_keyname] == {
        "lipinski_molecules": 4,
        "1_rule_broken": 0,
        "2_rules_broken": 0,
        "3_rules_broken": 0,
        "4_rules_broken": 0,
    }

def test_lipinski_plot():
    cider.lipinski_plot(testdict)
    # Assert that the picture is exported
    test_path = os.path.join(os.path.split(__file__)[0],
                             "output", "lipinski_plot.png")
    assert os.path.exists(test_path)
    # Assert that a new entry in self.figure_dict_keyname is generated
    assert (
        any(key == "lipinski_plot" for key in list(testdict[cider.figure_dict_keyname].keys()))
    )

def test_get_Murcko_scaffold():
    testset = testdict["set_A.sdf"][cider.import_keyname]
    scaffolds = cider._get_Murcko_scaffold(testset)
    # Assert that the correct scaffold list is generated (for default)
    expected_scaffold_list = ['c1ccccc1', 'c1ccccc1', 'c1ccccc1', 'O=c1ccoc2ccccc12']
    assert scaffolds[1] == expected_scaffold_list
    # Assert that the correct scaffold counts are calculated (for default)
    assert list(scaffolds[2]) == [0.75, 0.25]

def test_get_Murcko_scaffold_2():
    testset = testdict["set_A.sdf"][cider.import_keyname]
    scaffolds = cider._get_Murcko_scaffold(testset, framework=True)
    # Assert that the correct scaffold list is generated (for framework=True)
    expected_scaffold_list = ['c1ccccc1', 'c1ccccc1', 'c1ccccc1', 'C1=COc2ccccc2C1']
    assert scaffolds[1] == expected_scaffold_list
    # Assert that the correct scaffold counts are calculated (for framework=True)
    assert list(scaffolds[2]) == [0.75, 0.25]

def test_get_Murcko_scaffold_3():
    testset = testdict["set_A.sdf"][cider.import_keyname]
    scaffolds = cider._get_Murcko_scaffold(testset, skeleton=True)
    # Assert that the correct scaffold list is generated (for skeleton=True)
    expected_scaffold_list = ['C1CCCCC1', 'C1CCCCC1', 'C1CCCCC1', 'C1CCC2CCCCC2C1']
    assert scaffolds[1] == expected_scaffold_list
    # Assert that the correct scaffold counts are calculated (for skeleton=True)
    assert list(scaffolds[2]) == [0.75, 0.25]

def test_draw_most_frequent_scaffolds():
    cider.draw_most_frequent_scaffolds(testdict)
    # Assert that the function generates new keys in the dictionary
    assert (
        any(
            key == cider.scaffold_list_keyname
            for key in list(testdict["set_A.sdf"].keys())
        )
    )
    assert (
        any(
            key == cider.scaffold_summary_keyname
            for key in list(testdict["set_A.sdf"].keys())
        )
    )
    # Assert that the picture is exported
    test_path = os.path.join(os.path.split(__file__)[0],
                             "output", "scaffold_grid.png")
    assert os.path.exists(test_path)
    # Assert that a new entry in self.figure_dict_keyname is generated
    assert (
        any(key == "scaffold_grid" for key in list(testdict[cider.figure_dict_keyname].keys()))
    )

def test_chemical_space_visualization():
    cider.chemical_space_visualization(testdict, interactive=False)
    # Assert that the picture is exported
    test_path = os.path.join(os.path.split(__file__)[0],
                             "output", "chemical_space.png")
    assert os.path.exists(test_path)

def test_export_single_dict_values():
    cider.export_single_dict_values(testdict)
    # Assert that the picture is exported
    test_path_setA = os.path.join(os.path.split(__file__)[0],
                                  "output", "descriptor_values_set_A.csv")
    test_path_setB = os.path.join(os.path.split(__file__)[0],
                                  "output", "descriptor_values_set_B.csv")
    test_path_setD = os.path.join(os.path.split(__file__)[0],
                                  "output", "descriptor_values_set_D.csv")
    assert os.path.exists(test_path_setA)
    assert os.path.exists(test_path_setB)
    assert os.path.exists(test_path_setD)

def test_export_all_figures_pdf():
    cider.export_all_figures_pdf(testdict)
    # Assert that the picture is exported
    test_path = os.path.join(os.path.split(__file__)[0],
                             "output", "all_figures.pdf")
    assert os.path.exists(test_path)

def test_chemical_space_plot_tsne():
    cider.chemical_space_visualization(testdict,
                                       dimension_reduction = 'tsne',
                                       interactive = False,
                                       fp_bits = 1024,
                                       fp_radius = 2)
