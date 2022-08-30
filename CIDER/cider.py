"""MIT License

Copyright (c) 2022 Hannah Busch, Jonas Schaub, Otto Brinkhaus, Kohulan Rajan
and Christoph Steinbeck


Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE."""

# Section: Import Libraries

import os
import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Scaffolds import MurckoScaffold

import matplotlib.pyplot as plt

from matplotlib_venn import venn2
from matplotlib_venn import venn3

import chemplot as cp
from chemplot import descriptors

from fpdf import FPDF

from typing import List, Tuple, Dict
from itertools import count

class ChemicalDatasetComparator:
    """
    Wrapper class around all Cider functionalities

    ChemIcal DatasEt comparatoR (CIDER) is a Python package and ready-to-use Jupyter Notebook workflow which primarily utilizes RDKit to compare two or more chemical structure datasets (SD files) in different aspects, e.g. size, overlap, molecular descriptor distributions, chemical space clustering, etc., most of which can be visually inspected.
    """

    def __init__(self):
        """
        The class variables of CIDER function as keys for the dictionary in which all calculated and plotted data will be stored.
        """
        self.import_keyname = "SDMolSupplier_Object"
        self.dataset_length_keyname = "number_of_molecules"
        self.identifier_keyname = "identifier_list"
        self.duplicates_keyname = "number_of_duplicates"
        self.duplicates_id_keyname = "duplicates"
        self.shared_mols_keyname = "number_of_shared_molecules"
        self.shared_mols_id_keyname = "shared_molecules"
        self.lipinski_list_keyname = "number_of_broken_lipinski_rules"
        self.lipinski_summary_keyname = "lipinski_summary"
        self.database_id_keyname = "coconut_id"
        self.scaffold_list_keyname = "scaffold_list"
        self.scaffold_summary_keyname = "scaffold_summary"

    # Section: Import data and check for faulty SDFiles

    def _check_invalid_mols_in_SDF(self, all_dicts: dict, delete: bool = False) -> None:
        """
        This function checks if there are invalid entries in the SDFiles that can cause errors in the subsequent functions. At choice the invalid entries can be removed form the SDMolSupplier Object. The entry will remain in the original SDFile as it is. (private method)

        Args:
            all_dicts (dict): Dictionary with sub-dictionaries including SDMolSupplier Objects.
            delete (bool): Deleting invalid entries or not (default: False).
        """
        for single_dict in all_dicts:
            mol_index = -1
            invalid_index = []
            for mol in all_dicts[single_dict][self.import_keyname]:
                mol_index += 1
                if not mol:
                    print(
                        str(single_dict)
                        + " has invalid molecule at index "
                        + str(mol_index)
                    )
                    invalid_index.append(mol_index)
            if not invalid_index:
                print("No invalid molecules found in " + str(single_dict))
            elif delete and invalid_index:
                new_SDMol = list(all_dicts[single_dict][self.import_keyname])
                for index in sorted(invalid_index, reverse=True):
                    del new_SDMol[index]
                all_dicts[single_dict].update({self.import_keyname: new_SDMol})
                print(
                    str(len(invalid_index))
                    + " invalid molecule(s) are deleted from "
                    + str(single_dict)
                )
            elif not delete and invalid_index:
                print(
                    str(len(invalid_index))
                    + " invalid molecule(s) will remain in "
                    + str(single_dict)
                )
        return

    def import_as_data_dict(self, path_to_data: str, delete: bool = False) -> Dict:
        """
        This function creates a dictionary with the names of the imported file as keys. The values of each of these keys is a subdictionary. The first entry of every subdictionary is self.import_keyname (class variable, can be changed) as key and a SDMolSupplier Object of the SDFile as value. To find faulty molecules every entry of the SDMolSupplier Object will be parsed once. (Parsed molecules will not be stored in the dictionary to save memory.)

        Args:
            path_to_data (str): Path to the directory where the SDFiles are stored.
            delete (bool): Delete faulty molecules from imported dataset.

        Returns:
            all_dicts (dict): Dictionary with subdictionaries for every dataset, updated with the SDMolSupplier Objects.

        Raises:
            FileNotFoundError: if the data path is invalid
            Statement: if no SDFiles are found
        """
        all_dicts = {}
        data_dir = os.path.normpath(str(path_to_data))
        for dict_name in os.listdir(data_dir):
            if dict_name[-3:] == "sdf":
                single_dict = {}
                dict_path = os.path.join(data_dir, dict_name)
                single_dict[self.import_keyname] = Chem.SDMolSupplier(dict_path)
                all_dicts[dict_name] = single_dict
        if not all_dicts:
            print(
                "No SDFiles in "
                + str(data_dir)
                + " found!"
            )
        self._check_invalid_mols_in_SDF(all_dicts, delete)
        return all_dicts

    # Section: Get overview of the dataset size and molecules

    def get_number_of_molecules(self, all_dicts: dict) -> None:
        """
        This function updates the subdictionaries in the given dictionary (created from import_as_data_dict function) with the number of molecules in every dataset as new key-value pair. The key is the class variable 'cider.dataset_length_keyname'.

        Args:
            all_dicts (dict): Dictionary with subdictionaries with SDMolSupplier Objects.
        """
        for single_dict in all_dicts:
            number_of_molecules = len(all_dicts[single_dict][self.import_keyname])
            all_dicts[single_dict][self.dataset_length_keyname] = number_of_molecules
            print(
                "Number of molecules in "
                + single_dict
                + ": "
                + str(all_dicts[single_dict][self.dataset_length_keyname])
            )
        print("Updated dictionary with '" + self.dataset_length_keyname + "'")
        return

    def draw_molecules(
        self,
        all_dicts: dict,
        number_of_mols: int = 12,
        mols_per_row: int = 3,
        image_size: int = 200,
        data_type: str = "png",
        figsize: Tuple[float, float] = [20.0, 20.0],
        fontsize_title: int = 24,
        fontsize_subtitle: int = 20
    ):
        """
        This function creates an grid image of the first molecules of each dataset and exports the image to an output folder.

        Args:
            all_dicts (dict): dictionary with subdictionaries with SDMolSupplier Objects.
            number_of_mols (int): number of molecules form each dataset that will be displayed (default: 12).
            mols_per_row (int): number of molecules per row in the grid (default: 3).
            image_size (int): the size of the image for a single molecule (defaut: 200).
            data_type (str): data type for the exported files (e.g. png, jpg, pdf, default: png).
            figsize (float, float): Width, height of the image in inches (default: 20, 20)
            fontsize_title (int): Fontsize of the title (default: 24).
            fontsize_subtitle (int): Fontsize of the subtitles (default: 20).

        Returns:
            fig (matplotlib.figure): grid image of molecules
        """
        image_list = []
        title_list = []
        for single_dict in all_dicts:
            title_list.append(single_dict)
            to_draw = []
            if len(all_dicts[single_dict][self.import_keyname]) < number_of_mols:
                number_of_mols = len(all_dicts[single_dict][self.import_keyname])
            for i in range(number_of_mols):
                to_draw.append(all_dicts[single_dict][self.import_keyname][i])
            mol_grid = Draw.MolsToGridImage(
                to_draw,
                maxMols=number_of_mols,
                molsPerRow=mols_per_row,
                subImgSize=(image_size, image_size),
                returnPNG=False,
            )
            image_list.append(mol_grid)
        rows = len(image_list)
        fig = plt.figure(figsize=figsize)
        for j in range(0, rows):
            fig.add_subplot(rows, 1, j + 1)
            plt.axis("off")
            plt.imshow(image_list[j])
            plt.title(title_list[j], fontsize=fontsize_subtitle)
        fig.suptitle("Exemplary molecules from the datasets", fontsize=fontsize_title)
        if not os.path.exists("output"):
            os.makedirs("output")
        fig.savefig("output/mol_grid.%s" % (data_type))
        plt.close(fig)
        return fig

    # Section: Get database ID

    def get_database_id(self, all_dicts: dict, id_name: str) -> None:
        """
        This function updates subdictionaries of a given dictionary with a list of database IDs for the single molecules as new key-value pairs. Depending on which database the molecules are coming from, the key as a class variable can be changed accordingly.
        (To get the database ID the SDMolSupplier Objects needs to be parsed, this may take same time because no parsed molecules are saved in the dictionary to save memory.)

        Args:
            all_dicts (dict): dictionary with subdictionaries including SDMolSupplier Objects (self.import_keyname).
            id_name (str): ID name in the original SDFile.
        """
        for single_dict in all_dicts:
            database_id_list = []
            for mol in all_dicts[single_dict][self.import_keyname]:
                prop_dict = mol.GetPropsAsDict()
                database_id = prop_dict.get(id_name)
                database_id_list.append(database_id)
                all_dicts[single_dict][self.database_id_keyname] = database_id_list
        print("Updated dictionary with '" + self.database_id_keyname + "'")
        return

    # Section: Get string identifier

    def _get_identifier_list(
        self, moleculeset: Chem.SDMolSupplier, id_type: str = "inchi"
    ) -> Tuple[list, count]:
        """
        This function returns a list of InChI, InChIKey or canonical SMILES strings for all molecules in a given SDMolSupplier object. (private method)

        Args:
            moleculeset (rdkit.Chem.SDMolSupplier):
            id_type (str, optional): "inchi", "inchikey" or "smiles". Defaults to "inchi".

        Raises:
            ValueError: if ID_type is not "inchi," "inchikey" or "smiles".

        Returns:
            list[str]: List of identifiers based on given molecules.
            int: Counter of molecules for which no identifier could be determined
        """
        identifier_list = []
        failed_identifier_counter = 0
        for mol in moleculeset:
            if not mol:
                identifier = "Failed"
                failed_identifier_counter += 1
            elif id_type == "smiles":
                identifier = Chem.MolToSmiles(mol)
            elif id_type == "inchikey":
                identifier = Chem.MolToInchiKey(mol)
            elif id_type == "inchi":
                identifier = Chem.MolToInchi(mol)
            else:
                raise ValueError(
                    'id_type argument needs to be "smiles", "inchikey" or "inchi"!'
                )
            identifier_list.append(identifier)
        return identifier_list, failed_identifier_counter

    def get_identifier_list_key(self, all_dicts: dict, id_type: str = "inchi") -> None:
        """
        This function updates the subdictionaries in the given dictionary (created with the import_as_data_dict function) with a list of identifiers (InChI, InChIKey, canonical SMILES strings) as a new key-value pair using  _get_identifier_list on the SDMolSupplier Objects. The key self.identifier_keyname (class variable) can be changed.

        Args:
            all_dicts (dict): Dictionary with subdictionaries including SDMolSupplier Objects (self.import_keyname).
            id_type (str): Type of Identifier ("inchi", "inchikey" or "smiles")
        """
        for single_dict in all_dicts:
            identifier_list = self._get_identifier_list(
                all_dicts[single_dict][self.import_keyname], id_type
            )
            all_dicts[single_dict][self.identifier_keyname] = identifier_list[0]
            failed_identifier_counter = identifier_list[1]
            if failed_identifier_counter != 0:
                print(
                    str(single_dict)
                    + " failed to get "
                    + str(failed_identifier_counter)
                    + " identifier(s)!"
                )
        print("Updated dictionary with '" + self.identifier_keyname + "'")
        return

    # Section: Check for duplicates

    def get_duplicate_key(self, all_dicts: dict) -> None:
        """
        This function updates the subdictionaries in the given dictionary with the number of duplicates in the identifier list as a new key-value-Pair (key: self.duplicates_keyname) and a list of the duplicated identifier (key: self.duplicates_id_keyname).

        Args:
            all_dicts (dict): Dictionary with subdictionaries including a list of identifiers (self.identifier_keyname).
        """
        for single_dict in all_dicts:
            number_of_duplicates = len(
                all_dicts[single_dict][self.identifier_keyname]
            ) - len(set(all_dicts[single_dict][self.identifier_keyname]))
            all_dicts[single_dict][self.duplicates_keyname] = number_of_duplicates
            duplicates = []
            for mol in all_dicts[single_dict][self.identifier_keyname]:
                if all_dicts[single_dict][self.identifier_keyname].count(mol) > 1:
                    duplicates.append(mol)
            all_dicts[single_dict][self.duplicates_id_keyname] = set(duplicates)
            print(
                "Number of duplicates in "
                + single_dict
                + ": "
                + str(all_dicts[single_dict][self.duplicates_keyname])
                + ",  duplicates: "
                + str(all_dicts[single_dict][self.duplicates_id_keyname])
            )
        print(
            "Updated dictionary with '"
            + self.duplicates_keyname
            + "' and '"
            + self.duplicates_id_keyname
            + "'"
        )
        return

    # Section: Dataset comparison and visualization

    def get_shared_molecules_key(self, all_dicts: dict) -> None:
        """
        This function updates the subdictionaries in the given dictionary (created with the import_as_data_dict function) with the number of molecules that can be found in all of the given datasets (key: self.shared_mols_keyname) and an identifier list of these molecules (key: self.shared_mols_id_keyname) as two new key-value pairs (number of compared datasets can be any number).
        The comparison of the molecules is based on the identifiers (string representation), not the SDMolSupplier Object.

        Args:
            all_dicts (dict): Dictionary with subdictionaries including a lists of identifiers (self.identifier_keyname).
        """
        sets = []
        for single_dict in all_dicts:
            single_set = set(all_dicts[single_dict][self.identifier_keyname])
            sets.append(single_set)
        shared_molecules = set.intersection(*sets)
        for single_dict in all_dicts:
            all_dicts[single_dict][self.shared_mols_keyname] = len(shared_molecules)
            all_dicts[single_dict][self.shared_mols_id_keyname] = shared_molecules
        print(
            "Number of molecules that can be found in all datasets: "
            + str(len(shared_molecules))
            + ", identifiers: "
            + str(shared_molecules), '\n'
            "Updated dictionary with '"
            + self.shared_mols_keyname
            + "' and '"
            + self.shared_mols_id_keyname
            + "'"
        )
        return

    def visualize_intersection(self, all_dicts: dict, data_type: str = "png"):
        """
        This function returns a Venn diagram of the intersection between the molecules in the subdictionaries of the given dictionary. Every subdictionary is represented as a circle and the overlaps between the circles indicate the molecules present in more than one subdictionary. (The function only works with two or three subdictionaries.)
        The intersection is based on the identifiers (string representation).
        The diagram is saved in an output folder.

        Args:
            all_dicts (dict): Dictionary of dictionaries with identifier_keyname.
            data_type (str): Data type for the exported image (default: png).
        Returns:
            fig (matplotlib.figure): Venn diagram

        Raises:
            ValueError: If there is only one or more than three sets to be compared an error is raised.
        """
        sets = []
        for single_dict in all_dicts:
            single_set = set(all_dicts[single_dict][self.identifier_keyname])
            sets.append(single_set)
        fig = plt.figure(figsize=(10, 10))
        if len(sets) == 3:
            venn = venn3(sets, set_labels=(all_dicts.keys()), alpha=0.5)
        elif len(sets) == 2:
            venn = venn2(sets, set_labels=(all_dicts.keys()), alpha=0.5)
        else:
            raise ValueError(
                "Visualization only possible for two or three data sets!"
            )
        plt.title("Intersection as Venn diagram", fontsize=20)
        for text in venn.set_labels:
            text.set_fontsize(15)
        for x in range(len(venn.subset_labels)):
            if venn.subset_labels[x] is not None:
                venn.subset_labels[x].set_fontsize(15)
        if not os.path.exists("output"):
            os.makedirs("output")
        plt.savefig(
            "output/intersection.%s" % (data_type),
            bbox_inches="tight",
            transparent=True,
        )
        plt.close(fig)
        return fig

    # Section: Get descriptors and create plots

    def _get_descriptor_list(
        self, moleculeset: Chem.SDMolSupplier, descriptor: callable,
    ) -> List:
        """
        This function returns a list of descriptor values for all molecules in a given SDMolSupplier object and a callable descriptor (e.g Descriptors.MolWt or rdMolDescriptors.CalcExactMolWt). (private method)

        Args:
            moleculeset (rdkit.Chem.SDMolSupplier)
            descriptor (callable): RDKit method that returns a molecular descriptor for a given molecule.

        Returns:
            List[]: List of descriptor values
        """
        descriptor_list = []
        for mol in moleculeset:
            if mol:
                descriptor_list.append(descriptor(mol))
            else:
                descriptor_list.append(None)
        return descriptor_list

    def get_descriptor_list_key(
        self, all_dicts: dict, descriptor: callable, descriptor_list_keyname: str
    ) -> None:
        """
        This function updates the subdictionaries in the given dictionary with a list of descriptor values as a new key-value pair using _get_descriptor_list on the SDMolSupplier Objects in the subdictionaries.

        Args:
            all_dicts (dict): Dictionary of dictionaries with SDMolSupplier Object.
            descriptor (callable): RDKit method that returns a molecular descriptor for a given molecule
            descriptor_list_keyname (str): Key name for the dictionary entry (should match the descriptor)
        """
        for single_dict in all_dicts:
            descriptor_list = self._get_descriptor_list(
                all_dicts[single_dict][self.import_keyname], descriptor
            )
            all_dicts[single_dict][descriptor_list_keyname] = descriptor_list
        print("Updated dictionary with '" + descriptor_list_keyname + "'")
        return

    def _get_discrete_descriptor_counts(
        self, all_dicts: dict, descriptor_list_keyname: str
    ) -> None:
        """
        This function updates the subdictionaries in the given dictionary with the binned descriptor values for a given descriptor value list with discrete values (e.g. number of H-Bond donors or acceptors). (private method)

        Args:
            all_dicts (dict): Dictionary with subdictionaries including a discrete descriptor value list.
            descriptor_list_keyname (str): Name of the descriptor list.
        """
        binned_descriptor_list_keyname = str("binned_" + descriptor_list_keyname)
        find_max = []
        for single_dict in all_dicts:
            if None in all_dicts[single_dict][descriptor_list_keyname]:
                find_max.append(max([descriptor_value for descriptor_value in all_dicts[single_dict][descriptor_list_keyname] if descriptor_value is not None]))
            else:
                find_max.append(max(all_dicts[single_dict][descriptor_list_keyname]))
        maximum = max(find_max) + 1
        bins = pd.interval_range(start=0, end=maximum, freq=1, closed="left")
        for single_dict in all_dicts:
            counts = pd.value_counts(
                pd.cut(all_dicts[single_dict][descriptor_list_keyname], bins),
                sort=False,
            )
            all_dicts[single_dict][binned_descriptor_list_keyname] = counts
        return

    def _get_continuous_descriptor_counts(
        self, all_dicts: dict, descriptor_list_keyname: str, width_of_bins: float = 10.0
    ) -> None:
        """
        This function updates the subdictionaries in the given dictionary with the binned descriptor values for a given descriptor value list with continuous values (e.g. molecular weight or logP values). The interval size of the bins can be chosen. (private method)

        Args:
            all_dicts (dict): Dictionary with subdictionaries including a continuous descriptor value list.
            descriptor_list_keyname (str): name of the descriptor list.
            width_of_bins (int, optional): Interval size for the bins (default: 10)
        """
        binned_descriptor_list_keyname = str("binned_" + descriptor_list_keyname)
        find_min = []
        find_max = []
        for single_dict in all_dicts:
            find_min.append(min([descriptor_value for descriptor_value in all_dicts[single_dict][descriptor_list_keyname] if descriptor_value is not None]))
            find_max.append(max([descriptor_value for descriptor_value in all_dicts[single_dict][descriptor_list_keyname] if descriptor_value is not None]))
        if min(find_min) < round(min(find_min), ndigits=-1):
            lower = round(min(find_min), ndigits=-1) - 10
        else:
            lower = round(min(find_min), ndigits=-1)
        if max(find_max) > round(max(find_max), ndigits=-1):
            upper = round(max(find_max), ndigits=-1) + 10
        else:
            upper = round(max(find_max), ndigits=-1)
        if (upper - lower) % width_of_bins == 0:
            bins = pd.interval_range(start=lower, end=upper, freq=width_of_bins)
        else:
            to_add = width_of_bins - (upper - lower) % width_of_bins
            bins = pd.interval_range(
                start=lower, end=(upper + to_add), freq=width_of_bins
            )
        for single_dict in all_dicts:
            counts = pd.value_counts(
                pd.cut(all_dicts[single_dict][descriptor_list_keyname], bins),
                sort=False,
            )
            all_dicts[single_dict][binned_descriptor_list_keyname] = counts
        return

    def _discrete_descriptor_plot(
        self,
        all_dicts: dict,
        descriptor_list_keyname: str,
        data_type: str = "png",
        save_dataframe: bool = True,
        figsize: Tuple[float, float] = [15.0, 7.0],
        fontsize_tick_labels: int = 15,
        fontsize_legend: int = 15,
        fontsize_ylabel: int = 20,
        fontsize_xlabel: int = 20,
        fontsize_title: int = 24,
    ):
        """
        This function returns a bar-plot for a discrete descriptor with was previously binned.
        The plot is saved in an output folder as an image (data type can be chosen) and the data frame can also be saved as CSV file.

        args:
            all_dicts (dict): Dictionary with subdictionaries including a binned discrete descriptor.
            descriptor_list_keyname (str): Name of descriptor list for plotting.
            data_type (str): Data type for the exported image (default: png).
            save_dataframe (bool): Export dataframe as csv file or not (default: True).
            figsize (float, float): Width, height of the image in inches (default: 15, 7).
            fontsize_tick_labels (int): Fontsize of the labels on the ticks of the axis (default: 15).
            fontsize_legend (int): Fontsize of the legend (default: 15).
            fontsize_ylabel (int): Fontsize of the label of the y-axis (default: 20).
            fontsize_xlabel (int): Fontsize of the label of the x-axis (default: 20).
            fontsize_title (int): Fontsize of the title (default: 20).

        returns:
            fig (matplotlib.figure): Plot
        """
        binned_descriptor_list_keyname = str("binned_" + descriptor_list_keyname)
        first_dict = list(all_dicts.keys())[0]
        y_max = len(all_dicts[first_dict][binned_descriptor_list_keyname])
        descriptor_df_dict = {
            str("Number of " + descriptor_list_keyname): list(range(y_max))
        }
        for single_dict in all_dicts:
            header = single_dict
            descriptor_df_dict.update(
                {header: list(all_dicts[single_dict][binned_descriptor_list_keyname])}
            )
        descriptor_df = pd.DataFrame(descriptor_df_dict)
        if not os.path.exists("output"):
            os.makedirs("output")
        if save_dataframe:
            descriptor_df.to_csv("output/table_%s.csv" % (descriptor_list_keyname))
        descriptor_plot = descriptor_df.plot(
            x=str("Number of " + descriptor_list_keyname),
            kind="bar",
            stacked=False,
            rot=0,
            figsize=figsize,
            fontsize=fontsize_tick_labels,
        )
        descriptor_plot.legend(bbox_to_anchor=(1, 1), loc="upper left", fontsize=fontsize_legend)
        descriptor_plot.set_ylabel("Number of molecules", fontsize=fontsize_ylabel)
        descriptor_plot.set_xlabel(str(descriptor_list_keyname), fontsize=fontsize_xlabel)
        descriptor_plot.set_title(
            str("Distribution of " + descriptor_list_keyname), pad=20, fontsize=fontsize_title
        )
        fig = descriptor_plot.figure
        fig.savefig(
            "output/distribution_of_%s.%s" % (descriptor_list_keyname, data_type),
            bbox_inches="tight",
            transparent=True,
        )
        plt.close(fig)
        return fig

    def _continuous_descriptor_plot(
        self,
        all_dicts: dict,
        descriptor_list_keyname: str,
        data_type: str = "png",
        save_dataframe: bool = True,
        figsize: Tuple[float, float] = [15.0, 7.0],
        fontsize_tick_labels: int = 15,
        fontsize_legend: int = 15,
        fontsize_ylabel: int = 20,
        fontsize_xlabel: int = 20,
        fontsize_title: int = 24,
    ):
        """
        This function returns bar-plot for a continuous descriptor which was previously binned.
        The plot is saved in an output folder as an image (data type can be chosen) and the data frame can also be saved as CSV file.

        args:
            all_dicts (dict): Dictionary with subdictionaries including a binned continuous descriptor.
            descriptor_list_keyname (str): Name of descriptor list for plotting.
            data_type (str): Data type for the exported image (default: png).
            save_dataframe (bool): Export dataframe as csv file or not (default: True).
            fig_size (float, float): Width, height of the image in inches (default: 15, 7).
            fontsize_tick_labels (int): Fontsize of the labels on the ticks of the axis (default: 15).
            fontsize_legend (int): Fontsize of the legend (default: 15).
            fontsize_ylabel (int): Fontsize of the label of the y-axis (default: 20).
            fontsize_xlabel (int): Fontsize of the label of the x-axis (default: 20).
            fontsize_title (int): Fontsize of the title (default: 20).

        returns:
            fig (matplotlib.figure): Plot
        """
        binned_descriptor_list_keyname = str("binned_" + descriptor_list_keyname)
        first_dict = list(all_dicts.keys())[0]
        descriptor_df_dict = {
            str(descriptor_list_keyname + " Intervals"): all_dicts[first_dict][
                binned_descriptor_list_keyname
            ].keys()
        }
        for single_dict in all_dicts:
            header = single_dict
            descriptor_df_dict.update(
                {header: list(all_dicts[single_dict][binned_descriptor_list_keyname])}
            )
        descriptor_df = pd.DataFrame(descriptor_df_dict)
        if not os.path.exists("output"):
            os.makedirs("output")
        if save_dataframe:
            descriptor_df.to_csv("output/table_%s.csv" % (descriptor_list_keyname))
        descriptor_plot = descriptor_df.plot(
            x=str(descriptor_list_keyname + " Intervals"),
            kind="bar",
            stacked=False,
            figsize=figsize,
            fontsize=fontsize_tick_labels,
        )
        descriptor_plot.legend(bbox_to_anchor=(1, 1), loc="upper left", fontsize=fontsize_legend)
        descriptor_plot.set_ylabel("Number of molecules", fontsize=fontsize_ylabel)
        descriptor_plot.set_xlabel(
            str(descriptor_list_keyname + " Intervals"), fontsize=fontsize_xlabel
        )
        descriptor_plot.set_title(
            str("Distribution of " + descriptor_list_keyname), pad=20, fontsize=fontsize_title
        )
        fig = descriptor_plot.figure
        fig.savefig(
            "output/distribution_of_%s.%s" % (descriptor_list_keyname, data_type),
            bbox_inches="tight",
            transparent=True,
        )
        plt.close(fig)
        return fig

    def descriptor_counts_and_plot(
        self,
        all_dicts: dict,
        descriptor_list_keyname: str,
        width_of_bins: float = 10.0,
        data_type: str = "png",
        save_dataframe: bool = True,
        figsize: Tuple[float, float] = [15.0, 7.0],
        fontsize_tick_labels: int = 15,
        fontsize_legend: int = 15,
        fontsize_ylabel: int = 20,
        fontsize_xlabel: int = 20,
        fontsize_title: int = 24,
    ):
        """
        This function updates the subdictionaries in the given dictionary with the binned descriptor values for a given descriptor value list. The values can either be continuous (binning with _get_continuous_descriptor_counts and plotted with _continuous_descriptor_plot) or discrete (binning with _get_discrete_descriptor_counts and plotted with _discrete_descriptor_plot).
        The created plots are saved in an output folder and the data frame can also be exported as CSV.

        Args:
            all_dicts (dict): Dictionary with subdictionaries including a descriptor value list.
            descriptor_list_keyname (str): Name of the descriptor list for binning and plotting.
            width_of_bins (int, optional): interval size for the bins for continuous values (default: 10).
            data_type (str): Data type for the exported image (default: png).
            save_dataframe (bool): Export dataframe as csv file or not (default: True).
            fig_size (float, float): Width, height of the image in inches (default: 15, 7).
            fontsize_tick_labels (int): Fontsize of the labels on the ticks of the axis (default: 15).
            fontsize_legend (int): Fontsize of the legend (default: 15).
            fontsize_ylabel (int): Fontsize of the label of the y-axis (default: 20).
            fontsize_xlabel (int): Fontsize of the label of the x-axis (default: 20).
            fontsize_title (int): Fontsize of the title (default: 24).

        Returns:
            fig (matplotlib.figure): Plot
        """
        first_dict = list(all_dicts.keys())[0]
        if not (
            any(
                key == descriptor_list_keyname
                for key in list(all_dicts[first_dict].keys())
            )
        ):
            raise KeyError(
                "Descriptor ("
                + str(descriptor_list_keyname)
                + ") needs to be calculated before plotting!"
            )
        elif type(all_dicts[first_dict][descriptor_list_keyname][0]) == int:
            self._get_discrete_descriptor_counts(all_dicts, descriptor_list_keyname)
            fig = self._discrete_descriptor_plot(
                all_dicts, descriptor_list_keyname, data_type, save_dataframe, figsize, fontsize_tick_labels, fontsize_legend, fontsize_ylabel, fontsize_xlabel, fontsize_title
            )
        elif (
            type(all_dicts[first_dict][descriptor_list_keyname][0]) == float
            or type(all_dicts[first_dict][descriptor_list_keyname][0]) == np.float64
        ):
            self._get_continuous_descriptor_counts(
                all_dicts, descriptor_list_keyname, width_of_bins
            )
            fig = self._continuous_descriptor_plot(
                all_dicts, descriptor_list_keyname, data_type, save_dataframe, figsize, fontsize_tick_labels, fontsize_legend, fontsize_ylabel, fontsize_xlabel, fontsize_title
            )
        else:
            raise ValueError(
                'Descriptor values should be "int" or "float" (numpy.float64) to be binned!'
            )
        return fig

    # Section: Check Lipinski Rule of 5 and visualization

    def _test_for_lipinski(self, moleculeset: Chem.SDMolSupplier) -> List[int]:
        """
        This function returns a list with the number of Lipinski Rules broken for every molecule in the given molecule set.

        Args:
            moleculeset (Chem.SDMolSupplier): SDMolSupplier Objects

        Returns:
            list[int]: List of a number of broken Lipinski Rules.
        """
        num_of_break = []
        for mol in moleculeset:
            rule_break = 0
            if not mol:
                num_of_break.append(None)
                continue
            if Descriptors.MolLogP(mol) > 5:
                rule_break += 1
            if Descriptors.MolWt(mol) > 500:
                rule_break += 1
            if Descriptors.NumHAcceptors(mol) > 10:
                rule_break += 1
            if Descriptors.NumHDonors(mol) > 5:
                rule_break += 1
            num_of_break.append(rule_break)
        return num_of_break

    def get_lipinski_key(self, all_dicts: dict) -> None:
        """
        This function updates the subdictionaries in the given dictionary with the list of the number of broken Lipinski Rules for every molecule (lipinski_list_keyname) and a summary of the broken rules (lipinski_summary_keyname) using _test_for_lipinski.

        Args:
            all_dicts (dict): Dictionary with subdictionaries including the key 'self.import_keyname'.
        """
        for single_dict in all_dicts:
            lipinski_break_list = self._test_for_lipinski(
                all_dicts[single_dict][self.import_keyname]
            )
            all_dicts[single_dict][self.lipinski_list_keyname] = lipinski_break_list
            lipinski_summary = {
                "lipinski_molecules": lipinski_break_list.count(0),
                "1_rule_broken": lipinski_break_list.count(1),
                "2_rules_broken": lipinski_break_list.count(2),
                "3_rules_broken": lipinski_break_list.count(3),
                "4_rules_broken": lipinski_break_list.count(4),
            }
            all_dicts[single_dict][self.lipinski_summary_keyname] = lipinski_summary
        print(
            "Updated dictionary with '"
            + self.lipinski_summary_keyname
            + "' and '"
            + self.lipinski_list_keyname
            + "'"
        )
        return

    def lipinski_plot(
        self,
        all_dicts: dict,
        data_type: str = "png",
        save_dataframe: bool = True,
        figsize: Tuple[float, float] = [15.0, 7.0],
        fontsize_tick_labels: int = 15,
        fontsize_legend: int = 15,
        fontsize_ylabel: int = 20,
        fontsize_xlabel: int = 20,
        fontsize_title: int = 24,
    ):
        """
        This function returns a bar plot for the number of molecules in every subdictionary breaking 0 to 4 Lipinski rules using the 'lipinski_summary' key in the given dictionary. The plot is saved in an output folder (data type can be chosen) and the created data frame can also be exported as CSV.

        args:
            all_dicts (dict): Dictionary with subdictionaries including the key 'self.lipinski_summary_keyname'.
            data_type (str): Data type for the exported image (default: png).
            save_dataframe (bool): Export dataframe as csv or not (default: True).
            fig_size (float, float): Width, height of the image in inches (default: 15, 7).
            fontsize_tick_labels (int): Fontsize of the labels on the ticks of the axis (default: 15).
            fontsize_legend (int): Fontsize of the legend (default: 15).
            fontsize_ylabel (int): Fontsize of the label of the y-axis (default: 20).
            fontsize_xlabel (int): Fontsize of the label of the x-axis (default: 20).
            fontsize_title (int): Fontsize of the title (default: 24).

        returns:
            fig (matplotlib.figure): Plot
        """
        lipinski_df_dict = {"Number of broken rules": [0, 1, 2, 3, 4]}
        for single_dict in all_dicts:
            if not (
                any(
                    key == self.lipinski_summary_keyname
                    for key in list(all_dicts[single_dict].keys())
                )
            ):
                raise KeyError(
                    "Lipinski summary ("
                    + str(self.lipinski_summary_keyname)
                    + ") needs to be calculated before plotting! Use 'get_lipinski_key'!"
                )
            header = single_dict
            lipinski_df_dict.update(
                {
                    header: list(
                        all_dicts[single_dict][self.lipinski_summary_keyname].values()
                    )
                }
            )
        lipinski_df = pd.DataFrame(lipinski_df_dict)
        if not os.path.exists("output"):
            os.makedirs("output")
        if save_dataframe:
            lipinski_df.to_csv("output/table_lipinski_rules.csv")
        lipinski_plot = lipinski_df.plot(
            x="Number of broken rules",
            kind="bar",
            stacked=False,
            rot=0,
            figsize=figsize,
            fontsize=fontsize_tick_labels,
        )
        lipinski_plot.legend(bbox_to_anchor=(1, 1), loc="upper left", fontsize=fontsize_legend)
        lipinski_plot.set_ylabel("Number of molecules", fontsize=fontsize_ylabel)
        lipinski_plot.set_xlabel("Number of broken rules", fontsize=fontsize_xlabel)
        lipinski_plot.set_title(
            "Distribution of the number of broken Lipinski Rules", pad=20, fontsize=fontsize_title
        )
        fig = lipinski_plot.figure
        fig.savefig(
            "output/lipinski_rules_plot.%s" % (data_type),
            bbox_inches="tight",
            transparent=True,
        )
        plt.close(fig)
        return lipinski_plot.figure

    # Section: Scaffold analysis and plotting

    def _get_Murcko_scaffold(
        self,
        moleculeset: Chem.SDMolSupplier,
        number_of_scaffolds: int = 5,
        scaffolds_per_row: int = 5,
        image_size: int = 200
    ) -> Tuple[list, count]:  # what to use here?
        """
        This function creates a grid images of a chosen number of Murcko scaffolds for the molecules in a given SDMolSupplier Object. The scaffolds are sorted by their frequency.
        The relative number of occurrence of a scaffold in the dataset is shown below each image.

        args:
            moleculeset (Chem.SDMolSupplier): SDMolSupplier Objects
            number_of_scaffolds (int): Number of scaffolds displayed in the grid images (default: 5).
            scaffolds_per_row (int): Number of scaffolds in every row of the grid image (default: 5).
            image_size (int): Size of the image for a single molecule (default: 200).

        returns:
            scaffold_grid (PIL.PngImageFile): Grid image with most frequent scaffolds.
        """
        scaffold_list = []
        for mol in moleculeset:
            scaffold = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
            scaffold_list.append(scaffold)
        scaffold_counts = pd.Index(scaffold_list).value_counts(normalize=True)
        if len(scaffold_counts) < number_of_scaffolds:
            number_of_scaffolds = len(scaffold_counts)
        legend = [str(integer) for integer in (list(scaffold_counts)[: number_of_scaffolds])]
        smiles_list = list(scaffold_counts.keys())
        to_draw = []
        for mol in range(number_of_scaffolds):
            to_draw.append(Chem.MolFromSmiles(smiles_list[mol]))
        scaffold_grid = Draw.MolsToGridImage(
            to_draw,
            molsPerRow=scaffolds_per_row,
            subImgSize=(image_size, image_size),
            legends=legend,
            returnPNG=False,
        )
        return scaffold_grid, scaffold_list, scaffold_counts

    def draw_scaffolds(
        self,
        all_dicts: dict,
        number_of_scaffolds: int = 5,
        scaffolds_per_row : int = 5,
        image_size: int = 200,
        data_type: str = 'png',
        figsize: Tuple[float, float] = [20.0, 20.0],
        fontsize_title: int = 24,
        fontsize_subtitle: int = 20
    ):
        """
        This function creates a grid images of a chosen number of Murcko scaffolds for every subdictionary in the given dictionary and shows them together. The scaffolds in each gird image are sorted by their frequency. The relative number of occurrence of a scaffold in the dataset of the respective subdictionary is shown below each scaffold image.

        args:
            all_dicts (dict): Dictionary with subdictionaries including the key 'self.import_keyname'.
            number_of_scaffolds (int): Number of scaffolds displayed in the grid images (default: 5).
            scaffolds_per_row (int): Number of scaffolds in every row of the grid image (default: 5).
            image_size (int): Size of the image for a single molecule (default: 200).
            data_type (str): Data type for the exported image (default: png).
            figsize (float, float): Width, height of the image in inches (default: 20, 20)
            fontsize_title (int): Fontsize of the title (default: 24).
            fontsize_subtitle (int): Fontsize of the subtitles (default: 20).

        returns:
            scaffold_grid (PIL.PngImageFile): Grid image with most frequent scaffolds.
        """
        image_list = []
        title_list = []
        for single_dict in all_dicts:
            title_list.append(single_dict)
            scaffolds = self._get_Murcko_scaffold(
                all_dicts[single_dict][self.import_keyname],
                number_of_scaffolds,
                scaffolds_per_row,
                image_size,
            )
            image_list.append(scaffolds[0])
            all_dicts[single_dict][self.scaffold_list_keyname] = scaffolds[1]
            all_dicts[single_dict][self.scaffold_summary_keyname] = (scaffolds[2].to_frame('frequency')).rename_axis('scaffold SMILES')
        print(
            "Updated dictionary with '"
            + self.scaffold_list_keyname
            + "' and '"
            + self.scaffold_summary_keyname
            + "'"
        )
        rows = len(image_list)
        fig = plt.figure(figsize=figsize)
        for j in range(0, rows):
            fig.add_subplot(rows, 1, j + 1)
            plt.axis("off")
            plt.imshow(image_list[j])
            plt.title(title_list[j], fontsize=fontsize_subtitle)
        fig.suptitle("Most frequent Murcko scaffolds from the datasets", fontsize=fontsize_title)
        if not os.path.exists("output"):
            os.makedirs("output")
        fig.savefig("output/scaffold_grid.%s" % (data_type))
        plt.close(fig)
        return fig

    # Section: Chemical space visualization

    def chemical_space_visualization(
        self,
        all_dicts: dict,
        fp_radius: int = 2,
        fp_bits: int = 2048,
        dimension_reduction: str = "pca",
        interactive: bool = True,
    ):
        """
        This function returns a 2D visualization of the chemical space of the molecules in all datasets using the chemplot module.
        On basis of the calculated identifier (self.identifier_keyname) for every molecule a Extended Connectivity Fingerprint (ECFP) will be calculated with a definable fingerprint radius (fp_radius) and length (fp_size).
        Subsequent, the fingerprints are reduced to 2D for plotting. The dimension reduction can be done with PCA, UMAP or t-SNE and the plot can be interactive.

        Args:
            all_dicts (dict): Dictionary with subdictionaries including an identifier list (self.identifier_keyname).
            fp_radius (int): Radius of the Extended Connectivity Fingerprints (default: 2).
            fp_bits (int): Size of the Extended Connectivity Fingerprints (default: 2048).
            dimension_reduction (str): Method of dimension reduction (default: pca).
            interactive (bool): Creating an interactive plot or not (default: True).

        Returns:
            Chemical space visualization
        """
        all_mols_list = []
        target_list = []
        for single_dict in all_dicts:
            for mol in all_dicts[single_dict][self.identifier_keyname]:
                all_mols_list.append(mol)
                target_list.append(single_dict)
        if all_mols_list[0].startswith("InChI="):  # check if identifier is InChI
            chem_space = cp.Plotter.from_inchi(
                all_mols_list,  # list of inchi strings which are used to get Extended Connectivity Fingerprint (alternative: smiles),
                target=target_list,  # corresponding list for inchi_list, shows which dataset the molecules belong to
                target_type="C",  # classification (classes are the datasets listed in the target_list)
                sim_type="structural",  # similarity solely based on structure (no property is taken into account)
            )
        elif (
            len(all_mols_list[0]) == 27  # check if identifier is InChIKey
            and "-" in all_mols_list[0][14]
            and "-" in all_mols_list[0][25]
        ):
            all_mols_list.clear()
            for single_dict in all_dicts:
                for mol in all_dicts[single_dict][self.import_keyname]:
                    inchi = Chem.MolToInchi(mol)
                    all_mols_list.append(inchi)
            chem_space = cp.Plotter.from_inchi(
                all_mols_list,  # list of inchi strings which are used to get Extended Connectivity Fingerprint (alternative: smiles),
                target=target_list,  # corresponding list for inchi_list, shows which dataset the molecules belong to
                target_type="C",  # classification (classes are the datasets listed in the target_list)
                sim_type="structural",  # similarity solely based on structure (no property is taken into account)
            )
        else:
            chem_space = cp.Plotter.from_smiles(
                all_mols_list,  # list of smiles strings which are used to get Extended Connectivity Fingerprint (alternative: smiles),
                target=target_list,  # corresponding list for inchi_list, shows which dataset the molecules belong to
                target_type="C",  # classification (classes are the datasets listed in the target_list)
                sim_type="structural",  # similarity solely based on structure (no property is taken into account)
            )
        if fp_radius != 2 or fp_bits != 2048:
            if all_mols_list[0].startswith("InChI="):
                new_fingerprint = descriptors.get_ecfp_from_inchi(
                    all_mols_list, target_list, radius=fp_radius, nBits=fp_bits
                )
            else:
                new_fingerprint = descriptors.get_ecfp(
                    all_mols_list, target_list, radius=fp_radius, nBits=fp_bits
                )
            chem_space._Plotter__mols = new_fingerprint[0]
        if dimension_reduction == "pca":
            chem_space.pca()  # n_components, copy, whiten, svd_solver ...
        elif dimension_reduction == "tsne":
            chem_space.tsne()  # n_components, perplexity, learning_rate, n_iter, init, random_state ...
        elif dimension_reduction == "umap":
            chem_space.umap()  # n_neighbors, min_dist, pca, random_state ...
        else:
            raise ValueError('dimension_reduction should be "pca", "tsne" or "umap"!')
        if not os.path.exists("output"):
            os.makedirs("output")
        if not interactive:
            chem_space.visualize_plot(filename="output/chemical_space")
        else:
            chem_space.interactive_plot(show_plot=True)
        return chem_space

    # Section: Data export

    def export_single_dict_values(self, all_dicts: dict) -> None:
        """
        This function exports the descriptor values for each dictionary according to one imported SDFile as a single csv file in the output folder.

        Args:
            all_dicts (dict): Dictionary of dictionaries with calculated descriptor values.
        """
        for single_dict in all_dicts:
            new_dict = all_dicts[single_dict].copy()
            counter = 0
            for key in new_dict.copy():
                if type(new_dict[key]) == list:  # and len(new_dict[key]) == 100:
                    counter += 1
                else:
                    new_dict.pop(key)
            to_export = pd.DataFrame(new_dict)
            filename = single_dict[:-4]
            to_export.to_csv("output/descriptor_values_%s.csv" % (filename))
            print(single_dict + " : " + str(counter) + " exported descriptor values")
        return

    def export_all_figures_pdf(self) -> None:
        """
        This function returns a pdf including all created images from the output folder.

        Returns:
            pdf containing all images created.
        """
        pdf = FPDF()
        for image in os.listdir("output"):
            if image[-3:] == "png" or image[-3:] == "jpg":
                pdf.add_page()
                pdf.image("output/%s" % (image), 5, 5, 200)
            elif image[-3:] == "csv":
                pass
            else:
                print(image + " not included, due to unsupported image type.")
        pdf.output("output/all_figures.pdf")
        return
