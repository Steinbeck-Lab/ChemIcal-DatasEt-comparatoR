"""
MIT License

Copyright (c) 2022 Hannah Busch, Jonas Schaub, Otto Brinkhaus, Kohulan Rajan,
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
SOFTWARE.
"""

# Section: Import Libraries

import io
import logging
import os
import sys
import pandas as pd
from pandas.errors import ParserError
import numpy as np

from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    AtomValenceException,
    Descriptors,
    Draw,
    Lipinski,
    KekulizeException,
    rdchem,
)
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Scaffolds import MurckoScaffold

import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import matplotlib

import chemplot as cp
from bokeh.plotting import output_file, save

from typing import List, Tuple, Dict
from itertools import count
import PIL
from PIL import Image

from fpdf import FPDF, YPos, XPos
from fpdf.enums import Align
from datetime import date, datetime

logger = logging.getLogger("CIDER")  # needs to be here

# Section: Constructor


class ChemicalDatasetComparator:
    """
    Wrapper class around all Cider functionalities

    ChemIcal DatasEt comparatoR (CIDER) is a Python package and ready-to-use
    Jupyter Notebook workflow which primarily utilizes RDKit to compare two or
    more chemical structure datasets (SD files) in different aspects, e.g.
    size, overlap, molecular descriptor distributions, chemical space
    clustering, etc., most of which can be visually inspected.
    """

    def __init__(self):
        """
        The class variables of CIDER function as keys for the dictionary in
        which all calculated and plotted data will be stored. The keys will be
        generated when the CIDER method is executed.
        """
        # from cider.import_as_data_dict
        self.import_keyname = "rdkit_mol_Object"
        self.figure_dict_keyname = "figures"
        # from cider.get_number_of_molecules
        self.dataset_length_keyname = "number_of_molecules"
        # from cider.get_identifier_list_key
        self.identifier_keyname = "identifier_list"
        # cider.get_duplicate_key
        self.duplicates_keyname = "number_of_duplicates"
        self.duplicates_id_keyname = "duplicates"
        self.duplicates_index_keyname = "duplicates_index"
        # from cider.get_shared_molecules_key
        self.shared_mols_keyname = "number_of_shared_molecules"
        self.shared_mols_id_keyname = "shared_molecules"
        # from cider.get_lipinski_key
        self.lipinski_list_keyname = "number_of_broken_lipinski_rules"
        self.lipinski_summary_keyname = "lipinski_summary"
        # from cider.get_database_id
        self.database_id_keyname = "coconut_id"
        # from cider.draw_most_frequent_scaffolds
        self.scaffold_list_keyname = "scaffold_list"
        self.scaffold_summary_keyname = "scaffold_summary"

    # Section: Configuration for logging

    if not os.path.exists("output"):
        os.mkdir("output")
    if not os.path.exists("output/logs"):
        os.mkdir("output/logs")

    now = (str(datetime.now())[:-7]).replace(":", "-")

    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        level=logging.INFO,
        handlers=[
            logging.FileHandler("output/logs/%s_cider_logging.log" % (now)),
            logging.StreamHandler(sys.stdout),
        ],
    )
    logger = logging.getLogger("CIDER")

    def exception_handler_IP(self, exc_type, exc_value, exc_tb, tb_offset=None):
        """
        This function is a callable exception handler when working in IPython.
        It logs the exception including the traceback into the log file.

        Args:
            exc_type: Type of the exception.
            exc_value: Value of the exception.
            exc_tb: Traceback of the exception.
            tb_offset: Traceback offset.
        """
        logger.error("An Error occured while executing CIDER!", exc_info=True)
        self.showtraceback((exc_type, exc_value, exc_tb), tb_offset=tb_offset)
        return

    def exception_handler(exc_type, exc_value, exc_tb):
        """
        This function is a callable exception handler when NOT working in
        IPython. It logs the exception including the traceback into the log
        file.

        Args:
            exc_type: Type of the exception.
            exc_value: Value of the exception.
            exc_tb: Tracback of the exception.
            tb_offset: Traceback offset.
        """
        logger.exception(
            "An Error occured while executing CIDER!",
            exc_info=(exc_type, exc_value, exc_tb),
        )
        return

    """
    This function checks if CIDER is run in IPython or not and selects the
    corresponding exception handler for logging.
    """
    try:
        __IPYTHON__
        from IPython import get_ipython

        ip = get_ipython()
        ip.set_custom_exc((BaseException,), exception_handler_IP)
    except NameError:
        sys.excepthook = exception_handler

    # Section: Import data and check for faulty SDFiles

    def _check_invalid_mols_in_SDF(self, all_dicts: dict) -> None:
        """
        This function checks if there are invalid entries in the SDFiles/SMI
        files that can cause errors in the subsequent functions. Those invalid
        entries will be removed form the rdkit_mol_objects. The entry will
        remain in the original SDFile/SMI File as it is. Header lines from the
        SMI File will also be removed. (private method)

        Args:
            all_dicts (dict): Dictionary with sub-dictionaries including
            rdkit_mol_objects.
        """
        for single_dict in all_dicts:
            if single_dict == self.figure_dict_keyname:
                continue
            mol_index = -1
            invalid_index = []
            for mol in all_dicts[single_dict][self.import_keyname]:
                mol_index += 1
                if not mol:
                    logger.warning(
                        "%s has invalid molecule at index %d" % (single_dict, mol_index)
                    )
                    invalid_index.append(mol_index)
            if not invalid_index:
                logger.info("No faulty molecules found in %s" % (single_dict))
            else:
                new_SDMol = list(all_dicts[single_dict][self.import_keyname])
                for index in sorted(invalid_index, reverse=True):
                    del new_SDMol[index]
                all_dicts[single_dict].update({self.import_keyname: new_SDMol})
                if self.database_id_keyname in all_dicts[single_dict]:
                    new_id_list = list(all_dicts[single_dict][self.database_id_keyname])
                    for index in sorted(invalid_index, reverse=True):
                        del new_id_list[index]
                    all_dicts[single_dict].update(
                        {self.database_id_keyname: new_id_list}
                    )
                logger.info(
                    "%d invalid molecule(s) deleted from %s"
                    % (len(invalid_index), single_dict)
                )
        return

    def import_as_data_dict(self, path_to_data: str) -> Dict:
        """
        This function creates a dictionary with the names of the imported file
        as keys. The values of each of these keys is a subdictionary. The first
        entry of every subdictionary is self.import_keyname (class variable,
        can be changed) as key and rdkit_mol_objects (either a rdkit.Chem.
        rdmolfiles.SDMolSupplier Object or a list of rdkit.Chem.rdchem.mol
        Objects) of the SDFile as value. To find faulty molecules every entry
        of the rdkit_mol_objects will be parsed once. (Parsed molecules will
        not be stored in the dictionary to save memory.)

        Args:
            path_to_data (str): Path to the directory where the SDFiles are
            stored.

        Returns:
            all_dicts (dict): Dictionary with subdictionaries for every
            dataset, including the key 'self.import_keyname'.

        Raises:
            FileNotFoundError: if the data path is invalid.
            KeyError: if no SDFiles are in the given directory.
        """
        all_dicts = {}
        data_dir = os.path.abspath(str(path_to_data))
        for dict_name in os.listdir(data_dir):
            if dict_name[-3:] == "sdf" or dict_name[-3:] == "SDF":
                single_dict = {}
                dict_path = os.path.join(data_dir, dict_name)
                single_dict[self.import_keyname] = Chem.SDMolSupplier(dict_path)
                all_dicts[dict_name] = single_dict
        if not all_dicts:
            raise KeyError("No SDFiles found in the given directory %s!" % (data_dir))
        figure_dict = {}
        all_dicts[self.figure_dict_keyname] = figure_dict
        self._check_invalid_mols_in_SDF(all_dicts)
        logger.info("Created dictionary with keys: %s", list(all_dicts.keys()))
        os.chdir(os.path.dirname(data_dir))
        if not os.path.exists("output"):
            os.mkdir("output")
        else:
            if os.listdir("output"):
                logger.warning(
                    "Already existing output folder with files! Old data will be overwritten!"
                )
        return all_dicts

    def import_smi_as_data_dict(self, path_to_data: str):
        """
        This function creates a dictionary with the names of the imported file
        as keys. The values of each of these keys is a subdictionary. The first
        entry of every subdictionary is self.import_keyname (class variable,
        can be changed) as key and rdkit_mol Objects (list of rdkit.Chem.rdchem.
        mol objects) of the SMI File as value. To find faulty molecules, every
        entry of the rdkit_mol Objects will be parsed once. (Parsed molecules
        will not be stored in the dictionary to save memory.)

        Args:
            path_to_data (str): Path to the directory where the SMI Files are
            stored.

        Returns:
            all_dicts (dict): Dictionary with subdictionaries for every
            dataset, including the key 'self.import_keyname'.

        Raises:
            FileNotFoundError: if the data path is invalid.
            KeyError: if no SMI Files are in the given directory.
        """
        all_dicts = {}
        data_dir = os.path.abspath(str(path_to_data))
        for dict_name in os.listdir(data_dir):
            if (
                dict_name[-3:] == "smi"
                or dict_name[-3:] == "SMI"
                or dict_name[-3:] == "txt"
            ):
                single_dict = {}
                dict_path = os.path.join(data_dir, dict_name)
                try:
                    smi_table = pd.read_csv(
                        dict_path, sep=None, engine="python", header=None
                    )
                except ParserError:
                    smi_table = pd.read_csv(dict_path, header=None)
                for column in range(len(smi_table.columns)):
                    is_mol = []
                    for row in range(3):
                        is_mol.append((Chem.MolFromSmiles(smi_table[column][row])))
                    if any(is_mol):
                        smi_column = column
                        break
                rdkit_mol_list = []
                for mol in smi_table[smi_column]:
                    rdkit_mol_list.append(Chem.MolFromSmiles(mol))
                single_dict[self.import_keyname] = rdkit_mol_list
                all_dicts[dict_name] = single_dict
                if id:
                    try:
                        if smi_column == 1:
                            id_column = 0
                        elif smi_column == 0:
                            id_column = 1
                        id_list = list(smi_table[id_column])
                        single_dict[self.database_id_keyname] = id_list
                    except KeyError:
                        logger.info(
                            "Cannot find IDs for file %s! SMILES strings and database ID should be the first a d second entry of the files to import the ID."
                            % (dict_name)
                        )
                        continue
        if not all_dicts:
            raise KeyError("No SMI files found in the given directory %s!" % (data_dir))
        figure_dict = {}
        all_dicts[self.figure_dict_keyname] = figure_dict
        self._check_invalid_mols_in_SDF(all_dicts)
        logger.info("Created dictionary with keys: %s", list(all_dicts.keys()))
        os.chdir(os.path.dirname(data_dir))
        if not os.path.exists("output"):
            os.mkdir("output")
        else:
            if os.listdir("output"):
                logger.warning(
                    "Already existing output folder with files! Old data will be overwritten!"
                )
        return all_dicts

    # Section: Saving figures and images

    def _save_to_figure_dict(
        self, all_dicts: dict, keyname: str, fig, data_type: str = "png"
    ) -> None:
        """
        This function stores the images and figures created by CIDER in the
        'figure' (self.figure_dict_keyname) subdictionary and exports them to
        the output folder with a given data type. When the name for a figure is
        already used there will be a increasing number added. This only works
        within one analysis, after restarting CIDER the files in the output
        folder might be overwritten.

        Args:
            all_dicts (dict): Dictionary containing 'figures' subdictionary
            (self.figure_dict_keyname).
            keyname (str): Name of the image/figure.
            fig: Images or figure to be stored.
            data_type (str): Data type for the exported file. (Default: png)
        """
        if not any(
            key == keyname for key in list(all_dicts[self.figure_dict_keyname].keys())
        ):
            all_dicts[self.figure_dict_keyname][keyname] = fig
            plt.savefig(
                "output/%s.%s" % (keyname, data_type),
                bbox_inches="tight",
                transparent=True,
            )
            logger.info("Updated dictionary with '%s'", keyname)
        else:
            counter = 1
            new_keyname = keyname + "_" + str(counter)
            while any(
                key == new_keyname
                for key in list(all_dicts[self.figure_dict_keyname].keys())
            ):
                counter += 1
                new_keyname = keyname + "_" + str(counter)
            all_dicts[self.figure_dict_keyname][new_keyname] = fig
            plt.savefig(
                "output/%s.%s" % (new_keyname, data_type),
                bbox_inches="tight",
                transparent=True,
            )
            logger.info("Updated dictionary with '%s'", new_keyname)
        return

    # Section: Get overview of the dataset size and molecules

    def get_number_of_molecules(self, all_dicts: dict) -> None:
        """
        This function updates the subdictionaries in the given dictionary
        (created from import_as_data_dict function) with the number of
        molecules in every dataset as new key-value pair. The key is the class
        variable 'cider.dataset_length_keyname'.

        Args:
            all_dicts (dict): Dictionary with subdictionaries for every
            dataset, including the key 'self.import_keyname'.
        """
        for single_dict in all_dicts:
            if single_dict == self.figure_dict_keyname:
                continue
            number_of_molecules = len(all_dicts[single_dict][self.import_keyname])
            all_dicts[single_dict][self.dataset_length_keyname] = number_of_molecules
            logger.info(
                "Number of molecules in %s: %d"
                % (single_dict, all_dicts[single_dict][self.dataset_length_keyname])
            )
        logger.info("Updated dictionary with '%s'", self.dataset_length_keyname)
        return

    def draw_molecules(
        self,
        all_dicts: dict,
        number_of_mols: int = 10,
        mols_per_row: int = 5,
        image_size: int = 200,
        data_type: str = "png",
        figsize: Tuple[float, float] = [20.0, 20.0],
        fontsize_title: int = 24,
        fontsize_subtitle: int = 20,
    ) -> matplotlib.figure.Figure:
        """
        This function creates an grid image of the first molecules of each
        dataset and exports the image to an output folder.

        Args:
            all_dicts (dict): Dictionary with subdictionaries for every
            dataset, including the key 'self.import_keyname'.
            number_of_mols (int): number of molecules form each dataset that
            will be displayed (default: 12).
            mols_per_row (int): number of molecules per row in the grid
            (default: 3).
            image_size (int): the size of the image for a single molecule
            (default: 200).
            data_type (str): data type for the exported files (e.g. png, jpg,
            pdf, default: png).
            figsize (float, float): Width, height of the figure in inches
            (default: 20, 20)
            fontsize_title (int): Fontsize of the title (default: 24).
            fontsize_subtitle (int): Fontsize of the subtitles (default: 20).

        Returns:
            fig (matplotlib.figure): grid image of molecules
        """
        image_list = []
        title_list = []
        for single_dict in all_dicts:
            if single_dict == self.figure_dict_keyname:
                continue
            title_list.append(single_dict)
            to_draw = []
            if len(all_dicts[single_dict][self.import_keyname]) < number_of_mols:
                number_of_mols_final = len(all_dicts[single_dict][self.import_keyname])
            else:
                number_of_mols_final = number_of_mols
            for i in range(number_of_mols_final):
                to_draw.append(all_dicts[single_dict][self.import_keyname][i])
            for mol in to_draw:
                atom0_pos = [
                    mol.GetConformer().GetAtomPosition(0).x,
                    mol.GetConformer().GetAtomPosition(0).y,
                    mol.GetConformer().GetAtomPosition(0).z,
                ]
                atom1_pos = [
                    mol.GetConformer().GetAtomPosition(1).x,
                    mol.GetConformer().GetAtomPosition(1).y,
                    mol.GetConformer().GetAtomPosition(1).z,
                ]
                if atom0_pos == atom1_pos:
                    AllChem.Compute2DCoords(mol)
            mol_grid = Draw.MolsToGridImage(
                to_draw,
                maxMols=number_of_mols_final,
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
        self._save_to_figure_dict(all_dicts, "mol_grid", fig, data_type=data_type)
        plt.close(fig)
        return fig

    # Section: Get database ID

    def get_database_id(self, all_dicts: dict, id_name: str) -> None:
        """
        This function updates subdictionaries of a given dictionary with a list
        of database IDs for the single molecules as new key-value pairs.
        Depending on which database the molecules are coming from, the key as a
        class variable can be changed accordingly.
        (To get the database ID the rdkit_mol Objects (rdkit.Chem.rdmol.Mol or
        rdkit.Chem.rdmolfiles.SDMolSupplier) needs to be parsed, this may take
        same time because no parsed molecules are saved in the dictionary to
        save memory.)

        Args:
            all_dicts (dict): Dictionary with subdictionaries for every
            dataset, including the key 'self.import_keyname'.
            id_name (str): ID name in the original SDFile.
        """
        id_count = 0
        for single_dict in all_dicts:
            if single_dict == self.figure_dict_keyname:
                continue
            database_id_list = []
            for mol in all_dicts[single_dict][self.import_keyname]:
                prop_dict = mol.GetPropsAsDict()
                database_id = prop_dict.get(id_name)
                if database_id is not None:
                    id_count += 1
                database_id_list.append(database_id)
                all_dicts[single_dict][self.database_id_keyname] = database_id_list
        if id_count == 0:
            logger.info(
                "No database IDs with '%s' found. (Maybe check for spelling mistakes)"
                % (id_name)
            )
        logger.info("Updated dictionary with '%s'", self.database_id_keyname)
        return

    # Section: Get string identifier

    def _get_identifier_list(
        self, moleculeset, id_type: str = "inchi"
    ) -> Tuple[list, count]:
        """
        This function returns a list of InChI, InChIKey or canonical SMILES
        strings for all molecules in the given rdkit_mol Objects (rdkit.Chem.
        rdmol.Mol or rdkit.Chem.rdmolfile.SDMolSupplier). (private method)

        Args:
            moleculeset (rdkit.Chem.rdmolfile.SDMolSupplier or list[rdkit.Chem.
            rdmol.Mol]):
            id_type (str, optional): "inchi", "inchikey" or "smiles". Defaults
            to "inchi".

        Raises:
            ValueError: if ID_type is not "inchi," "inchikey" or "smiles".

        Returns:
            list[str]: List of identifiers based on given molecules.
            int: Counter of molecules for which no identifier could be
            determined
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
        This function updates the subdictionaries in the given dictionary
        (created with the import_as_data_dict function) with a list of
        identifiers (InChI, InChIKey, canonical SMILES strings) as a new
        key-value pair using  _get_identifier_list on the rdkit_mol Objects
        (rdkit.Chem.rdmol.Mol or rdkit.Chem.rdmolfiles.SDMolSupplier). The key
        self.identifier_keyname (class variable) can be changed.

        Args:
            all_dicts (dict): Dictionary with subdictionaries for every
            dataset, including the key 'self.import_keyname'.
            id_type (str): Type of Identifier ("inchi", "inchikey" or "smiles")
        """
        for single_dict in all_dicts:
            if single_dict == self.figure_dict_keyname:
                continue
            identifier_list = self._get_identifier_list(
                all_dicts[single_dict][self.import_keyname], id_type
            )
            all_dicts[single_dict][self.identifier_keyname] = identifier_list[0]
            failed_identifier_counter = identifier_list[1]
            if failed_identifier_counter != 0:
                logger.warning(
                    "%s failed to get %d identifier(s)!",
                    single_dict,
                    failed_identifier_counter,
                )
        logger.info("Updated dictionary with '%s'", self.identifier_keyname)
        return

    # Section: Check for duplicates

    def get_duplicate_key(self, all_dicts: dict) -> None:
        """
        This function updates the subdictionaries in the given dictionary with
        the number of duplicates in the identifier list as a new key-value-Pair
        (key: self.duplicates_keyname), a list of the duplicated identifier
        (key: self.duplicates_id_keyname) and a list of the indices of the
        duplicates in the rdkit_mol Object (key self.duplicates_index_keyname).

        Args:
            all_dicts (dict): Dictionary with subdictionaries including a list
            of identifiers (self.identifier_keyname).

        Raises:
            KeyError: if there is no identifier list.
        """
        for single_dict in all_dicts:
            if single_dict == self.figure_dict_keyname:
                continue
            if not any(
                key == self.identifier_keyname
                for key in list(all_dicts[single_dict].keys())
            ):
                raise KeyError(
                    "A identifier list is needed, please run 'get_identifier_list_key'!"
                )
            mol_id_dict = {}
            duplicates = set()
            index = -1
            for mol in all_dicts[single_dict][self.identifier_keyname]:
                index += 1
                if mol not in mol_id_dict.keys():
                    mol_id_dict[mol] = [index]
                else:
                    duplicates.add(mol)
                    mol_id_dict[mol].append(index)
            all_dicts[single_dict][self.duplicates_keyname] = len(duplicates)
            all_dicts[single_dict][self.duplicates_id_keyname] = duplicates
            all_dicts[single_dict][self.duplicates_index_keyname] = []
            for mol in duplicates:
                all_dicts[single_dict][self.duplicates_index_keyname].append(mol_id_dict[mol])
            logger.info(
                "Number of duplicates in %s: %d, duplicate identifier(s): %s, duplicate index: %s",
                single_dict,
                all_dicts[single_dict][self.duplicates_keyname],
                all_dicts[single_dict][self.duplicates_id_keyname],
                all_dicts[single_dict][self.duplicates_index_keyname]
            )
        logger.info(
            "Updated dictionary with '%s', '%s' and '%s'",
            self.duplicates_keyname,
            self.duplicates_id_keyname,
            self.duplicates_index_keyname
        )
        return

    # Section: Dataset comparison and visualization

    def get_shared_molecules_key(self, all_dicts: dict) -> None:
        """
        This function updates the subdictionaries in the given dictionary
        (created with the import_as_data_dict function) with the number of
        molecules that can be found in all of the given datasets (key:
        self.shared_mols_keyname) and an identifier list of these molecules
        (key: self.shared_mols_id_keyname) as two new key-value pairs (number
        of compared datasets can be any number).
        The comparison of the molecules is based on the identifiers (string
        representation), not the rdkit_mol Object (rdkit.Chem.rdmol.Mol or
        rdkit.Chem.rdmolfiles.SDMolSupplier).

        Args:
            all_dicts (dict): Dictionary with subdictionaries including a lists
            of identifiers (self.identifier_keyname).

        Raises:
            KeyError: if there is no identifier list.
            ValueError: if there is only one dataset.
        """
        sets = []
        if len(all_dicts.keys()) <= 2:
            raise ValueError(
                "Only one dataset is given. Shared molecules can only calculated when comparing at least 2 datasets!"
            )
        for single_dict in all_dicts:
            if single_dict == self.figure_dict_keyname:
                continue
            if not any(
                key == self.identifier_keyname
                for key in list(all_dicts[single_dict].keys())
            ):
                raise KeyError(
                    "A identifier list is needed, please run 'get_identifier_list_key'!"
                )
                # try:
                #     raise KeyError(
                #         "A identifier list is needed, please run 'get_identifier_list_key'!"
                #     )
                # except KeyError as e:
                #     logger.error(str(e), exc_info=True)
                #     raise
            single_set = set(all_dicts[single_dict][self.identifier_keyname])
            sets.append(single_set)
        shared_molecules = set.intersection(*sets)
        for single_dict in all_dicts:
            if single_dict == self.figure_dict_keyname:
                continue
            all_dicts[single_dict][self.shared_mols_keyname] = len(shared_molecules)
            all_dicts[single_dict][self.shared_mols_id_keyname] = shared_molecules
        logger.info(
            "Number of molecules found in all datasets: %d, identifier(s): %s",
            len(shared_molecules),
            shared_molecules,
        )
        logger.info(
            "Updated dictionary with '%s' and '%s'",
            self.shared_mols_keyname,
            self.shared_mols_id_keyname,
        )
        return

    def visualize_intersection(
        self, all_dicts: dict, data_type: str = "png"
    ) -> matplotlib.figure.Figure:
        """
        This function returns a Venn diagram of the intersection between the
        molecules in the subdictionaries of the given dictionary. Every
        subdictionary is represented as a circle and the overlaps between the
        circles indicate the molecules present in more than one subdictionary.
        (The function only works with two or three subdictionaries.)
        The intersection is based on the identifiers (string representation).
        The diagram is saved in an output folder.

        Args:
            all_dicts (dict): Dictionary of dictionaries with
            identifier_keyname.
            data_type (str): Data type for the exported image (default: png).
        Returns:
            fig (matplotlib.figure): Venn diagram

        Raises:
            ValueError: If there is only one or more than three sets to be
            compared.
            KeyError: If there is no identifier list.
        """
        sets = []
        for single_dict in all_dicts:
            if single_dict == self.figure_dict_keyname:
                continue
            if not any(
                key == self.identifier_keyname
                for key in list(all_dicts[single_dict].keys())
            ):
                raise KeyError(
                    "A identifier list is needed, please run 'get_identifier_list_key'!"
                )
            single_set = set(all_dicts[single_dict][self.identifier_keyname])
            sets.append(single_set)
        fig = plt.figure(figsize=(10, 10))
        if len(sets) == 3:
            venn = venn3(sets, set_labels=(all_dicts.keys()), alpha=0.5)
        elif len(sets) == 2:
            venn = venn2(sets, set_labels=(all_dicts.keys()), alpha=0.5)
        else:
            raise ValueError("Visualization only possible for two or three data sets!")
        plt.title("Intersection as Venn diagram", fontsize=20)
        for text in venn.set_labels:
            text.set_fontsize(15)
        for x in range(len(venn.subset_labels)):
            if venn.subset_labels[x] is not None:
                venn.subset_labels[x].set_fontsize(15)
        self._save_to_figure_dict(all_dicts, "intersection", fig, data_type=data_type)
        plt.close(fig)
        return fig

    # Section: Get descriptors and create plots

    def _get_descriptor_list(
        self,
        moleculeset: Chem.SDMolSupplier,
        descriptor: callable,
    ) -> List:
        """
        This function returns a list of descriptor values for all molecules in
        the given rdkit_mol objects (rdkit.Chem.rdmol.Mol or rdkit.Chem.
        rdmolfiles.SDMolSupplier) and a callable descriptor (e.g Descriptors.
        MolWt or rdMolDescriptors.CalcExactMolWt). (private method)

        Args:
            moleculeset (rdkit.Chem.rdmolfiles.SDMolSupplier or list[rdkit.Chem.
            rdmol.Mol])
            descriptor (callable): RDKit method that returns a molecular
            descriptor for a given molecule.

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
        This function updates the subdictionaries in the given dictionary with
        a list of descriptor values as a new key-value pair using
        _get_descriptor_list on the rdkit_mol Objects (rdkit.Chem.rdmol.Mol or
        rdkit.Chem.rdmolfiles.SDMolSupplier) in the subdictionaries.

        Args:
            all_dicts (dict): Dictionary with subdictionaries for every
            dataset, including the key 'self.import_keyname'.
            descriptor (callable): RDKit method that returns a molecular
            descriptor for a given molecule
            descriptor_list_keyname (str): Key name for the dictionary entry
            (should match the descriptor)
        """
        for single_dict in all_dicts:
            if single_dict == self.figure_dict_keyname:
                continue
            descriptor_list = self._get_descriptor_list(
                all_dicts[single_dict][self.import_keyname], descriptor
            )
            all_dicts[single_dict][descriptor_list_keyname] = descriptor_list
        logger.info("Updated dictionary with '%s'", descriptor_list_keyname)
        return

    def _get_discrete_descriptor_counts(
        self, all_dicts: dict, descriptor_list_keyname: str
    ) -> None:
        """
        This function updates the subdictionaries in the given dictionary with
        the binned descriptor values for a given descriptor value list with
        discrete values (e.g. number of H-Bond donors or acceptors). (private
        method)

        Args:
            all_dicts (dict): Dictionary with subdictionaries including a
            discrete descriptor value list.
            descriptor_list_keyname (str): Name of the descriptor list.
        """
        binned_descriptor_list_keyname = str("binned_" + descriptor_list_keyname)
        find_max = []
        for single_dict in all_dicts:
            if single_dict == self.figure_dict_keyname:
                continue
            if None in all_dicts[single_dict][descriptor_list_keyname]:
                find_max.append(
                    max(
                        [
                            descriptor_value
                            for descriptor_value in all_dicts[single_dict][
                                descriptor_list_keyname
                            ]
                            if descriptor_value is not None
                        ]
                    )
                )
            else:
                find_max.append(max(all_dicts[single_dict][descriptor_list_keyname]))
        maximum = max(find_max) + 1
        bins = pd.interval_range(start=0, end=maximum, freq=1, closed="left")
        for single_dict in all_dicts:
            if single_dict == self.figure_dict_keyname:
                continue
            counts = pd.value_counts(
                pd.cut(all_dicts[single_dict][descriptor_list_keyname], bins),
                sort=False,
            )
            all_dicts[single_dict][binned_descriptor_list_keyname] = counts
        logger.info("Updated the dictionary with '%s'", binned_descriptor_list_keyname)
        return

    def _get_continuous_descriptor_counts(
        self, all_dicts: dict, descriptor_list_keyname: str, width_of_bins: float = 10.0
    ) -> None:
        """
        This function updates the subdictionaries in the given dictionary with
        the binned descriptor values for a given descriptor value list with
        continuous values (e.g. molecular weight or logP values). The interval
        size of the bins can be chosen. (private method)

        Args:
            all_dicts (dict): Dictionary with subdictionaries including a
            continuous descriptor value list.
            descriptor_list_keyname (str): name of the descriptor list.
            width_of_bins (int, optional): Interval size for the bins (default:
            10)
        """
        binned_descriptor_list_keyname = str("binned_" + descriptor_list_keyname)
        find_min = []
        find_max = []
        for single_dict in all_dicts:
            if single_dict == self.figure_dict_keyname:
                continue
            find_min.append(
                min(
                    [
                        descriptor_value
                        for descriptor_value in all_dicts[single_dict][
                            descriptor_list_keyname
                        ]
                        if descriptor_value is not None
                    ]
                )
            )
            find_max.append(
                max(
                    [
                        descriptor_value
                        for descriptor_value in all_dicts[single_dict][
                            descriptor_list_keyname
                        ]
                        if descriptor_value is not None
                    ]
                )
            )
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
            if single_dict == self.figure_dict_keyname:
                continue
            counts = pd.value_counts(
                pd.cut(all_dicts[single_dict][descriptor_list_keyname], bins),
                sort=False,
            )
            all_dicts[single_dict][binned_descriptor_list_keyname] = counts
        logger.info("Updated the dictionary with '%s'", binned_descriptor_list_keyname)
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
    ) -> matplotlib.figure.Figure:
        """
        This function returns a bar-plot for a discrete descriptor with was
        previously binned.
        The plot is saved in an output folder as an image (data type can be
        chosen) and the data frame can also be saved as CSV file.

        args:
            all_dicts (dict): Dictionary with subdictionaries including a
            binned discrete descriptor.
            descriptor_list_keyname (str): Name of descriptor list for plotting.
            data_type (str): Data type for the exported image (default: png).
            save_dataframe (bool): Export dataframe as csv file or not
            (default: True).
            figsize (float, float): Width, height of the image in inches
            (default: 15, 7).
            fontsize_tick_labels (int): Fontsize of the labels on the ticks of
            the axis (default: 15).
            fontsize_legend (int): Fontsize of the legend (default: 15).
            fontsize_ylabel (int): Fontsize of the label of the y-axis
            (default: 20).
            fontsize_xlabel (int): Fontsize of the label of the x-axis
            (default: 20).
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
            if single_dict == self.figure_dict_keyname:
                continue
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
        descriptor_plot.legend(
            bbox_to_anchor=(1, 1), loc="upper left", fontsize=fontsize_legend
        )
        descriptor_plot.set_ylabel("Number of molecules", fontsize=fontsize_ylabel)
        descriptor_plot.set_xlabel(
            str(descriptor_list_keyname), fontsize=fontsize_xlabel
        )
        descriptor_plot.set_title(
            str("Distribution of " + descriptor_list_keyname),
            pad=20,
            fontsize=fontsize_title,
        )
        fig = descriptor_plot.figure
        self._save_to_figure_dict(
            all_dicts,
            keyname=("distribution_of_" + str(descriptor_list_keyname)),
            fig=fig,
            data_type=data_type,
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
    ) -> matplotlib.figure.Figure:
        """
        This function returns bar-plot for a continuous descriptor which was
        previously binned.
        The plot is saved in an output folder as an image (data type can be
        chosen) and the data frame can also be saved as CSV file.

        args:
            all_dicts (dict): Dictionary with subdictionaries including a
            binned continuous descriptor.
            descriptor_list_keyname (str): Name of descriptor list for plotting.
            data_type (str): Data type for the exported image (default: png).
            save_dataframe (bool): Export dataframe as csv file or not
            (default: True).
            fig_size (float, float): Width, height of the image in inches
            (default: 15, 7).
            fontsize_tick_labels (int): Fontsize of the labels on the ticks of
            the axis (default: 15).
            fontsize_legend (int): Fontsize of the legend (default: 15).
            fontsize_ylabel (int): Fontsize of the label of the y-axis
            (default: 20).
            fontsize_xlabel (int): Fontsize of the label of the x-axis
            (default: 20).
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
            if single_dict == self.figure_dict_keyname:
                continue
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
        descriptor_plot.legend(
            bbox_to_anchor=(1, 1), loc="upper left", fontsize=fontsize_legend
        )
        descriptor_plot.set_ylabel("Number of molecules", fontsize=fontsize_ylabel)
        descriptor_plot.set_xlabel(
            str(descriptor_list_keyname + " Intervals"), fontsize=fontsize_xlabel
        )
        descriptor_plot.set_title(
            str("Distribution of " + descriptor_list_keyname),
            pad=20,
            fontsize=fontsize_title,
        )
        fig = descriptor_plot.figure
        self._save_to_figure_dict(
            all_dicts,
            keyname=("distribution_of_" + str(descriptor_list_keyname)),
            fig=fig,
            data_type=data_type,
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
    ) -> matplotlib.figure.Figure:
        """
        This function updates the subdictionaries in the given dictionary with
        the binned descriptor values for a given descriptor value list. The
        values can either be continuous (binning with
        _get_continuous_descriptor_counts and plotted with
        _continuous_descriptor_plot) or discrete (binning with
        _get_discrete_descriptor_counts and plotted with
        _discrete_descriptor_plot).
        The created plots are saved in an output folder and the data frame can
        also be exported as CSV.

        Args:
            all_dicts (dict): Dictionary with subdictionaries including a
            descriptor value list.
            descriptor_list_keyname (str): Name of the descriptor list for
            binning and plotting.
            width_of_bins (int, optional): interval size for the bins for
            continuous values (default: 10).
            data_type (str): Data type for the exported image (default: png).
            save_dataframe (bool): Export dataframe as csv file or not
            (default: True).
            figsize (float, float): Width, height of the image in inches
            (default: 15, 7).
            fontsize_tick_labels (int): Fontsize of the labels on the ticks of
            the axis (default: 15).
            fontsize_legend (int): Fontsize of the legend (default: 15).
            fontsize_ylabel (int): Fontsize of the label of the y-axis
            (default: 20).
            fontsize_xlabel (int): Fontsize of the label of the x-axis
            (default: 20).
            fontsize_title (int): Fontsize of the title (default: 24).

        Raises:
            KeyError: if there is not the needed descriptor list.
            ValueError: if the descriptor values are not int or float and can
            therefore not be plotted.

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
                "A descriptor list ("
                + str(descriptor_list_keyname)
                + ") is needed for plotting. Please run 'get_descriptor_list_key'!"
            )
        elif type(all_dicts[first_dict][descriptor_list_keyname][0]) == int:
            self._get_discrete_descriptor_counts(all_dicts, descriptor_list_keyname)
            fig = self._discrete_descriptor_plot(
                all_dicts,
                descriptor_list_keyname,
                data_type,
                save_dataframe,
                figsize,
                fontsize_tick_labels,
                fontsize_legend,
                fontsize_ylabel,
                fontsize_xlabel,
                fontsize_title,
            )
        elif (
            type(all_dicts[first_dict][descriptor_list_keyname][0]) == float
            or type(all_dicts[first_dict][descriptor_list_keyname][0]) == np.float64
        ):
            self._get_continuous_descriptor_counts(
                all_dicts, descriptor_list_keyname, width_of_bins
            )
            fig = self._continuous_descriptor_plot(
                all_dicts,
                descriptor_list_keyname,
                data_type,
                save_dataframe,
                figsize,
                fontsize_tick_labels,
                fontsize_legend,
                fontsize_ylabel,
                fontsize_xlabel,
                fontsize_title,
            )
        else:
            raise ValueError(
                'Descriptor values should be "int" or "float" (numpy.float64) to be binned!'
            )
        return fig

    # Section: Check Lipinski Rule of 5 and visualization

    def _test_for_lipinski(self, moleculeset) -> List[int]:
        """
        This function returns a list with the number of Lipinski Rules broken
        for every molecule in the given molecule set.

        Args:
            moleculeset (rdkit.Chem.rdmolfiles.SDMolSupplier or list[rdkit.Chem.
            rdmol.Mol]): rdkit_mol Objects

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
            if Lipinski.NumHAcceptors(mol) > 10:
                rule_break += 1
            if Lipinski.NumHDonors(mol) > 5:
                rule_break += 1
            num_of_break.append(rule_break)
        return num_of_break

    def get_lipinski_key(self, all_dicts: dict) -> None:
        """
        This function updates the subdictionaries in the given dictionary with
        the list of the number of broken Lipinski Rules for every molecule
        (lipinski_list_keyname) and a summary of the broken rules
        (lipinski_summary_keyname) using _test_for_lipinski.

        Args:
            all_dicts (dict): Dictionary with subdictionaries including the key
            'self.import_keyname'.
        """
        for single_dict in all_dicts:
            if single_dict == self.figure_dict_keyname:
                continue
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
        logger.info(
            "Updated dictionary with '%s' and '%s'",
            self.lipinski_list_keyname,
            self.lipinski_summary_keyname,
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
    ) -> matplotlib.figure.Figure:
        """
        This function returns a bar plot for the number of molecules in every
        subdictionary breaking 0 to 4 Lipinski rules using the
        'lipinski_summary' key in the given dictionary. The plot is saved in an
        output folder (data type can be chosen) and the created data frame can
        also be exported as CSV.

        args:
            all_dicts (dict): Dictionary with subdictionaries including the key
            'self.lipinski_summary_keyname'.
            data_type (str): Data type for the exported image (default: png).
            save_dataframe (bool): Export dataframe as csv or not (default:
            True).
            fig_size (float, float): Width, height of the image in inches
            (default: 15, 7).
            fontsize_tick_labels (int): Fontsize of the labels on the ticks of
            the axis (default: 15).
            fontsize_legend (int): Fontsize of the legend (default: 15).
            fontsize_ylabel (int): Fontsize of the label of the y-axis
            (default: 20).
            fontsize_xlabel (int): Fontsize of the label of the x-axis
            (default: 20).
            fontsize_title (int): Fontsize of the title (default: 24).

        Raises:
            KeyError: if the Lipinski key is missing and therefore no plot can
            be generated.

        Returns:
            fig (matplotlib.figure.Figure): Plot
        """
        lipinski_df_dict = {"Number of broken rules": [0, 1, 2, 3, 4]}
        for single_dict in all_dicts:
            if single_dict == self.figure_dict_keyname:
                continue
            if not (
                any(
                    key == self.lipinski_summary_keyname
                    for key in list(all_dicts[single_dict].keys())
                )
            ):
                raise KeyError(
                    "Lipinski summary ("
                    + str(self.lipinski_summary_keyname)
                    + ") is needed for plotting! Please run 'get_lipinski_key'!"
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
        lipinski_plot.legend(
            bbox_to_anchor=(1, 1), loc="upper left", fontsize=fontsize_legend
        )
        lipinski_plot.set_ylabel("Number of molecules", fontsize=fontsize_ylabel)
        lipinski_plot.set_xlabel("Number of broken rules", fontsize=fontsize_xlabel)
        lipinski_plot.set_title(
            "Distribution of the number of broken Lipinski Rules",
            pad=20,
            fontsize=fontsize_title,
        )
        fig = lipinski_plot.figure
        self._save_to_figure_dict(all_dicts, "lipinski_plot", fig, data_type=data_type)
        plt.close(fig)
        return fig

    # Section: Scaffold analysis and plotting

    def _get_scaffold(
        self,
        moleculeset,
        number_of_structures: int = 5,
        structures_per_row: int = 5,
        image_size: int = 200,
        framework: bool = False,
        graph_framework: bool = False,
        normalize: bool = True,
    ) -> Tuple[PIL.PngImagePlugin.PngImageFile, list, pd.core.series.Series]:
        """
        This function creates a grid images of a chosen number of scaffolds/
        frameworks/graph framework for the molecules in a given rdkit_mol
        Object (rdkit.Chem.rdmol.Mol or rdkit.Chem.rdmolfiles.SDMolSupplier).
        The scaffolds/frameworks/graph framework are sorted by their frequency.
        The relative or absolute number of occurrence of a scaffold/framework/
        graph framework in the dataset is shown below each image.

        args:
            moleculeset (rdkit.Chem.rdmolfiles.SDMolSupplier or list[rdkit.Chem.
            rdmol.Mol]): rdkit_mol Objects.
            number_of_structures (int): Number of structures displayed in the
            grid images (default: 5).
            structures_per_row (int): Number of structures in every row of the
            grid image (default: 5).
            image_size (int): Size of the image for a single molecule (default:
            200).
            framework (bool): Remove terminal atoms with double bond (default:
            False).
            graph_framework (bool): Creating graph framework (default: False).
            normalize (bool): Using relative numbers for scaffold analysis
            (default: True).

        returns:
            structure_grid (PIL.PngImageFile): Grid image with most frequent
            scaffolds/frameworks/graph framework.
            structure_list (list): List of scaffolds/frameworks/graph framework
            for every molecule.
            structure_counts (pandas.Series): Absolute or relative frequency of
            each scaffold/frameworks/graph framework.

        raises:
            KekulizeException: If molecule of the rdkit_mol Object cannot be
            kekulized.
        """
        structure_list = []
        scaffold_list = []
        for mol in moleculeset:
            try:
                Chem.Kekulize(mol)
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            except KekulizeException:
                logger.info(
                    "Molecule can not be kekulized and will be excluded from scaffold analysis!"
                )
                continue
            # scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            scaffold_list.append(scaffold)
        if not framework and not graph_framework:
            for mol in scaffold_list:
                structure_list.append(Chem.MolToSmiles(mol))
        if framework or graph_framework:
            framework_list = []
            for mol in scaffold_list:
                to_remove = []
                Chem.Kekulize(mol)
                edit_mol = rdchem.RWMol(mol)
                for atom in edit_mol.GetAtoms():
                    if len(atom.GetNeighbors()) == 1:
                        to_remove.append(atom.GetIdx())
                for index in sorted(to_remove, reverse=True):
                    edit_mol.RemoveAtom(index)
                new_mol = edit_mol.GetMol()
                framework = Chem.RemoveHs(new_mol)
                framework_list.append(framework)
            for mol in framework_list:
                structure_list.append(Chem.MolToSmiles(mol))
        if graph_framework:
            structure_list.clear()
            graph_framework_list = []
            for mol in framework_list:
                try:
                    graph_framework_list.append(MurckoScaffold.MakeScaffoldGeneric(mol))
                except AtomValenceException:
                    index = framework_list.index(mol)
                    identifier = moleculeset[index]
                    logger.info(
                        "Graph framework can not be generated, molecule (%s, index %d) will be excluded from scaffold analysis!",
                        identifier,
                        index,
                    )
                continue
            for mol in graph_framework_list:
                structure_list.append(Chem.MolToSmiles(mol))
        structure_list = ["*" if mol == "" else mol for mol in structure_list]
        structure_counts = pd.Index(structure_list).value_counts(normalize=normalize)
        if len(structure_counts) < number_of_structures:
            number_of_structures = len(structure_counts)
        legend = [
            str(integer) for integer in (list(structure_counts)[:number_of_structures])
        ]
        smiles_list = list(structure_counts.keys())
        to_draw = []
        for index in range(number_of_structures):
            to_draw.append(Chem.MolFromSmiles(smiles_list[index]))
            if smiles_list[index] == "*":
                legend[index] = legend[index] + " (No rings/scaffolds)"
        structure_grid = Draw.MolsToGridImage(
            to_draw,
            maxMols=number_of_structures,
            molsPerRow=structures_per_row,
            subImgSize=(image_size, image_size),
            legends=legend,
            returnPNG=False,
        )
        return structure_grid, structure_list, structure_counts

    def draw_most_frequent_scaffolds(
        self,
        all_dicts: dict,
        number_of_structures: int = 5,
        structures_per_row: int = 5,
        image_size: int = 200,
        framework: bool = False,
        graph_framework: bool = False,
        normalize: bool = True,
        data_type: str = "png",
        figsize: Tuple[float, float] = [20.0, 20.0],
        fontsize_title: int = 24,
        fontsize_subtitle: int = 20,
    ) -> matplotlib.figure.Figure:
        """
        This function creates a grid images of a chosen number of scaffolds/
        frameworks/graph framework for every subdictionary in the given
        dictionary and shows them together. The scaffolds/frameworks/graph
        framework in each gird image are sorted by their frequency. The
        relative or absolute number of occurrence of a scaffold/framework/graph
        framework in the dataset of the respective subdictionary is shown below
        each image.

        args:
            all_dicts (dict): Dictionary with subdictionaries including the key
            'self.import_keyname'.
            number_of_scaffolds (int): Number of scaffolds displayed in the
            grid images (default: 5).
            scaffolds_per_row (int): Number of scaffolds in every row of the
            grid image (default: 5).
            image_size (int): Size of the image for a single molecule (default:
            200).
            framework (bool): Remove terminal atoms with double bond (default:
            False).
            graph_framework (bool): Creating graph framework (default: False).
            normalize (bool): Using relative numbers for scaffold analysis
            (default: True).
            data_type (str): Data type for the exported image (default: png).
            figsize (float, float): Width, height of the image in inches
            (default: 20, 20)
            fontsize_title (int): Fontsize of the title (default: 24).
            fontsize_subtitle (int): Fontsize of the subtitles (default: 20).

        returns:
            fig (PIL.PngImageFile): Grid images with most frequent scaffolds/
            frameworks/graph framework for each subdictionary.
        """
        image_list = []
        title_list = []
        for single_dict in all_dicts:
            if single_dict == self.figure_dict_keyname:
                continue
            title_list.append(single_dict)
            scaffolds = self._get_scaffold(
                all_dicts[single_dict][self.import_keyname],
                number_of_structures,
                structures_per_row,
                image_size,
                framework,
                graph_framework,
                normalize,
            )
            image_list.append(scaffolds[0])
            all_dicts[single_dict][self.scaffold_list_keyname] = scaffolds[1]
            all_dicts[single_dict][self.scaffold_summary_keyname] = (
                scaffolds[2].to_frame("frequency")
            ).rename_axis("scaffold SMILES")
        logger.info(
            "Updated dictionary with '%s' and '%s'",
            self.scaffold_list_keyname,
            self.scaffold_summary_keyname,
        )
        rows = len(image_list)
        fig = plt.figure(figsize=figsize)
        for j in range(0, rows):
            fig.add_subplot(rows, 1, j + 1)
            plt.axis("off")
            plt.imshow(image_list[j])
            plt.title(title_list[j], fontsize=fontsize_subtitle)
        fig.suptitle(
            "Most frequent scaffolds from the datasets", fontsize=fontsize_title
        )
        self._save_to_figure_dict(all_dicts, "scaffold_grid", fig, data_type=data_type)
        plt.close(fig)
        return fig

    # Section: Chemical space visualization

    def chemical_space_visualization(
        self,
        all_dicts: dict,
        fp_radius: int = 2,
        fp_bits: int = 512,
        dimension_reduction: str = "pca",
        interactive: bool = True,
    ):
        """
        This function returns a 2D visualization of the chemical space of the
        molecules in all datasets using the chemplot module.
        On basis of the calculated identifier (self.identifier_keyname) for
        every molecule a Extended Connectivity Fingerprint (ECFP) will be
        calculated with a definable fingerprint radius (fp_radius) and length
        (fp_size).
        Subsequent, the fingerprints are reduced to 2D for plotting. The
        dimension reduction can be done with PCA, UMAP or t-SNE and the plot
        can be interactive.

        Args:
            all_dicts (dict): Dictionary with subdictionaries including an
            identifier list (self.identifier_keyname).
            fp_radius (int): Radius of the Extended Connectivity Fingerprints
            (default: 2).
            fp_bits (int): Size of the Extended Connectivity Fingerprints
            (default: 2048).
            dimension_reduction (str): Method of dimension reduction (default:
            pca).
            interactive (bool): Creating an interactive plot or not (default:
            True).
        Raises:
            KeyError: if the identifier list is missing.
            ValueError: if the dimension reduction is not pca, tsne or umap.

        Returns:
            Chemical space visualization
        """
        all_mols_list = []
        target_list = []
        for single_dict in all_dicts:
            if single_dict == self.figure_dict_keyname:
                continue
            if not any(
                key == self.identifier_keyname
                for key in list(all_dicts[single_dict].keys())
            ):
                raise KeyError(
                    "A identifier list is needed, please run 'get_identifier_list_key'!"
                )
            for mol in all_dicts[single_dict][self.identifier_keyname]:
                all_mols_list.append(mol)
                target_list.append(single_dict)
        if all_mols_list[0].startswith("InChI="):  # check if identifier is InChI
            chem_space = cp.Plotter.from_inchi(
                all_mols_list,  # list of inchi strings which are used to get Extended Connectivity Fingerprint (alternative: smiles),
                target=target_list,  # corresponding list for inchi_list, shows which dataset the molecules belong to
                target_type="C",  # classification (classes are the datasets listed in the target_list)
                sim_type="structural",  # similarity solely based on structure (no property is taken into account)
                radius=fp_radius,
                nBits=fp_bits,
            )
        elif (
            len(all_mols_list[0]) == 27  # check if identifier is InChIKey
            and "-" in all_mols_list[0][14]
            and "-" in all_mols_list[0][25]
        ):
            all_mols_list.clear()
            for single_dict in all_dicts:
                if single_dict == self.figure_dict_keyname:
                    continue
                for mol in all_dicts[single_dict][self.import_keyname]:
                    inchi = Chem.MolToInchi(mol)
                    all_mols_list.append(inchi)
            chem_space = cp.Plotter.from_inchi(
                all_mols_list,  # list of inchi strings which are used to get Extended Connectivity Fingerprint (alternative: smiles),
                target=target_list,  # corresponding list for inchi_list, shows which dataset the molecules belong to
                target_type="C",  # classification (classes are the datasets listed in the target_list)
                sim_type="structural",  # similarity solely based on structure (no property is taken into account)
                radius=fp_radius,
                nBits=fp_bits,
            )
        else:
            chem_space = cp.Plotter.from_smiles(
                all_mols_list,  # list of smiles strings which are used to get Extended Connectivity Fingerprint (alternative: smiles),
                target=target_list,  # corresponding list for inchi_list, shows which dataset the molecules belong to
                target_type="C",  # classification (classes are the datasets listed in the target_list)
                sim_type="structural",  # similarity solely based on structure (no property is taken into account)
                radius=fp_radius,
                nBits=fp_bits,
            )
        if dimension_reduction == "pca":
            chem_space.pca()  # n_components, copy, whiten, svd_solver ...
        elif dimension_reduction == "tsne":
            chem_space.tsne(
                learning_rate=200.0, init="random"
            )  # n_components, perplexity, learning_rate, n_iter, init, random_state ...
        elif dimension_reduction == "umap":
            chem_space.umap()  # n_neighbors, min_dist, pca, random_state ...
        else:
            raise ValueError('dimension_reduction should be "pca", "tsne" or "umap"!')
        if not os.path.exists("output"):
            os.makedirs("output")
        if not interactive:
            fig = chem_space.visualize_plot().figure
            self._save_to_figure_dict(all_dicts, "chemical_space", fig)
            plt.close(fig)
        else:
            fig = chem_space.interactive_plot(show_plot=True)
            if not os.path.exists("output/interactive_chemical_space.html"):
                output_file("output/interactive_chemical_space.html")
            else:
                counter = 1
                file_name = str("interactive_chemical_space_")
                while os.path.exists("output/%s%d.html" % (file_name, counter)):
                    counter += 1
                output_file("output/%s%d.html" % (file_name, counter))
            save(fig)
        return fig

    # Section: Data export

    def export_single_dict_values(self, all_dicts: dict) -> None:
        """
        This function exports only the (not-binned) descriptor values for each
        dictionary according to the imported SDFile as a single csv file in the
        output folder.

        Args:
            all_dicts (dict): Dictionary with subdictionaries containing the
            calculated descriptor values.
        """
        for single_dict in all_dicts:
            if single_dict == self.figure_dict_keyname:
                continue
            new_dict = all_dicts[single_dict].copy()
            counter = 0
            for key in new_dict.copy():
                if key == self.import_keyname or self.duplicates_index_keyname:
                    new_dict.pop(key)
                elif type(new_dict[key]) == list:
                    counter += 1
                else:
                    new_dict.pop(key)
            to_export = pd.DataFrame(new_dict)
            filename = single_dict[:-4]
            to_export.to_csv(
                "output/descriptor_values_%s.csv" % (filename), index=False
            )
            logger.info("%s: %d exported descriptor values", single_dict, counter)
        return

    def export_figure_report(self, all_dicts: dict) -> None:
        """
        This function exports a pdf including some general information and all
        created figures form the figures subdictionary.

        Args:
            all_dicts (dict): Dictionary with subdictionary 'self.
            figure_dict_keyname' containing all created plots
        """
        first_dict = list(all_dicts.keys())[0]
        today = date.today()
        data = (
            ("Date:", str(today)),
            ("Data:", str(list(all_dicts.keys())[:-1])[1:-1]),
            ("Generated keys:", str(list(all_dicts[first_dict].keys()))[1:-1]),
            (
                "Generated figures:",
                str(list(all_dicts[self.figure_dict_keyname].keys()))[1:-1],
            ),
        )
        pdf = FPDF()
        pdf.set_font("Helvetica")
        pdf.add_page()
        pdf.set_font_size(30)
        pdf.cell(txt="Report CIDER", align="C", w=pdf.epw, border="B")
        pdf.set_font_size(12)
        pdf.set_y(50)
        for row in data:
            pdf.multi_cell(40, 7, row[0], new_y=YPos.LAST, new_x=XPos.RIGHT)
            pdf.multi_cell(150, 7, row[1], new_y=YPos.NEXT, new_x=XPos.LMARGIN)
        for fig in all_dicts["figures"]:
            pdf.add_page()
            reader = io.BytesIO()
            all_dicts["figures"][fig].savefig(reader, bbox_inches="tight", format="png")
            fig = Image.open(reader)
            if fig.height <= fig.width:
                pdf.image(reader, w=pdf.epw)
            if fig.height > fig.width:
                if fig.height / fig.width < 1.4:
                    pdf.image(reader, w=pdf.epw)
                if fig.height / fig.width > 1.4:
                    pdf.image(reader, h=pdf.eph)
            fig.close()
        pdf.output("output/cider_report.pdf")
        logger.info("'cider_report.pdf' exported")
        return
