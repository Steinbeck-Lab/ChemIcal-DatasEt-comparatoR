"""
CIDER Python Package.
============================
ChemIcal DatasEt compaRator (CIDER) is a Python package which primarily utilizes RDKit
to compare and visualize different chemical compounds from two different datasets.
"""

__version__ = "0.0.1-dev"

__all__ = [
    "CIDER",
]


from .cider import (
    import_as_data_dict,
    get_number_of_molecules,
    draw_molecules,
    get_database_id,
    get_identifier_list,
    get_identifier_list_key,
    get_duplicate_key,
    get_shared_molecules_key,
    visualize_intersection,
    get_descriptor_list,
    get_descriptor_list_key,
    get_value_from_id,
    get_discrete_descriptor_counts,
    get_continuous_descriptor_counts,
    discrete_descriptor_plot,
    continuous_descriptor_plot,
    descriptor_counts_and_plot,
    test_for_lipinski,
    get_lipinski_key,
    lipinski_plot,
    chemical_space_visualization,
    export_single_dict_values,
    export_all_picture_pdf,
)
