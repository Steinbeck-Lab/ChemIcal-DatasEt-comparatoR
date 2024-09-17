"""
CIDER Python Package.

============================
ChemIcal DatasEt comparatoR (CIDER) is a Python package and ready-to-use Jupyter Notebook workflow
which primarily utilizes RDKit to compare two or more chemical structure datasets (SD files) in different aspects,
e.g. size, overlap, molecular descriptor distributions, chemical space clustering, etc., most of which can be visually inspected in the notebook.
"""

__version__ = "1.0.0"

__all__ = [
    "CIDER",
]


from .cider import ChemicalDatasetComparator
