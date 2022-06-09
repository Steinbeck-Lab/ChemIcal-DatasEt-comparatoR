#!/usr/bin/env python

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="CIDER",
    version="0.0.1-dev",
    author="Hannah Busch",
    author_email="hannah.busch@uni-jena.de",
    maintainer="Jonas Schaub, Otto Brinkhaus, Kohulan Rajan",
    maintainer_email="jonas.schaub@uni-jena.de, henning.brinkhaus@uni-jena.de, kohulan.rajan@uni-jena.de",
    description="ChemIcal DatasEt compaRator (CIDER) is a Python package which primarily utilizes RDKit to compare and visualize different chemical compounds from two different datasets.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/hannbus/ChemIcal_DatasEt_compaRator",
    packages=setuptools.find_packages(),
    license="MIT",
    install_requires=[
        "notebook",
        "seaborn==0.11.2",
        "matplotlib==3.5.1",
        "chemplot==1.2.0",
        "matplotlib_venn==0.11.6",
        "FPDF==1.7.2",
        "rdkit-pypi",
    ],
    package_data={"CIDER": ["assets/*.*", "assets/*/*.*", "assets/*/*/*.*"]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.5",
)