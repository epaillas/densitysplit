### Redshift-space distortions with split densities

A package for running the redshift-space distortions with split densities pipeline presented in https://arxiv.org/abs/2101.09854.

The interface of the code is written in Python 3, while the engine that takes care of intensive calculations is written in the Julia programming language, so you must have it installed along with the packages listed below.

### Julia requirements

  - `CellListMap` 
  - `StaticArrays` 
  - `LinearAlgebra` 

### Python requirements

  - `pyjulia` >= 0.5.7
  - `numpy` >= 1.20.1

### Installation

Under the main directory, install the package with `python setup.py install`. 

### Usage

The easiest way to run the pipeline is through the stand-alone code provided in the `bin` directory, which takes a `params.yaml` configuration file as an input. Alternatively, you can import the code as a Python package and run individual parts of the pipeline as you wish. Examples of this are shown in the `examples` directory.
