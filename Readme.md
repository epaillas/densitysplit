### DensitySplit

A package for running the redshift-space distortions with split densities pipeline presented in https://arxiv.org/abs/2101.09854.

### Installation:

Under the main directory, type `make` to compile the source files.

Requirements: 
  - `Python` >= 3.0.0
  - `GCC` >= 4.2.1
  - `Scipy` >= 1.5.0

Earlier versions of those software packages might work as well, but they have not been tested.

### Usage:

The package is divided into two separate codes, depending on whether you want to run the pipeline on a periodic simulation or a galaxy survey. You can specify the input parameters through a YAML configuration file (examples can be found under the `parameters` directory), and then run the code as

`python ds_box.py -p configuration_file.yaml`

if you want to run it on a cosmological simulation, or

`python ds_survey -p configuration_file.yaml`

if you are working with a galaxy survey. You can also overwrite the YAML configuration parameters by providing specific command-line arguments, which is useful if you want to run the code many times with different input parameters. You can find example scripts of this under the `examples` directory.

