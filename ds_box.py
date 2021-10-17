import argparse
from py.split import DensitySplit
from py.cosmology import Cosmology
import yaml

parser = argparse.ArgumentParser()
parser.add_argument(
    '-p', '--parameter_file', help='YAML parameter file',
    default='default_params.yaml')

args = parser.parse_args()

# load parameter from YAML file
parms = yaml.full_load(args.parameter_file)

# rewrite parameters with command-line arguments
for name in vars(args):
    if name is not None:
        parms[name] = args.__dict__[name]


pipe = DensitySplit(parms)

pipe.generated_seeds()
pipe.split_densities()
