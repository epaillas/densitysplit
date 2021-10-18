import argparse
from survey.split import DensitySplit
from survey.cosmology import Cosmology
import yaml
from os import path

parser = argparse.ArgumentParser()
parser.add_argument(
    '-p', '--parameter_file', help='YAML parameter file',
    default='parameters/ds_survey.yaml')
parser.add_argument('-f', '--filter_radius')
parser.add_argument('-H', '--handle')
parser.add_argument('-t', '--tracers_filename')
parser.add_argument('-r', '--randoms_filename')
parser.add_argument('-q', '--quantiles')
parser.add_argument('-s', '--nseeds')
parser.add_argument('-S', '--seeds_filename')
parser.add_argument('-w', '--use_weights')
parser.add_argument('-T', '--nthreads')
parser.add_argument('--sampling_filename')


# TODO add all required input CL parameters

args = parser.parse_args()

# load parameter from YAML file
with open(args.parameter_file) as file:
    params = yaml.full_load(file)

# rewrite parameters with command-line arguments
for name in vars(args):
    if args.__dict__[name] is not None:
        params[name] = args.__dict__[name]


params['seeds_filename'] = f'{params["handle"]}_seeds.dat'

pipe = DensitySplit(params)

if not path.isfile(params['seeds_filename']):
    pipe.generate_seeds()

pipe.split_densities()
