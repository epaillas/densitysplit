import argparse
from survey.split import DensitySplit
import yaml
from os import path

parser = argparse.ArgumentParser()
parser.add_argument(
    '-p', '--parameter_file', help='YAML parameter file',
    default='parameters/ds_survey.yaml')
parser.add_argument('-f', '--smoothing_radius')
parser.add_argument('-H', '--handle')
parser.add_argument('-t', '--data_fn')
parser.add_argument('-r', '--randoms_fn')
parser.add_argument('-q', '--nquantiles')
parser.add_argument('-s', '--nseeds')
parser.add_argument('-S', '--seeds_fn')
parser.add_argument('-T', '--nthreads')

args = parser.parse_args()

# load parameter from YAML file
with open(args.parameter_file) as file:
    params = yaml.full_load(file)

# command-line arguments take precedence over configuration file
if args.data_fn is not None:
    params['input']['data_fn'] = args.data_fn

if args.randoms_fn is not None:
    params['input']['randoms_fn'] = args.randoms_fn

if args.handle is not None:
    params['output']['handle'] = args.handle

if args.seeds_fn is not None:
    params['output']['seeds_fn'] = args.seeds_fn

if args.nthreads is not None:
    params['algorithm']['nthreads'] = args.nthreads

if args.nquantiles is not None:
    params['algorithm']['nquantiles'] = args.nquantiles

# # rewrite parameters with command-line arguments
# for name in vars(args):
#     if args.__dict__[name] is not None:
#         params[name] = args.__dict__[name]

pipe = DensitySplit(params)

if not path.isfile(params['seeds_filename']):
    pipe.generate_seeds()

pipe.split_densities()
