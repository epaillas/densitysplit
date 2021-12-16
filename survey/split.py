import sys
from os import path
import re
import numpy as np
import subprocess
from astropy.io import fits
from scipy.io import FortranFile
from .cosmology import Cosmology
from .utilities import (sky_to_cartesian, cartesian_to_sky,
    save_as_unformatted, read_unformatted)


class DensitySplit:
    """
    Class to perform the density split algorithm on
    a cosmological survey, as in 
    https://arxiv.org/abs/1810.02864.
    """

    def __init__(self, parameters: dict):
        """
        Initialize class.

        Args:
            parameter_file: YAML file containing configuration parameters.
        """

        self.params = parameters

        allowed_formats = ['rdz', 'xyz']
        self.catalogs = {}
        self.positions = {}
        self.weights = {}
        columns = []
        input_format = {}
        input_fns, output_fns = {}, {}
        self.input_fns_unf = {}

        # Whether nbar is provided by randoms catalog,
        # or nbar is assumed uniform
        allowed_selection_functions = ['uniform', 'randoms', '']
        selection_function = self.params['algorithm'].pop(
                'selection_function', '').lower()
        if selection_function not in allowed_selection_functions:
            raise Exception('Unknown input selection function {}. Choices are {}'.format(selection_function, allowed_selection_functions))
        # First check what we have in input/output
        input_fns, output_fns = {}, {}
        for name in ['data', 'randoms']:
            tmp_fn = self.params['input'].get('{}_fn'.format(name), None)
            if tmp_fn is None:
                if name == 'randoms':
                    if selection_function == 'randoms':
                        raise Exception('Please provide randoms catalog.')
                    # No randoms provided and no instruction on selection function, defaults to uniform nbar
                    if not selection_function:
                        # logger.info('No randoms provided.')
                        selection_function = 'uniform'
                else:
                    raise Exception('Please provide data catalog.')
            else: # We've got a file name!
                input_fns[name] = tmp_fn
            tmp_fn = self.params['output'].get('{}_fn'.format(name), None)
            if tmp_fn is not None:
                # Check that requested catalog can be supplied given input
                if name not in input_fns:
                    raise Exception('Cannot output {} catalog if not provided as input.'.format(name))
                output_fns[name] = tmp_fn
        # Randoms catalog provided and no instruction on selection function, defaults to nbar from randoms
        if not selection_function:
            selection_function = 'randoms'
        # logger.info('Using {} selection function.'.format(selection_function))
        self.selection_function = selection_function

        for name in ['data', 'randoms']:
            input_fns[name] = self.params['input'].get('{}_fn'.format(name), None)

        for name in input_fns:
            input_format[name] = None
            for format in allowed_formats:
                cols = self.params['input'].get('{}_{}'.format(format,name),None)
                if cols is not None:
                    if input_format[name] is not None:
                        raise Exception('Cannot use two different input formats 1')
                    input_format[name] = format
                    position_columns = cols 

            if input_format[name] is None:
                for format in allowed_formats:
                    cols = self.params['input'].get(format,None)
                    if cols is not None:
                        # Check whether e.g. 'rdz' but 'xyz' has been specified previously
                        if input_format[name] is not None:
                            raise Exception('Cannot use two different input formats 2')
                        input_format[name] = format
                        position_columns = cols
            format = input_format[name]
            # No format 'xyz', 'xyz_data', 'rdz_randoms', ... found
            if format is None:
                raise Exception('Unknown input format. Choices are {}'.format(allowed_formats))
            position_columns = self.make_list(position_columns)

            mask_str = self.params['input'].get(
                'mask_{}'.format(name), self.params['input'].get('mask', None))
            mask_str, mask_columns = self.decode_eval_str(mask_str)
            weight_str = self.params['input'].get(
                 'weights_{}'.format(name), self.params['input'].get('weights', None))
            weight_str, weight_columns = self.decode_eval_str(weight_str)

            for cols in [position_columns, mask_columns, weight_columns]:
                columns += cols

            columns = self.remove_duplicates(columns)

            with fits.open(input_fns[name]) as hdul:
                tmp = hdul[1].data
        
            catalog = {}
            for col in columns:
                catalog[col] = tmp[col]

            mask = eval(mask_str, catalog)
            
            for col in columns:
                catalog[col] = catalog[col][mask]

            # prepare Cartesian positions from input columns
            if format == 'rdz':
                if not len(position_columns) == 3:  # RA, DEC, Z
                    raise Exception('Format rdz requires 3 position columns')
                omega_m = self.params['cosmology']['omega_m']
                cosmo = Cosmology(omega_m=omega_m)
                self.positions[name] = sky_to_cartesian(
                    catalog[position_columns[0]], catalog[position_columns[1]],
                    catalog[position_columns[2]], cosmo
                )

            else:  # format == 'xyz'
                if len(position_columns) == 3:  # X, Y, Z
                    self.positions[name] = np.array(
                        [catalog[col] for col in position_columns]
                    ).T
                elif len(position_columns) == 1:  # single array of shape (N, 3)
                    self.positions[name] = catalog[position_columns[0]]
                else:
                    raise Exception('Format xyz requires 1 or 3 position columns')

            # calculate weights
            self.weights[name] = eval(weight_str, catalog)

            # prepare unformatted files for Fortran routines
            self.input_fns_unf[name] = input_fns[name].split('.fits')[0] + '.unf'
            cout = np.c_[self.positions[name], self.weights[name]]
            save_as_unformatted(cout, self.input_fns_unf[name])

            self.catalogs[name] = catalog


    def generate_seeds(self):

        if self.selection_function == 'randoms'
            # sample from randoms file
            sampling_data = self.positions['randoms']
            nseeds = self.params['algorithm']['nseeds']
            idx = np.random.choice(
                len(sampling_data), size=nseeds, replace=False
            )
            seeds = sampling_data[idx]
            save_as_unformatted(seeds, self.params['output']['seeds_fn'])
        else:
            raise Exception('Uniform selection function, but boxsize not provided')

    def split_densities(self):
        """
        Split the random seeds according to the local
        galaxy density.
        """

        estimator = 'DP'  # hard-coded for now
        use_weights = '1'

        handle = self.params['output']['handle']
        filter_type = 'tophat' 
        smoothing_radius = self.params['algorithm']['smoothing_radius']

        filter_fn = f'{handle}_{filter_type}{smoothing_radius}.unf'

        gridmin = -5000  # hard coded for now
        gridmax = 5000  # hard coded for now
        rmin = 0.0
        rmax = smoothing_radius
        tracer_format = 'unformatted'

        binpath = path.join(path.dirname(__file__),
                            'bin', 'tophat_filter.exe')

        cmd = [
            binpath,
            self.params['output']['seeds_fn'],
            self.input_fns_unf['data'],
            self.input_fns_unf['randoms'],
            self.input_fns_unf['randoms'],
            filter_fn,
            str(rmin),
            str(rmax),
            str(smoothing_radius),
            str(self.params['mesh']['nmesh']),
            str(gridmin),
            str(gridmax),
            estimator,
            str(self.params['algorithm']['nthreads']),
            use_weights,
            tracer_format
            ]

        subprocess.call(cmd)

        # open filter file
        with FortranFile(filter_fn, 'r') as f:
            smoothed_delta = f.read_ints()[0]
            smoothed_delta = f.read_reals(dtype=np.float64)

        # split seeds in quantiles
        quantiles = self.params['algorithm']['nquantiles']
        seeds = read_unformatted(self.params['output']['seeds_fn'])
        cosmo = Cosmology(omega_m=self.params['cosmology']['omega_m'])
        sky_seeds = cartesian_to_sky(seeds[:,0], seeds[:,1], seeds[:,2], cosmo)
        seeds[:, :3] = sky_seeds
        nseeds = len(seeds)
        idx = np.argsort(smoothed_delta)
        sorted_seeds = seeds[idx]

        binned_seeds = {}
        for i in range(1, quantiles + 1):
            binned_seeds['DS{}'.format(i)] = sorted_seeds[
                int((i - 1) * nseeds / quantiles):int(i * nseeds / quantiles)]

            weights = np.ones(len(binned_seeds[f'DS{i}']))
            cout = np.c_[binned_seeds[f'DS{i}'], weights]
            quantiles_filename = f'{handle}_DS{i}_seeds_sky.dat'
            np.savetxt(quantiles_filename, cout, fmt=4*'%10.5f')



    def decode_eval_str(self, s):
        # Change ${col} => col, and return list of columns
        if s is None:
            return '', []
        toret = str(s)
        columns = []
        for replace in re.finditer('(\${.*?})',s):
            value = replace.group(1)
            col = value[2:-1]
            toret = toret.replace(value,col)
            if col not in columns: columns.append(col)
        return toret, columns


    def remove_duplicates(self, cols):
        # Remove duplicate column names
        toret = []
        for col in cols:
            if col not in toret: toret.append(col)
        return toret

    def make_list(self, cols):
        # Turn single column name to list of column names
        if cols is None: return []
        if isinstance(cols, str): cols = [cols]
        return cols

 
