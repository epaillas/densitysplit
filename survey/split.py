import sys
from os import path
import numpy as np
import subprocess
from scipy.io import FortranFile
from .cosmology import Cosmology
from .utilities import sky_to_cartesian


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
    
    def generate_seeds(self):
        """
        Generates the initial random points for the
        density split algorithm. These can either be
        subsampled from an input file, or be randomly
        selected from a uniform distribution within
        a rectangular boundary.
        """
        np.random.seed(0)  # set random seed for reproducibility

        handle = self.params['handle']
        seeds_filename = f'{handle}_seeds.dat'
        self.params['seeds_filename'] = seeds_filename
        if not path.isfile(seeds_filename):
            if self.params['seeds_method'] == 'subsampling':
                sampling_filename = self.params['sampling_filename']
                if not path.isfile(sampling_filename):
                    raise FileNotFoundError(f'{sampling_filename} not found.')

                # only text file supported for now
                sampling_data = np.genfromtxt(sampling_filename)
                idx = np.random.choice(
                    len(sampling_data), size=self.params['nseeds'],
                    replace=False
                )
                seeds = sampling_data[idx]

                # convert to comoving coordinates if necessary
                if self.params['convert_seeds']:
                    if self.params['omega_m'] is None:
                        raise ValueError('If convert_seeds is True, '
                                         'omega_m needs to be specified.')
                cosmo = Cosmology(omega_m=self.params['omega_m'])
                seeds = sky_to_cartesian(seeds, cosmo)

            elif self.params['seeds_method'] == 'uniform':
                x = np.random.uniform(
                    self.params['seeds_xmin'],
                    self.params['seeds_xmax'],
                    self.params['nseeds']
                )
                y = np.random.uniform(
                    self.params['seeds_ymin'],
                    self.params['seeds_ymax'],
                    self.params['nseeds']
                )
                z = np.random.uniform(
                    self.params['seeds_zmin'],
                    self.params['seeds_zmax'],
                    self.params['nseeds']
                )

                seeds = np.c_[x, y, z]
            else:
                sys.exit('Sampling method not recognized')

            # add constant weights
            weights = np.ones(len(seeds))
            seeds = np.c_[seeds, weights]

            # save to file
            handle = self.params['handle']
            seeds_filename = f'{handle}_seeds.dat'
            np.savetxt(seeds_filename, seeds)

    def split_densities(self):
        """
        Split the random seeds according to the local
        galaxy density.
        """

        self.params['cf_estimator'] = 'DP'  # hard-coded for now

        handle = self.params['handle']
        filter_type = 'tophat' 
        filter_radius = self.params['filter_radius']

        # check if files exists
        for filename in [
            self.params['seeds_filename'],
            self.params['tracers_filename'],
            self.params['randoms_filename']
        ]:
            if not path.isfile(filename):
                raise FileNotFoundError(f'{filename} not found')

        if self.params['use_weights'] is True:
            use_weights = 1
        else:
            use_weights = 0

        if self.params['convert_tracers']:
            if self.params['omega_m'] is None:
                raise ValueError('If convert_tracers is True, '
                                 'omega_m needs to be specified.')

            cosmo = Cosmology(omega_m=self.params['omega_m'])

            # convert tracer coordinates
            converted_file = f'{handle}_tracers_sky.dat'
            self.params['tracers_filename'] = converted_file 
            if not path.isfile(converted_file):
                tracers = np.genfromtxt(self.params['tracers_filename'])
                ra = tracers[:, 0]
                dec = tracers[:, 1]
                z = tracers[:, 2]
                sky_tracers = np.c_[ra, dec, z]
                comoving_tracers = sky_to_cartesian(sky_tracers, cosmo)
                tracers[:, :3] = comoving_tracers
                np.savetxt(self.params['tracers_filename'], tracers)

            # convert random coordinates
            converted_file = f'{handle}_randoms_sky.dat'
            self.params['randoms_filename'] = converted_file 
            if not path.isfile(converted_file):
                randoms = np.genfromtxt(self.params['randoms_filename'])
                ra = randoms[:, 0]
                dec = randoms[:, 1]
                z = randoms[:, 2]
                sky_randoms = np.c_[ra, dec, z]
                comoving_randoms = sky_to_cartesian(sky_randoms, cosmo)
                randoms[:, :3] = comoving_randoms
                np.savetxt(self.params['randoms_filename'], randoms)

        # figure out bounding box in comoving coordinates

        filter_filename = f'{handle}_{filter_type}{filter_radius}.dat'
        if not path.isfile(filter_filename):

            filter_unformatted = f'{handle}_{filter_type}{filter_radius}.unf'

            gridmin = -5000  # hard coded for now
            gridmax = 5000  # hard coded for now
            rmin = 0.0
            rmax = filter_radius
            tracer_format = 'ascii'

            binpath = path.join(path.dirname(__file__),
                                'bin', 'tophat_filter.exe')

            cmd = [
                binpath,
                self.params['seeds_filename'],
                self.params['tracers_filename'],
                self.params['randoms_filename'],
                self.params['randoms_filename'],
                filter_unformatted,
                str(rmin),
                str(rmax),
                str(filter_radius),
                str(self.params['ngrid']),
                str(gridmin),
                str(gridmax),
                self.params['cf_estimator'],
                str(self.params['nthreads']),
                str(use_weights),
                tracer_format
                ]

            subprocess.call(cmd)

            # open filter file
            with FortranFile(filter_unformatted, 'r') as f:
                smoothed_delta = f.read_ints()[0]
                smoothed_delta = f.read_reals(dtype=np.float64)

            np.savetxt(filter_filename, smoothed_delta)

        # split seeds in quantiles
        smoothed_delta = np.genfromtxt(filter_filename)
        quantiles = self.params['quantiles']
        seeds = np.genfromtxt(self.params['seeds_filename'])
        nseeds = len(seeds)
        idx = np.argsort(smoothed_delta)
        sorted_seeds = seeds[idx]

        binned_seeds = {}
        for i in range(1, quantiles + 1):
            binned_seeds['DS{}'.format(i)] = sorted_seeds[
                int((i - 1) * nseeds / quantiles):int(i * nseeds / quantiles)]

            cout = binned_seeds[f'DS{i}']
            quantiles_filename = f'{handle}_DS{i}_seeds.dat'
            np.savetxt(quantiles_filename, cout)
