import numpy as np
from pyrecon import RealMesh
from pandas import qcut
import sys


class DensitySplit:
    """
    Class to split a set of points into quantiles based on the local 
    density field.

    Parameters
    ----------
    data_positions : array_like
        Positions of the data points.
    boxsize : float, optional
        Size of the box. If not provided, randoms are required.
    data_weights : array_like, optional
        Weights of the data points. If not provided, all points are 
        assumed to have the same weight.
    randoms_positions : array_like, optional
        Positions of the random points. If not provided, boxsize must 
        be provided.
    randoms_weights : array_like, optional
        Weights of the random points. If not provided, all points are 
        assumed to have the same weight.

    Attributes
    ----------
    data_positions : array_like
        Positions of the data points.
    boxsize : float
        Size of the box.
    data_weights : array_like
        Weights of the data points.
    randoms_positions : array_like
        Positions of the random points.
    randoms_weights : array_like
        Weights of the random points.
    mesh : CatalogMesh
        CatalogMesh object.
    density : array_like
        Density field at the sampling positions.
    sampling_positions : array_like
        Positions where the density field was sampled.
    """
    def __init__(self, data_positions, boxsize=None, boxcenter=None,
        data_weights=None, randoms_positions=None, randoms_weights=None,
        cellsize=None, wrap=False, nthreads=None):
        self.data_positions = data_positions
        self.boxsize = boxsize
        self.boxcenter = boxcenter
        self.cellsize = cellsize
        self.wrap = wrap
        self.nthreads = nthreads

        if data_weights is not None:
            self.data_weights = data_weights
        else:
            self.data_weights = np.ones(len(data_positions))

        if boxsize is None:
            if randoms_positions is None:
                raise ValueError(
                    'boxsize is set to None, but randoms were not provided.')
            if randoms_weights is None:
                self.randoms_weights = np.ones(len(randoms_positions))
            else:
                self.randoms_weights = randoms_weights
            self.randoms_positions = randoms_positions


    def get_density_mesh(self, sampling_positions, smoothing_radius):
        """
        Get the overdensity field.

        Parameters
        ----------
        smooth_radius : float
            Radius of the smoothing filter.
        sampling_positions : array_like
            Positions where the density field should be sampled.
        Returns
        -------
        density : array_like
            Density field at the sampling positions.
        """
        self.data_mesh = RealMesh(boxsize=self.boxsize, cellsize=self.cellsize,
                                  boxcenter=self.boxcenter, nthreads=self.nthreads)
        self.data_mesh.assign_cic(self.data_positions, wrap=self.wrap)
        self.data_mesh.smooth_gaussian(smoothing_radius, engine='fftw', save_wisdom=True)
        if self.boxsize is None:
            self.randoms_mesh = RealMesh(boxsize=self.boxsize, cellsize=self.cellsize,
                                         boxcenter=self.boxcenter, nthreads=self.nthreads)
            self.randoms_mesh.assign_cic(self.randoms_positions, wrap=self.wrap)
            self.randoms_mesh.smooth_gaussian(smoothing_radius, engine='fftw',)
            sum_data, sum_randoms = np.sum(self.data_mesh.value), np.sum(self.randoms_mesh.value)
            alpha = sum_data * 1. / sum_randoms
            self.delta_mesh = self.data_mesh - alpha * self.randoms_mesh
            mask = self.randoms_mesh > 0
            self.delta_mesh[mask] /= alpha * self.randoms_mesh[mask]
            self.delta_mesh[~mask] = 0.0
            del self.data_mesh
            del self.randoms_mesh
        else:
            self.delta_mesh = self.data_mesh / np.mean(self.data_mesh) - 1.
            del self.data_mesh
        self.delta = self.delta_mesh.read_cic(sampling_positions) 
        self.sampling_positions = sampling_positions
        return self.delta


    def get_density_paircount(self, smooth_radius, sampling_positions, filter_shape='Tophat', nthreads=1):
        """
        Get the overdensity field using pair counting.

        Parameters
        ----------
        smooth_radius : float
            Radius of the smoothing filter.
        sampling_positions : array_like
            Positions where the density field should be sampled.
        filter_shape : str, optional
            Shape of the smoothing filter.
        nthreads : int, optional
            Number of threads to use.
        
        Returns
        -------
        density : array_like
            Density field at the sampling positions.
        """
        from julia.api import Julia
        from os import path, environ
        environ['JULIA_NUM_THREADS'] = f'{nthreads}'
        jl = Julia(compiled_modules=False)
        from julia import Main
        module_path = path.join(path.dirname(__file__),
            'fastmodules', 'count_pairs.jl')
        jl.eval(f'include("{module_path}")')
        if self.boxsize is None:
            Main.positions1 = sampling_positions.T
            Main.weights1 = np.ones(len(sampling_positions))
            Main.positions2 = self.data_positions.T
            Main.weights2 = self.data_weights
            Main.smooth_radius = smooth_radius
            D1D2 = jl.eval(f"count_pairs_survey_{filter_shape.lower()}(positions1, positions2, weights1, weights2, smooth_radius)")
            Main.positions2 = self.randoms_positions.T 
            Main.weights2 = self.randoms_weights 
            D1R2 = jl.eval(f"count_pairs_survey_{filter_shape.lower()}(positions1, positions2, weights1, weights2, smooth_radius)")
            D1R2 *= np.sum(self.data_weights) / np.sum(self.randoms_weights)
        else:
            Main.positions1 = sampling_positions.T
            Main.weights1 = np.ones(len(sampling_positions)) 
            Main.positions2 = self.data_positions.T
            Main.weights2 = self.data_weights
            Main.boxsize = self.boxsize
            Main.smooth_radius = smooth_radius
            D1D2 = jl.eval(f"count_pairs_box_{filter_shape.lower()}(positions1, positions2, weights1, weights2, boxsize, smooth_radius)")
            bin_volume = 4/3 * np.pi * smooth_radius ** 3
            mean_density = np.sum(self.data_weights) / (self.boxsize[0]*self.boxsize[1]*self.boxsize[2])
            D1R2 = bin_volume * mean_density * np.ones(len(sampling_positions))
        self.delta = D1D2 - D1R2
        mask = D1R2 > 0
        self.delta[mask] /= D1R2[mask]
        self.delta[~mask] = 0
        masked_fraction = 1 - np.sum(mask) / len(mask)
        self.sampling_positions = sampling_positions
        return self.delta


    def get_quantiles(self, nquantiles, return_idx=False):
        """
        Get the quantiles of the density field.

        Parameters
        ----------
        nquantiles : int
            Number of quantiles.
        return_idx : bool, optional
            Whether to return index of the quantile of each query point.

        Returns
        -------
        quantiles : array_like
            Quantiles of the density field.
        quantiles_idx : array_like, optional
            Index of the quantile of each query point.
        """
        quantiles_idx = qcut(self.delta, nquantiles, labels=False)
        quantiles = []
        for i in range(nquantiles):
            quantiles.append(self.sampling_positions[quantiles_idx == i])
        self.quantiles = quantiles
        if return_idx:
            return quantiles, quantiles_idx
        return quantiles
