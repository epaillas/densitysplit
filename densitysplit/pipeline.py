import numpy as np
from pypower import CatalogMesh
from pandas import qcut
from densitysplit import filters
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
    def __init__(self, data_positions, boxsize=None,
        data_weights=None, randoms_positions=None,
        randoms_weights=None):
        self.data_positions = data_positions
        self.boxsize = boxsize

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

    def get_mesh(self):
        """
        Get the CatalogMesh object.

        Returns
        -------
        mesh : CatalogMesh
            CatalogMesh object.
        """
        if self.boxsize is None:
            mesh  = CatalogMesh(data_positions=self.data_positions,
                data_weights=self.data_weights, cellsize=self.cellsize,
                resampler=self.resampler, position_type='pos',
                randoms_positions=self.randoms_positions,
                randoms_weights=self.randoms_weights, boxpad=self.boxpad,
                interlacing=0)
        else:
            mesh  = CatalogMesh(data_positions=self.data_positions,
                data_weights=self.data_weights, boxsize=self.boxsize,
                cellsize=self.cellsize, resampler=self.resampler, 
                position_type='pos', interlacing=0)
        return mesh


    def get_box_randoms(self, nrandoms, seed=42):
        """
        Get random points inside the box.
        
        Parameters
        ----------
        nrandoms : int
            Number of random points to generate.
        seed : int, optional
            Random seed.

        Returns
        -------
        randoms : array_like
            Random points inside the box.
        """
        np.random.seed(seed)
        randoms = np.random.rand(nrandoms, 3) * self.boxsize
        return randoms


    def get_density(self, smooth_radius, cellsize, compensate=True,
        resampler='cic', sampling='randoms', sampling_positions=None,
        boxpad=2.0, filter_shape='tophat', ran_min=0.01):
        """
        Get the density field.

        Parameters
        ----------
        smooth_radius : float
            Radius of the smoothing filter.
        cellsize : float
            Size of the cells in the mesh.
        compensate : bool, optional
            Compensate for the shot noise.
        resampler : str, optional
            Resampling method.
        sampling : str, optional
            Sampling method.
        sampling_positions : array_like, optional
            Positions where the density field should be sampled. If not 
            provided, the randoms are used.
        boxpad : float, optional
            Padding of the box.
        filter_shape : str, optional
            Shape of the smoothing filter.

        Returns
        -------
        density : array_like
            Density field at the sampling positions.
        sampling_positions : array_like
            Positions where the density field was sampled.
        """
        self.cellsize = cellsize
        self.boxpad = boxpad
        self.resampler = resampler
        self.mesh = self.get_mesh()

        data_mesh = self.mesh.to_mesh(field='data', compensate=compensate)
        data_mesh = data_mesh.r2c().apply(
            getattr(filters, filter_shape)(r=smooth_radius))
        data_mesh = data_mesh.c2r()

        if self.boxsize is None:
            randoms_mesh = self.mesh.to_mesh(field='data-normalized_randoms',
                compensate=compensate)
            randoms_mesh = randoms_mesh.r2c().apply(
                getattr(filters, filter_shape)(r=smooth_radius))
            randoms_mesh = randoms_mesh.c2r()
            threshold = ran_min * np.mean(randoms_mesh)
            mask = randoms_mesh > threshold

            density_mesh = data_mesh - randoms_mesh
            density_mesh[mask] /= randoms_mesh[mask]
            density_mesh[~mask] = 0.0
            shift = self.mesh.boxsize / 2 - self.mesh.boxcenter
        else:
            nmesh = self.mesh.nmesh[0]
            norm = sum(self.data_weights)
            density_mesh = data_mesh/(norm/(nmesh**3)) - 1
            shift = 0

        if sampling_positions is not None:
            self.density = density_mesh.readout(sampling_positions + shift)
            self.sampling_positions = sampling_positions
        else:
            if sampling == 'randoms':
                if self.boxsize is None:
                    self.density = density_mesh.readout(self.randoms_positions + shift,
                        resampler=resampler)
                    self.sampling_positions = self.randoms_positions
                else:
                    nrandoms = 5 * len(self.data_positions)
                    randoms = self.get_box_randoms(nrandoms)
                    self.density = density_mesh.readout(randoms + shift, resampler=resampler)
                    self.sampling_positions = randoms
            elif sampling == 'data':
                self.density = density_mesh.readout(self.data_positions + shift,
                    resampler=resampler)
                self.sampling_positions = self.data_positions
            else:
                raise ValueError('Invalid sampling method. If sampling_positions '
                    f'is not provided, need to set sampling to "data" or "randoms"')
            
        return self.density

    def get_quantiles(self, nquantiles, return_density=False):
        """
        Get the quantiles of the density field.

        Parameters
        ----------
        nquantiles : int
            Number of quantiles.
        return_density : bool, optional
            Return the density field at the quantiles.

        Returns
        -------
        quantiles : array_like
            Quantiles of the density field.
        density_quantiles : array_like
            Density field at the quantiles.
        """
        quantiles_idx = qcut(self.density, nquantiles, labels=False)
        quantiles = []
        for i in range(nquantiles):
            quantiles.append(self.sampling_positions[quantiles_idx == i])
        self.quantiles = quantiles
        if return_density:
            density_quantiles = []
            for i in range(nquantiles):
                density_quantiles.append(
                    np.mean(self.density[quantiles_idx == i])
                )
            density_quantiles = np.asarray(density_quantiles, dtype=float)
            return quantiles, density_quantiles
        return quantiles
