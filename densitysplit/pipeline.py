import numpy as np
from pypower import CatalogMesh
from pandas import qcut


class DensitySplit:
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
        if self.boxsize is None:
            mesh  = CatalogMesh(data_positions=self.data_positions,
                data_weights=self.data_weights, nmesh=self.nmesh,
                resampler=self.resampler, position_type='pos',
                randoms_positions=self.randoms_positions,
                randoms_weights=self.randoms_weights, boxpad=1.0,
                interlacing=0)
        else:
            mesh  = CatalogMesh(data_positions=self.data_positions,
                data_weights=self.data_weights, boxsize=self.boxsize,
                nmesh=self.nmesh, resampler=self.resampler, position_type='pos',
                interlacing=0)
        return mesh


    def get_box_randoms(self, nrandoms, seed=42):
        np.random.seed(seed)
        self.randoms_positions = np.random.rand(nrandoms, 3) * self.boxsize


    def get_density(self, smooth_radius, nmesh, compensate=True,
        resampler='cic', sampling='randoms'):
        self.nmesh = nmesh
        self.resampler = resampler
        self.mesh = self.get_mesh()

        if self.boxsize is None:
            data_mesh = self.mesh.to_mesh(field='data', compensate=compensate)
            data_mesh = data_mesh.r2c().apply(TopHat(r=smooth_radius))
            data_mesh = data_mesh.c2r()

            randoms_mesh = self.mesh.to_mesh(field='data-normalized_randoms',
                compensate=compensate)
            is_empty = randoms_mesh == 0.0 
            randoms_mesh = randoms_mesh.r2c().apply(TopHat(r=smooth_radius))
            randoms_mesh = randoms_mesh.c2r()
            
            density_mesh = data_mesh / randoms_mesh - 1
            density_mesh[is_empty] = np.nan

        else:
            data_mesh = self.mesh.to_mesh(field='data', compensate=compensate)
            data_mesh = data_mesh.r2c().apply(TopHat(r=smooth_radius))
            data_mesh = data_mesh.c2r()

            norm = sum(self.data_weights)
            density_mesh = data_mesh/(norm/(nmesh**3)) - 1

        if sampling == 'randoms':
            if not hasattr(self, 'randoms_positions'):
                nrandoms = 5 * len(self.data_positions)
                self.get_box_randoms(nrandoms)
            self.density = density_mesh.readout(self.randoms_positions)
        elif sampling == 'data':
            self.density = density_mesh.readout(self.data_positions)
        self.sampling = sampling
        return self.density

    def get_quantiles(self, nquantiles):
        quantiles_idx = qcut(self.density, nquantiles, labels=False)
        quantiles = []
        for i in range(nquantiles):
            if self.sampling == 'randoms':
                quantiles.append(self.randoms_positions[quantiles_idx == i])
            elif self.sampling == 'data':
                quantiles.append(self.data_positions[quantiles_idx == i])
        # self.quantiles = np.array(quantiles, dtype=object)
        self.quantiles = quantiles
        return quantiles

class TopHat(object):
    # adapted from https://github.com/bccp/nbodykit/
    def __init__(self, r):
        self.r = r
    def __call__(self, k, v):
        r = self.r
        k = sum(ki ** 2 for ki in k) ** 0.5
        kr = k * r
        with np.errstate(divide='ignore', invalid='ignore'):
            w = 3 * (np.sin(kr) / kr ** 3 - np.cos(kr) / kr ** 2)
        w[k == 0] = 1.0
        return w * v


