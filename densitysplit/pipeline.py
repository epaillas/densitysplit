import numpy as np
from nbodykit.lab import ArrayCatalog
from nbodykit.filters import TopHat
from pandas import qcut


class DensitySplit:
    def __init__(self, data_positions, boxsize, nmesh):

        self.data_positions = data_positions
        self.boxsize = boxsize
        self.nmesh = nmesh

    def get_data_mesh(self):
        cat_size = len(self.data_positions)
        dset = np.empty(cat_size, dtype=[('Position', ('f8', 3))])
        dset['Position'][:, 0] = self.data_positions[:, 0]
        dset['Position'][:, 1] = self.data_positions[:, 1]
        dset['Position'][:, 2] = self.data_positions[:, 2]
        cat = ArrayCatalog(dset)
        self.data_mesh = cat.to_mesh(BoxSize=self.boxsize,
            Nmesh=self.nmesh, position='Position')

    def get_randoms_positions(self, nrandoms, seed=42):
        np.random.seed(seed)
        self.random_positions = np.random.rand(nrandoms, 3) * self.boxsize
        return self.random_positions

    def get_density(self, smooth_radius):
        self.get_data_mesh()
        filtered_mesh = self.data_mesh.apply(TopHat(r=smooth_radius))

        painted_mesh = filtered_mesh.paint(mode='real')
        density_mesh = painted_mesh - 1
        self.density = density_mesh.readout(self.randoms_positions)
        return self.density

    def get_quantiles(self, nquantiles):
        quantiles_idx = qcut(self.density, nquantiles, labels=False)
        quantiles = []
        for i in range(nquantiles):
            quantiles.append(self.randoms_positions[quantiles_idx == i])
        self.quantiles = quantiles
        return self.quantiles


