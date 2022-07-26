import numpy as np
from nbodykit.lab import ArrayCatalog
from nbodykit.filters import TopHat
from pandas import qcut


class DensitySplit:
    def __init__(self, data_positions, boxsize,
        data_weights=None):

        self.data_positions = data_positions
        self.boxsize = boxsize

        if data_weights is not None:
            self.data_weights = data_weights
        else:
            self.data_weights = np.ones(len(data_positions))

    def get_data_mesh(self):
        cat_size = len(self.data_positions)
        dset = np.empty(cat_size, dtype=[('Position',
            ('f8', 3)), ('Weight', ('f8'))])
        dset['Position'][:, 0] = self.data_positions[:, 0]
        dset['Position'][:, 1] = self.data_positions[:, 1]
        dset['Position'][:, 2] = self.data_positions[:, 2]
        dset['Weight'] = self.data_weights
        cat = ArrayCatalog(dset)
        self.data_mesh = cat.to_mesh(BoxSize=self.boxsize,
            Nmesh=self.nmesh, position='Position', weight='Weight',
            compensated=self.compensated, window=self.window)

    def get_randoms_positions(self, nrandoms, seed=42):
        np.random.seed(seed)
        self.randoms_positions = np.random.rand(nrandoms, 3) * self.boxsize
        return self.randoms_positions

    def get_density(self, smooth_radius, nmesh, compensated=False,
        window='cic', sampling='randoms'):
        self.nmesh = nmesh
        self.window = window
        self.compensated = compensated
        self.get_data_mesh()
        filtered_mesh = self.data_mesh.apply(TopHat(r=smooth_radius))
        painted_mesh = filtered_mesh.paint(mode='real')
        density_mesh = painted_mesh - 1
        if sampling == 'randoms':
            if not hasattr(self, 'randoms_positions'):
                self.get_randoms_positions()
            self.density = density_mesh.readout(self.randoms_positions)
        elif sampling == 'data':
            self.density = density_mesh.readout(self.data_positions)
        self.sampling = sampling
        return self.density

    def get_quantiles(self, nquantiles):
        quantiles_idx = qcut(self.density, nquantiles, labels=False)
        quantiles = []
        for i in range(nquantiles):
            if self.sampling = 'randoms':
                quantiles.append(self.randoms_positions[quantiles_idx == i])
            elif self.sampling == 'data':
                quantiles.append(self.data_positions[quantiles_idx == i])
        self.quantiles = quantiles
        return self.quantiles


