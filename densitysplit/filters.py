import numpy as np


class TopHat(object):
    '''Top-hat filter in Fourier space
    adapted from https://github.com/bccp/nbodykit/

    Parameters
    ----------
    r : float
        the radius of the top-hat filter
    '''
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


class Gaussian(object):
    '''Gaussian filter in Fourier space

    Parameters
    ----------
    r : float
        the radius of the Gaussian filter
    '''
    def __init__(self, r):
        self.r = r

    def __call__(self, k, v):
        r = self.r
        k2 = sum(ki ** 2 for ki in k)
        return np.exp(- 0.5 * k2 * r**2) * v
