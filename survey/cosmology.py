import numpy as np
from scipy.integrate import quad, simps
from scipy.special import hyp2f1, legendre, eval_legendre
from scipy.interpolate import InterpolatedUnivariateSpline

class Cosmology:
    '''
    Class for cosmological calculations.
    '''
    def __init__(self,
                omega_m=0.308,
                h=0.676,
                c=299792.458):

        omega_lambda = 1.0 - omega_m
        H0 = 100.0

        # we generate a mapping between distance and redshift
        # over which we can interpolate later to speed up 
        # calculations
        zgrid = np.linspace(0, 4, 1000)
        rgrid = np.zeros_like(zgrid)
        for i in range(len(zgrid)):
            rgrid[i] = quad(lambda x: (c/H0) / np.sqrt(omega_m * (1 + x) ** 3
                      + omega_lambda), 0, zgrid[i])[0]

        self.H0 = H0
        self.h = h
        self.c = c
        self.omega_m = omega_m
        self.omega_lambda = omega_lambda
        self.zgrid = zgrid
        self.rgrid = rgrid

    def Ez(self, z):
        ez = np.sqrt(self.omega_m * (1 + z) ** 3 + self.omega_lambda)
        return ez

    def HubbleParameter(self, z):
        return self.H0 * self.Ez(z)

    # comoving distance in Mpc/h
    def ComovingDistance(self, z):
        return np.interp(z, self.zgrid, self.rgrid)

    # angular diameter distance in Mpc/h
    def AngularDiameterDistance(self, z):
        return np.interp(z, self.zgrid, self.rgrid) / (1 + z)

    # redshift at a given comoving distance
    def Redshift(self, r):
        return np.interp(r, self.rgrid, self.zgrid)

    # growth factor at a given redshift
    def GrowthFactor(self, z):
        az = 1. / (1 + z)
        growth = az ** 2.5 * np.sqrt(self.omega_lambda + self.omega_m * az ** (-3.)) * \
                hyp2f1(5. / 6, 3. / 2, 11. / 6, -(self.omega_lambda * az ** 3.) / self.omega_m) / \
                hyp2f1(5. / 6, 3. / 2, 11. / 6, -self.omega_lambda / self.omega_m)
        return growth

    # linear growth rate at a given redshift
    def GrowthRate(self, z):
        f = ((self.omega_m * (1 + z)**3.) / (self.omega_m * (1 + z)**3 + self.omega_lambda))**0.55
        return f