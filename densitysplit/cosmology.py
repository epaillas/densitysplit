import numpy as np
from scipy.integrate import quad, simps
from scipy.special import hyp2f1, legendre, eval_legendre
from scipy.interpolate import InterpolatedUnivariateSpline

class Cosmology:
    '''Class for cosmology calculations.

    Parameters
    ----------
    omega_m : float, optional
        Matter density parameter.
    h : float, optional
        Hubble parameter.
    c : float, optional
        Speed of light in km/s.
    
    Attributes
    ----------
    H0 : float
        Hubble constant in km/s/Mpc.
    h : float
        Dimensionless Hubble parameter.
    c : float
        Speed of light in km/s.
    omega_m : float
        Matter density parameter.
    omega_lambda : float
        Dark energy density parameter.
    zgrid : array_like
        Redshift grid used for interpolation.
    rgrid : array_like
        Comoving distance grid used for interpolation.
        '''
    def __init__(self,
                omega_m=0.308,
                h=0.676,
                c=299792.458):

        omega_lambda = 1.0 - omega_m
        H0 = 100.0

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
        '''Calculate Hubble parameter as a function of redshift.
        
        Parameters
        ----------
        z : float
            Redshift.
            
            Returns
            -------
            float
            Hubble parameter in units of 100 km/s/Mpc.
            '''
        return self.H0 * self.Ez(z)

    def ComovingDistance(self, z):
        '''Calculate comoving distance as a function of redshift.

        Parameters
        ----------
        z : float
            Redshift.

        Returns
        -------
        float
            Comoving distance in Mpc/h.
        '''
        return np.interp(z, self.zgrid, self.rgrid)

    def AngularDiameterDistance(self, z):
        '''Calculate angular diameter distance as a function of redshift.

        Parameters
        ----------
        z : float
            Redshift.

        Returns
        -------
        float
            Angular diameter distance in Mpc/h.
        '''
        return np.interp(z, self.zgrid, self.rgrid) / (1 + z)

    def Redshift(self, r):
        '''Calculate redshift as a function of comoving distance.

        Parameters
        ----------
        r : float
            Comoving distance in Mpc/h.

        Returns
        -------
        float
            Redshift.
        '''
        return np.interp(r, self.rgrid, self.zgrid)

    def GrowthFactor(self, z):
        '''Calculate linear growth factor as a function of redshift.

        Parameters
        ----------
        z : float
            Redshift.

        Returns
        -------
        float
            Linear growth factor.
        '''
        az = 1. / (1 + z)
        growth = az ** 2.5 * np.sqrt(self.omega_lambda + self.omega_m * az ** (-3.)) * \
                hyp2f1(5. / 6, 3. / 2, 11. / 6, -(self.omega_lambda * az ** 3.) / self.omega_m) / \
                hyp2f1(5. / 6, 3. / 2, 11. / 6, -self.omega_lambda / self.omega_m)
        return growth

    def GrowthRate(self, z):
        '''Calculate linear growth rate as a function of redshift.

        Parameters
        ----------
        z : float
            Redshift.

        Returns
        -------
        float
            Linear growth rate.
        '''
        f = ((self.omega_m * (1 + z)**3.) / (self.omega_m * (1 + z)**3 + self.omega_lambda))**0.55
        return f