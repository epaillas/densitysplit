import numpy as np
from scipy.io import FortranFile
from astropy.io import fits
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.special import eval_legendre
from scipy.integrate import simps


def sky_to_cartesian(data, cosmology):
    '''Converts ra, dec, redshift to cartesian coordinates.

    Parameters
    ----------
    data : array_like
        Array of ra, dec, redshift.
    cosmology : Cosmology
        Cosmology object.

    Returns
    -------
    cout : array_like
        Array of x, y, z coordinates.
    '''
    ra = data[:, 0]
    dec = data[:, 1]
    redshift = data[:, 2]

    if np.any(dec > 90):
        dec = 90 - dec

    dist = cosmology.ComovingDistance(redshift)
    x = dist * np.cos(dec * np.pi / 180) * np.cos(ra * np.pi / 180)
    y = dist * np.cos(dec * np.pi / 180) * np.sin(ra * np.pi / 180)
    z = dist * np.sin(dec * np.pi / 180)

    cout = np.c_[x, y, z]
    return cout


def cartesian_to_sky(data, cosmology):
    '''Converts cartesian coordinates to ra, dec, redshift.

    Parameters
    ----------
    data : array_like
        Array of x, y, z coordinates.
    cosmology : Cosmology
        Cosmology object.

    Returns
    -------
    cout : array_like
        Array of ra, dec, redshift.
    '''
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    dist = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    dec = 90 - np.degrees(np.arccos(z / dist))
    ra = np.degrees(np.arctan2(y, x))
    ra[ra < 0] += 360
    redshift = cosmology.Redshift(dist)
    cout = np.c_[ra, dec, redshift]
    return cout