import numpy as np
from scipy.io import FortranFile
from astropy.io import fits
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.special import eval_legendre
from scipy.integrate import simps


def sky_to_cartesian(data, cosmology):
    '''
    Converts a catalogue in sky coordinates
    to comoving cartesian coordinates. This
    assumes the usual convention where
    RA = [0, 360] and DEC = [-90, 90].

    Parameters:  data: 2D array_like
                 Array containing the catalogue in sky coordinates.

                 cosmology: Cosmology object
                 Object describing the cosmology for coordinate conversion.
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
    '''
    Converts a catalogue in comoving cartesian
    coordinates to sky coordinates. The returned
    array is in the usual convetion RA = [0, 360]
    and DEC = [-90, 90].

    Parameters:  data: 2D array_like
                 Array containing the catalogue in cartesian coordinates.

                cosmology: Cosmology object
                Object describing the cosmology for coordinate conversion.
    '''

    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]

    dis = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    dec = np.arctan2(np.sqrt(x ** 2 + y ** 2), z) * 180 / np.pi
    ra = np.arctan2(y, x) * 180 / np.pi
    redshift = cosmology.Redshift(dis)

    ind = ra > 360
    ra[ind] -= 360
    ind = ra < 0
    ra[ind] += 360

    dec = 90 - dec

    cout = np.c_[ra, dec, redshift]
    return cout


def save_as_unformatted(data, filename):
    '''
    Saves a numpy array as an unformatted
    Fortran 90 file that can be handled by
    this package's numerical routines.

    Parameters:  data: ND array_like
                 Array to be saved.

                 filename: str
                 Name of the output file.
    '''
    data = np.asarray(data)

    nrows, ncols = np.shape(data)
    f = FortranFile(filename, 'w')
    nrows, ncols = np.shape(data)
    f.write_record(nrows)
    f.write_record(ncols)
    f.write_record(data)
    f.close()


def read_boss_fits(filename, columns=['RA', 'DEC', 'Z']):
    '''
    Reads a BOSS DR12 clustering catalogue in fits format
    and returns it as a numpy array in the desired order.

    Parameters:  filename: str
                 Name of the input file.

                 columns: list of strings
                 List containing the strings corresponding to the
                 columns that want to be extracted.
    '''

    with fits.open(filename) as hdul:
        cat = hdul[1].data

    cout = cat[columns[0]]

    for column in columns[1:]:
        cout = np.c_[cout, cat[column]]

    return cout


def compute_boss_weights(
    weight_fkp, weight_cp=None, weight_noz=None,
    weight_systot=None, is_random=False
):
    '''
    Computes the total weight of each galaxy from the
    BOSS sample by combining individual weights.

    Parameters:  weight_fkp: 1D array_like
                 FKP weights.

                 weight_cp: 1D array_like
                 Close pair weights.

                 weight_noz: 1D array_like
                 Redshift failure weights.

                 weight_systot: 1D array_like
                 Total imaging systematic weights.

                 is_random: boolean
                 Set to True if the weights are from a
                 random caatalogue.

    '''
    weight = weight_fkp

    if not is_random:
        weight *= weight_systot * (weight_cp + weight_noz - 1)

    return weight


def read_array_2d(filename):
    '''
    Read a two dimensional from an ascii
    file.

    Parameters:  filename: str
                 Name of the ascii file containing
                 the array
    '''
    data = np.genfromtxt(filename)
    dim1 = np.unique(data[:, 0])
    dim2 = np.unique(data[:, 1])

    vary_dim2 = False
    if data[0, 0] == data[1, 0]:
        vary_dim2 = True

    result = np.zeros([len(dim1), len(dim2)])
    counter = 0
    if vary_dim2:
        for i in range(len(dim1)):
            for j in range(len(dim2)):
                result[i, j] = data[counter, 2]
                counter += 1
    else:
        for i in range(len(dim2)):
            for j in range(len(dim1)):
                result[j, i] = data[counter, 2]
                counter += 1
    return dim1, dim2, result


def get_multipole(ell, s, mu, xi_smu):
    '''
    Calculate the ell multipole moment
    from the redshift-space correlation
    function xi_smu.

    Parameters:  ell: int
                 Multipole moment to calculate.

                 s: 1D array_like
                 Radial bins of the correlation function.

                 mu: 1D array_like
                 Mu bins of the correlation function, where
                 mu is the cosine of the angle with respect
                 to the line of sight.

                 xi_smu: 2D array_like


    '''
    multipole = np.zeros(xi_smu.shape[0])
    if mu.min() < 0:
        factor = 2
        mumin = -1
    else:
        factor = 1
        mumin = 0
    for i in range(xi_smu.shape[0]):
        mufunc = InterpolatedUnivariateSpline(mu, xi_smu[i, :], k=3, ext=3)
        xaxis = np.linspace(mumin, 1, 1000)
        lmu = eval_legendre(ell, xaxis)
        yaxis = mufunc(xaxis) * (2 * ell + 1) / factor * lmu
        multipole[i] = simps(yaxis, xaxis)
    return s, multipole
