"""
Helper functions

author: Max Briel

"""
import numba
import numpy as np
from kea.constants import *
from scipy import integrate
from scipy.optimize import fminbound

def calculateLB(z):
    """Calculates the lookback time from the redshift.

    Parameters
    ----------
    z : float
        The redshift

    Returns
    -------
    float
        The lookback time at the redshift

    """
    def func(x):
        E = np.sqrt(MILLENIUM_OMEGAM*(1+x)**3 +MILLENIUM_OMEGAK*(1+x)**2 + MILLENIUM_OMEGAL)
        return 1/((1+x)*E)
    return t_HUBBLE *integrate.quad(func, 0, z)[0]/(60*60*24*365.2388526*1e9)


def approximateZ(LB):
    """Approximate the Redshift from a given lookback.

    Parameters
    ----------
    LB : float
        Lookback Time

    Returns
    -------
    float
        Approximated Redshift

    """
    zmin = 1e-8
    zmax = 1000
    ztol = 1e-8
    maxfun = 500
    f = lambda z: abs(calculateLB(z)- LB)
    zbest, resval, ierr, ncall = fminbound(f,zmin, zmax, maxfun=maxfun, full_output=1, xtol=ztol)
    return zbest




@numba.njit
def _numba_get_bin(x, edges):
    """
    The numba version to get the bin number.
    Used to speed up the calculation.

    Parameters
    ----------
    x : float
        value where you want to know the bin number
    edges: array
        An array with the edges of the histogram
    Returns
    -------
    int
        Position

    """
    if x < edges[0] or x > edges[-1]:
        raise Exception("x outside of range")
    out = 0
    for i,val in enumerate(edges):
        if val > x:
            out = i-1
            break
    if x == edges[-1]:
        out = len(edges)-2
    return out



@numba.njit
def _numba_integral(x1, x2, edges, values, bin_width):
    """The numba wrapper around a basic integration for the histogram.

    Parameters
    ----------
    x1 : float
        lower bound of the integration
    x2 : float
        upper bound of the integration
    edges : array
        The histogram bin edges
    values : array
        The values in each bin of the histogram
    bin_width : array
        The width of each bin in the historgam

    Returns
    -------
    float
        The integral between **x1** and **x2**

    """
    lower_bin = _numba_get_bin(x1, edges)
    upper_bin = _numba_get_bin(x2, edges)

    total = 0
    if lower_bin == upper_bin:
        total = values[lower_bin] * (x2-x1)
    else:

        ledge = lower_bin+1
        total += (values[lower_bin] * (edges[ledge]-x1))
        total += (values[upper_bin] * (x2 - edges[upper_bin]))

        if ledge < upper_bin:
            total += np.sum(bin_width[ledge:upper_bin]*values[ledge:upper_bin])

    return total



@numba.njit
def _numba_calculate_2D_SFRD(snapnum, stellar_mass, Z_stellar_mass, stellar_formation_rate):
    """
    Numba function to calculate the metallicity and calculate the sfr in each bin
    over time and metallicity.

    Parameters
    ----------
    snapnum : numpy array
        An array containing the snapshot numbers of the galaxies
    stellar_mass : numpy array
        An array containing the stellar mass of the galaxies
    Z_stellar_mass : numpy array
        An array containing the metals in stellar mass of the galaxies
    stellar_formation_rate : numpy array
        An array containing the stellar formation rate of the galaxies

    Returns
    -------
    2D numpy array
        A 2D numpy array with the metallicity and the snapnumbers (length: 64).
    """
    np_SFR = np.zeros((len(BPASS_NUM_METALLICITIES),64), dtype=np.float64)

    for snap, sm, Zsm, sfr in zip(snapnum, stellar_mass, Z_stellar_mass, stellar_formation_rate):
        if sm == 0:
            continue
        metallicity = Zsm/sm

        Zbin = np.argmin(np.abs(BPASS_NUM_METALLICITIES - metallicity))
        np_SFR[Zbin,snap] += sfr

    return np_SFR/((62.5/0.73)**3)



@numba.njit
def _numba_event_rate_loop(Z_values, mass_per_bin, edges, np_rates, DTD_width):
    """
    Numba function to calculate the event rates for the given rates, SFR,
    and mass per bin

    Parameters
    ----------
    Z_values: numpy array
        An array containing the metallicity values at the edges of the
        final binning.
    mass_per_bin : numpy array
        An array containig the amount of mass per bin in the final binning.
    np_rates : 2D numpy array
        A 2D array containig the different metallicities over time
        in BPASS binning. Format np_rates[metallicity][time]
    DTD_width : numpy array
        An array containing the bin widths from the Delay Time Distribution.

    Returns
    -------
    An non-normalised event rate. Number of events per bin.

    """
    Z = (Z_values[1:] + Z_values[:1])/2
    Z_per_bin = np.array([np.argmin(np.abs(i - BPASS_NUM_METALLICITIES)) for i in Z])
    event_rate = np.zeros(len(mass_per_bin))
    for count in range(len(mass_per_bin)):
            t = edges[count+1]
            for j in range(0,count+1):
                p1 = t - edges[j]
                p2 = t - edges[j+1]
                bin_events = _numba_integral(p2, p1, BPASS_LINEAR_TIME_EDGES, np_rates[Z_per_bin[count]], DTD_width)
                event_rate[j] += bin_events*mass_per_bin[count]
    return event_rate



@numba.njit()
def _property_summation(data):
    """
    Perform a summation over the galxies in the data set.

    Parameters
    ----------
    Data : 2D numpy array
        A 2D numpy array, where each row is a galaxy containing its snapshot
        number, stellar formation rate, stellar mass, and metals stellar mass.

    Returns
    -------
    tuple with numpy arrays
        A tuple containing 3 numpy arrays: sfr, SM, ZSM over the snapshots.
    """
    sfr = np.zeros(64)
    sm = np.zeros(64)
    Zsm = np.zeros(64)

    for i in range(len(data)):
        row = data[i]
        i = int(row[0])
        sfr[i] += row[1]
        sm[i] += row[2]
        Zsm[i] += row[3]
    return (sfr, sm, Zsm)
