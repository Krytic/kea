#
# Functions to help extract the star formation rate & get an event rate
#
# Author: Max Briel
#
import numpy as np
from kea.hist import histogram
from scipy import interpolate, integrate
from kea.constants import *
from scipy.optimize import fminbound
import kea.helpers
import numba
import pandas as pd

def getSFRD(cosmological_simulation, time_relations, length, h):
    """Extracts the total star formation rate density from a cosmological model
    without taking the metallicity into account.

    Parameters
    ----------
    cosmological_simulation : pandas DataFrame
        The cosmological simulation where to extract the SFR from. Needs to
        contain a column with 'sfr'.
    time_relations : pandas DataFrame
        contains the relation between snapshot number (snapnum) and lookback time
    length : float
        the length of the simulation
    h : float
        the Hubble parameter

    Returns
    -------
    numpy array
        An array with a length of 64 with the stellar formation rate.
    """
    SFR = []
    for i in time_relations["snapNum"]:
        data = cosmological_simulation.loc[
                                    cosmological_simulation["snapnum"] == i]
        if len(data) == 0:
            SFR.append(0)
        else:
            SFR.append(data["sfr"].values.sum())
    return np.array(SFR)/((length/h)**3)


def calculate_mass_per_bin(time_edges, interpolated_SFRD):
    """
    Calculates the mass created per bin given.

    Parameters
    ----------
    time_edges: float array
        An array containing the bin edges in Gyrs.
    interpolated_SFRD : scipy spline interpolation
        An intepolation of the SFRD in yrs.

    Returns
    -------
    numpy array
        A numpy array containing the mass generated in each bin in :math:`M_\odot`

    """
    mass_per_bin = np.zeros(len(time_edges)-1)
    for i in range(1, len(time_edges)):
        t1 = time_edges[i-1]
        t2 = time_edges[i]
        mass_per_bin[i-1] = interpolate.splint(t1*1e9, t2*1e9, interpolated_SFRD)
    return mass_per_bin



@numba.njit()
def calculate_event_rates(edges, mass_per_bin, DTD):
    """ Calculates the events rates by combining BPASS models and the
    cosmological simulations.

    Parameters
    ----------
    edges: numpy array
        An array containing the edges of the mass_per_bin binning. length: N
    mass_per_bin: numpy array
        An array containing the amout of mass per bin. length: N-1
    DTD : numpy array
        An array containing the Delay Time Distribution extracted from BPASS in
        normal BPASS binning.
    Returns
    -------
    numpy array
        A numpy array containing the event rates
        with units: :math:`\#events/yr/Gpc^3`.

    """
    event_rate = np.zeros(len(mass_per_bin))
    DTD_width = np.diff(BPASS_LINEAR_TIME_EDGES)
    for b, mass in enumerate(mass_per_bin):
        t = edges[b+1]
        for j in range(0, b+1):
            p1 = t - edges[j]
            p2 = t - edges[j+1]
            bin_events = kea.helpers._numba_integral(p2, p1, BPASS_LINEAR_TIME_EDGES, DTD, DTD_width)
            event_rate[j] += mass * bin_events # return Events/Mpc3 per bin

    return event_rate/np.diff(edges) * 1e9 # To units Events/yr/Gpc3


def calculate_2D_SFRD(data_file):
    """
    Extracts the Stellar Formation Rate Density from the Milenium Simulation
    and puts it into a 2D grid of time and metallicity (pandas DataFrame).

    Parameters
    ----------

    Returns
    -------
    pandas DataFrame
        A 2D Stellar Formation Rate Density over snapshot number and metallicity.
    """

    data = pd.read_csv(data_file,
                       comment="#",
                       usecols=["snapnum", "stellarMass",
                                "metalsStellarMass", "sfr"]
                       )

    np_SFRD = kea.helpers._numba_calculate_2D_SFRD(data["snapnum"].to_numpy(),
                                       data["stellarMass"].to_numpy(),
                                       data["metalsStellarMass"].to_numpy(),
                                       data["sfr"].to_numpy())
    return pd.DataFrame(np_SFRD.T, columns=BPASS_NUM_METALLICITIES)
