#
# Functions to help extract the star formation rate & get an event rate
#
# Author: Max Briel
#
import numpy as np
from kea.hist import histogram, BPASS_hist
from scipy import interpolate, integrate
from kea.constants import *
from scipy.optimize import fminbound

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


def getEventRates(SFR, DTDs, sampling_rate, now):
    """ Calculates the events rates by combining BPASS models and the
    cosmological simulations.

    Parameters
    ----------
    SFR : scipy.interpolate spline function
        A scipy.interpolate spline of the stellar formation rates (in :math:`M/yr/Mpc^3`)
    DTDs : dictionary of BPASS_hists
        The Delay Time Distributions extracted in BPASS ordered in a histogram
        based on the event type
    sampling_rate : int
        The sampling rate for the new histogram (number of bins)
    now : float
        The current age of the universe in Gyrs.

    Returns
    -------
    dictionary of histograms
        A dictionary containing histograms with the event rates per event type
        with units :math:`\#events/yr/Gpc^3`.

    """
    events = {i: histogram(0, now, sampling_rate) for i in DTDs}
    item = list(events.values())[0]
    lookback = item.getBinEdges()
    Nbins = item.getNBins()
    mass = 0

    for i in range(1, Nbins+1):
        t1 = lookback[i-1]
        t2 = lookback[i]
        mass = interpolate.splint(t1*1e9, t2*1e9, SFR)  # Total Mass/Mpc3 in bin
        if mass == 0:
            continue
        for j in range(0, i):
            p1 = t2 - lookback[j]
            p2 = t2 - lookback[j+1]
            for d in DTDs:
                bin_events = DTDs[d].integral(p2, p1) # Event in the bin per M
                events[d].Fill(lookback[j], bin_events*mass)

    # normalise to on a per yr basis
    bins = np.array([item.getBinWidth(i)*1e9 for i in range(0, Nbins)])
    for i in events:
        events[i] = events[i]/bins
    return events


def getIndividualEventRates(SFR, Zfunc, DTDs, sampling_rate, now):

    events = {i: histogram(0, now, sampling_rate) for i in list(DTDs.values())[0]}
    item = list(events.values())[0]
    lookback = item.getBinEdges()
    Nbins = item.getNBins()

    for i in range(1, Nbins+1):
        t1 = lookback[i-1]
        t2 = lookback[i]

        # calculate the metallicity and find nearest
        Zvalues = interpolate.splev([t1, t2], Zfunc)
        Z = np.abs(Zvalues[0] - Zvalues[1])/2
        location = np.where(Z <= BPASS_METALLICITY_EDGES)[0]
        metallicity = 0
        if len(location) == 0:
            metallicity = BPASS_NUM_METALLICITIES[0]
        else:
            metallicity = BPASS_NUM_METALLICITIES[location[0]]

        mass = interpolate.splint(t1*1e9, t2*1e9, SFR)  # Total Mass/Mpc3 in bin

        if mass == 0:
            continue
        for j in range(0, i):
            p1 = t2 - lookback[j]
            p2 = t2 - lookback[j+1]
            M = DTDs[metallicity]
            for d in M:
                bin_events = M[d].integral(p2, p1)
                events[d].Fill(lookback[j], bin_events*mass)

    # normalise to on a per yr basis
    bins = np.array([item.getBinWidth(i)*1e9 for i in range(0, Nbins)])
    for i in events:
        events[i] = events[i]/bins
    return events


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
    print(calculateLB(zmin), LB)
    zbest, resval, ierr, ncall = fminbound(f,zmin, zmax, maxfun=maxfun, full_output=1, xtol=ztol)
    return zbest
