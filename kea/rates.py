#
# Functions to help extract the star formation rate & get an event rate
#
# Author: Max Briel
#
import numpy as np
from kea.hist import histogram, BPASS_hist
from scipy import interpolate

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
