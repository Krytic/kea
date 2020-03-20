#
# Progams to help extract the star formation rate & get an event rate
#
#
import numpy as np
from kea.hist import histogram, BPASS_hist
from scipy import interpolate

def getSFRD(cosmological_simulation, time_relations, length, h):
    """
    cosmological_simulation is cosmological simulation data with a column 'sfr'
    time_relations has the snapshot number and the lookback time together
    """
    SFR = []
    for i in time_relations["snapNum"]:
        data = cosmological_simulation.loc[cosmological_simulation["snapnum"] == i]
        if len(data) == 0:
            SFR.append(0)
        else:
            SFR.append(data["sfr"].values.sum())
    return np.array(SFR)/((length/h)**3)


def getEventRates(SFR, DTDs, sampling_rate, now):
    """
    SFR is a spline interpolation of the star formation rate (given in M/yr[/Mpc])
    DTDs are the Delay-Time-Distributions for each type of transient event
    you want to look at in dictionary form
    sampling_rate is the sampling rate at which the final distribution will be.
    now is the time for the here and now

    Outputs: a dictionary with histograms for each given DTD type in #events/yr/Gpc^3
    """
    events = {i: histogram(0, now, sampling_rate) for i in DTDs}
    item = list(events.values())[0]
    lookback = item.getBinEdges()
    Nbins = item.getNBins()
    mass = 0


    for i in range(1, Nbins+1):
        t1 = lookback[i-1]
        t2 = lookback[i]
        mass = interpolate.splint(t1*1e9, t2*1e9, SFR)  # Total Mass/Mpc in bin
        if mass == 0:
            continue
        for j in range(0, i):
            p1 = t2 - lookback[j]
            p2 = t2 - lookback[j+1]
            for i in DTDs:
                bin_events = DTDs[i].integral(p2, p1) # Event in the bin per M
                events[i].Fill(lookback[j], bin_events*mass)

    # normalise to on a per year basis
    bins = np.array([item.getBinWidth(i) for i in range(0, Nbins)])
    for i in events:
        events[i] = events[i]/bins
    return events
