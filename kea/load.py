"""
Functions to import data from BPASS

Author: Max Briel
"""
import gzip
import os
from hoki import load
import kea.hist
from kea.constants import *
import pandas as pd
import numpy as np
import pickle

def gunzip(source_filepath, dest_filepath, block_size=65536):
    """Unpacks a zipped file.

    Parameters
    ----------
    source_filepath : string
        The path to the file to unpack
    dest_filepath : string
        The filename to unpack to
    block_size : int
        The size of data to keep in memory at once

    """
    with gzip.open(source_filepath, 'rb') as s_file, \
            open(dest_filepath, 'wb') as d_file:
        while True:
            block = s_file.read(block_size)
            if not block:
                break
            else:
                d_file.write(block)
        d_file.write(block)


def packnload(file):
    """Load the data from a BPASS zipped file.

    Parameters
    ----------
    file : string
        A string pointing to a zipped version of BPASS data.

    Returns
    -------
    pandas DataFrame
        A pandas DataFrame containing the model data from the BPASS file,
        which is a *hoki* output.

    """
    gunzip(file+".gz", file)
    out = load.model_output(file)
    os.remove(file)
    return out


def loadBPASS(file, types):
    """load BPASS rates.

    Parameters
    ----------
    file : string
        A string pointing to a file to be loaded
    types : array of strings
        An array of types of events to load from the BPASS file

    Returns
    -------
    dict of BPASS histograms
        A dictionary with BPASS histograms of the supernovae event rates
        from the file with an entry of each give type in **types**.
        The event rates are in #events/yr/:math:`M_\odot`.

    """
    SNe_rates = packnload(file)
    rates = {i:kea.hist.BPASS_hist() for i in types}

    val = list(rates.values())[0]
    log_bins = val.getLogBins()
    for t in rates:
        if t == "ccsn":
            rates[t].Fill(log_bins, SNe_rates[["IIP", "II", "Ib", "Ic"]].sum(axis=1))
        else:
            rates[t].Fill(log_bins, SNe_rates[t].values)

    # normalise the rates to #Events/yr/M_sun
    bin_widths = np.array([val.getBinWidth(i)*1e9 for i in range(0, val.getNBins())])
    for i in rates:
        rates[i] = rates[i]/1e6/bin_widths

    return rates


def loadGW(file, types):
    """ Load the given **types** from the Graviational wave events file

    Parameters
    ----------
    file : string
        The path to a gravitational wave file from BPASS
    types : array of strings
        the type of events to extract from the file

    Returns
    -------
    dict of BPASS histograms
        A dictionary with BPASS histograms of the gravitational wave event
        rates from the file with an entry of each give type in **types**.
        The event rates are in #events/yr/:math:`M_\odot`.

    """
    data = pd.read_csv(file,
                    sep= "\s+",
                    names=["log_age", "NSNS", "BHNS", "BHBH", "age_yrs"],
                    engine="python")
    rates = {i:kea.hist.BPASS_hist() for i in types}
    val = list(rates.values())[0]
    log_bins = val.getLogBins()
    for i in rates:
        rates[i].Fill(log_bins, data[i])

    # normalise the rates to #Events/yr/M_sun
    bin_widths = np.array([val.getBinWidth(i)*1e9 for i in range(0, val.getNBins())])
    if "v2.2bray" in file:
        for i in rates:
            rates[i] = rates[i]/1e6/bin_widths
    elif "v2.1hobbs" in file:
        for i in rates:
            rates[i] = rates[i]/bin_widths
    elif "v2.2hobbs" in file:
        for i in rates:
            rates[i] = rates[i]/1e6/bin_widths
    else:
        print("NO Normalisation applied")
    return rates


def loadAllGW(folder):
    """Loads all the GW event metallicities in the given folder.

    Parameters
    ----------
    folder : str
        FOlder containing the GW events

    Returns
    -------
    dict
        A dictionary containing per metallicity the BPASS rates per event type.

    """
    compact_types = ["BHBH", "BHNS", "NSNS"]
    num_Z = BPASS_NUM_METALLICITIES

    rates = {}
    for x, i in enumerate(BPASS_METALLICITIES):
        rates[num_Z[x]] = loadGW(folder+"gwmergerdata.z"+i+".dat", compact_types)

    return rates


def loadAllRates(data_folder):
    """ Loads the SNe & compact merger rates for all available
    metallicities in BPASS.

    Loads v2.2 Hobbs for the Gravitational Wave events.

    Parameters
    ----------
    data_folder : string
        Folder containing the BPASS & GW models

    Returns
    -------
    dict["metallicity"]["event type"]
        A dictionary in a dictionary containing BPASS histograms. All BPASS
        metalicities are used and all BPASS & GW event types or used.

        The event rate are in #events/yr/:math:`M_\odot`.
    """
    SNe_types = ["IIP", "II", "Ib", "Ic", "Ia", "LGRB", "PISNe"]
    compact_types = ["BHBH", "BHNS", "NSNS"]
    num_Z = BPASS_NUM_METALLICITIES

    rates = {}
    for x, i in enumerate(BPASS_METALLICITIES):
        rates[num_Z[x]] = loadBPASS(data_folder+"bpass_v2.2.1_imf135_300/supernova-bin-imf135_300.z"+i+".dat", SNe_types)
        rates[num_Z[x]].update(loadGW(data_folder+"GWrates/v2.2hobbs/gwmergerdata.z"+i+".dat", compact_types))

    return rates


def calculate2dSFRD(data_file, time_file):
    """Extracts the Stellar Formation Rate Density from the Milenium Simulation
    and puts it into a 2D grid of time and metallicity (pandas DataFrame).

    Parameters
    ----------
    data_file : str
        The filename of the cosmological simulation.
    time_file : str
        The filename of the time relation in the cosmological simulation.

    Returns
    -------
    pandas DataFrame
        A 2D Stellar Formation Rate Density over snapshot number and metallicity.
    """
    # edges are defined by BPASS metallicity binning
    edges = [0.00005, 0.0005, 0.0015, 0.0025, 0.0035, 0.005, 0.007, 0.009, 0.012, 0.017, 0.025, 0.035]
    text_metals = ["em5", "em4", "001","002", "003", "004", "006", "008", "010", "014", "020", "030", "040"]
    num_metals = [0.00001, 0.0001, 0.001, 0.002, 0.003, 0.004, 0.006, 0.008, 0.010, 0.014, 0.020, 0.030, 0.040]

    data = pd.read_csv(data_file, comment="#")
    tr = pd.read_csv(time_file, comment="#")

    SFR = pd.DataFrame(0.0, columns=num_metals, index=tr["snapNum"])

    for index, row in data.iterrows():
        snapnum = row["snapnum"]
        stellarMass = row["stellarMass"]  #both in 10**10 solar masses
        if stellarMass == 0.0:
            continue
        metallicity = row["metalsStellarMass"]/stellarMass

        sfr = row["sfr"]
        location = np.where(metallicity <= edges)[0]
        if len(location) == 0:
            Zbin = num_metals[0]
        else:
            Zbin = num_metals[location[0]]
        SFR.at[int(snapnum), Zbin] += sfr

    SFR = SFR/((62.5/0.73)**3)  #Msun/yr/Mpc3
    return SFR


def load2dSFRD(file):
    """Loads the 2d Stellar Formation Rate Density from a pickle file.

    Parameters
    ----------
    file : str
        The filename of the pickled 2d SFRD

    Returns
    -------
    pandas DataFrame
        An unpickled pandas DataFrame containing the Stellar Formation Rate
        Density over snapshot number of the Millenium Simulation and metallicity
        of BPASS.
    """
    SFRS = pickle.load(open(file, 'rb'))
    return SFRS
