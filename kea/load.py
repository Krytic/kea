#
# Functions to import data from BPASS
#
# Author: Max Briel
import gzip
import os
from hoki import load
import kea.hist
import pandas as pd
import numpy as np

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
                    names=["log_age", "BHBH", "BHNS", "NSNS", "age_yrs"],
                    engine="python")
    rates = {i:kea.hist.BPASS_hist() for i in types}
    val = list(rates.values())[0]
    log_bins = val.getLogBins()
    for i in rates:
        rates[i].Fill(log_bins, data[i])

    # normalise the rates to #Events/yr/M_sun
    bin_widths = np.array([val.getBinWidth(i)*1e9 for i in range(0, val.getNBins())])
    for i in rates:
        rates[i] = rates[i]/1e6/bin_widths

    return rates


def loadAllRates(data_folder):
    """ Loads the SNe & compact merger rates for all available
    metallicities in BPASS

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
    metallicities = ["em5", "em4", "001","002", "003", "004", "006", "008", "010", "014", "020", "030", "040"]
    SNe_types = ["ccsn", "Ia", "LGRB", "PISNe"]
    compact_types = ["BHBH", "BHNS", "NSNS"]

    rates = {}
    for i in metallicities:
        rates[i] = loadBPASS(data_folder+"bpass_v2.2.1_imf135_300/supernova-bin-imf135_300.z"+i+".dat", SNe_types)
        rates[i].update(loadGW(data_folder+"GWrates/v2.2hobbs/gwmergerdata.z"+i+".dat", compact_types))

    return rates
