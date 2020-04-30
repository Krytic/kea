"""
Functions to import data from BPASS

Author: Max Briel
"""
import gzip
import os
import hoki.load
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
    out = hoki.load.model_output(file)
    os.remove(file)
    out.set_index('log_age', inplace=True)
    return out


def loadBPASS(file):
    """load BPASS rates.

    Parameters
    ----------
    file : string
        A string pointing to a file to be loaded

    Returns
    -------
    pandas DataFrame
        A pandas DataFrame of normalised supernovae event rates
        from the file with an entry of each give type in the
        The event rates are in #events/yr/:math:`M_\odot`.

    """
    SNe_rates = packnload(file)
    bins = SNe_rates['age_yrs']
    SNe_rates = SNe_rates.div(1e6*SNe_rates['age_yrs'], axis=0)
    SNe_rates.drop(['age_yrs',
                    'low_mass',
                    'e_Ia',
                    'e_IIP',
                    'e_II',
                    'e_Ib',
                    'e_Ic',
                    'e_LGRB',
                    'e_PISNe',
                    'e_low_mass'], axis=1, inplace=True)
    return SNe_rates


def loadGW(file):
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

                    names=["log_age", "NSNS", "BHNS", "BHBH", "age_yrs"],
                    engine="python")
    data.set_index('log_age', inplace=True)

    bins = data['age_yrs']
    names = {i:i for i in data.columns.values}

    # normalise the rates to #Events/yr/M_sun
    if "v2.2bray" in file:
        data = data.div(1e12*data['age_yrs'], axis=0)
        names["BHBH"] = "BHBHv22bray"
        names["BHNS"] = "BHNSv22bray"
        names["NSNS"] = "NSNSv22bray"
    elif "v2.1hobbs" in file:
        data = data.div(1e6*data['age_yrs'], axis=0)
        names["BHBH"] = "BHBHv21hobbs"
        names["BHNS"] = "BHNSv21hobbs"
        names["NSNS"] = "NSNSv21hobbs"
    elif "v2.2hobbs" in file:
        data = data.div(1e12*data['age_yrs'], axis=0)
        names["BHBH"] = "BHBHv22hobbs"
        names["BHNS"] = "BHNSv22hobbs"
        names["NSNS"] = "NSNSv22hobbs"
    else:
        print("NO Normalisation applied")
    data.drop('age_yrs', axis=1, inplace=True)
    data.rename(names, axis=1, inplace=True)
    return data


def loadAllGW(folder):
    """Loads all the GW event metallicities in the given folder.

    Parameters
    ----------
    folder : str
        Folder containing the GW events

    Returns
    -------
    dict
        A dictionary containing per metallicity the BPASS rates per event type.

    """
    rates = kea.load.loadGW(folder
                    + "gwmergerdata.z"
                    + BPASS_METALLICITIES[0]
                    + ".dat")

    columns = rates.columns.values
    arrays = [BPASS_NUM_METALLICITIES,columns]

    cols = pd.MultiIndex.from_product(arrays,
                                        names=['Metallicity', 'Event Type'])

    rates = pd.DataFrame(columns=np.linspace(6,11,51),
                         index=cols,
                         dtype=np.float64)

    for num, name in enumerate(BPASS_METALLICITIES):
        # Add GW events
        data = kea.load.loadGW(folder + "gwmergerdata.z"+name+".dat")
        rates.loc[(BPASS_NUM_METALLICITIES[num], slice(None))][columns[0]:columns[-1]] = data.T

    return rates.swaplevel(0,1).T


def load_all_GW(data_folder):
    v21_hobbs = data_folder+"GWrates/v2.1hobbs/"
    v22_hobbs = data_folder+"GWrates/v2.2hobbs/"
    v22_bray = data_folder+"GWrates/v2.2bray/"
    v21h_rates =kea.load.loadAllGW(v21_hobbs)
    v22b_rates = kea.load.loadAllGW(v22_bray)
    v22h_rates = kea.load.loadAllGW(v22_hobbs)
    return pd.concat([v21h_rates, v22b_rates, v22h_rates], axis=1, sort=False)


def load_all_rates(data_folder):
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
    arrays = [BPASS_NUM_METALLICITIES,
                BPASS_SUPERNOVA_TYPES
                +["NSNSv22hobbs", "BHNSv22hobbs", "BHBHv22hobbs"]]

    columns = pd.MultiIndex.from_product(arrays,
                                         names=['Metallicity', 'Event Type'])

    rates = pd.DataFrame(columns=np.linspace(6,11,51), index=columns, dtype=np.float64)
    rates.T.index.name = "log_age"

    for num, name in enumerate(BPASS_METALLICITIES):
        # Add SNe Events:
        data = kea.load.loadBPASS(data_folder+"bpass_v2.2.1_imf135_300/supernova-bin-imf135_300.z"+name+".dat")
        rates.loc[(BPASS_NUM_METALLICITIES[num], slice(None))]["IIP":"PISNe"] = data.T

        # Add GW events
        data = kea.load.loadGW(data_folder+"GWrates/v2.2hobbs/gwmergerdata.z"+name+".dat")
        rates.loc[(BPASS_NUM_METALLICITIES[num], slice(None))]["NSNSv22hobbs":"BHBHv22hobbs"] = data.T

    return rates.swaplevel(0,1).T


def load_all_colours(data_folder):
    """load all BPASS colout files

    Parameters
    ----------
    data_folder: string
        A string pointing to the BPASS folder to load


    Returns
    -------
    pandas DataFrame
        A pandas DataFrame of normalised colours for all metallicities
        The event rates are in absolute Magnitude (Vega) for a starburst at
        t=0 with :math:`10^6 M_\odot`.

    """
    x = packnload(f"{data_folder}/colours-bin-imf135_300.z002.dat")

    arrays = [BPASS_NUM_METALLICITIES, x.columns.values]

    columns = pd.MultiIndex.from_product(arrays,
                                         names=['Metallicity', 'colours'])

    colours = pd.DataFrame(columns=columns, index=np.linspace(6,11,51), dtype=np.float64)
    colours.index.name = "log_age"
    for num, name in enumerate(BPASS_METALLICITIES):
        data = packnload(f"{data_folder}/colours-bin-imf135_300.z{name}.dat")
        colours[BPASS_NUM_METALLICITIES[num]] = data.to_numpy()

    return colours
