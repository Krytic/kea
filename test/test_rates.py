"""
Tests of the rates submodule

"""

import numpy as np
import pandas as pd
import kea.rates
import pickle
from kea.constants import *
from scipy import interpolate


data_folder = "../data/"


def test_getSFRD():

    h= 0.73
    cosmod = pd.read_csv(data_folder+"TEST_DATA/Millenium_Simulation.dat", comment="#")
    tr = pd.read_csv(data_folder+"cosmological_simulation/timerel.dat", comment="#")
    SFRD = kea.rates.getSFRD(cosmod, tr, 62.5, h)
    SFRD_check = np.loadtxt(data_folder+"TEST_DATA/SFRD_output.txt", delimiter=",")
    assert np.allclose(SFRD_check, SFRD)



def test_calculate_event_rates():
    nr_bins = 1000
    edges = np.linspace(0,NOW, nr_bins)

    cosmod = pd.read_csv(data_folder+"TEST_DATA/Millenium_Simulation.dat", comment="#")
    tr = pd.read_csv(data_folder+"cosmological_simulation/timerel.dat", comment="#")
    SNe_rates = kea.load.loadBPASS(data_folder+"bpass_v2.2.1_imf135_300/supernova-bin-imf135_300.z002.dat")

    SFRD = kea.rates.getSFRD(cosmod, tr, 62.5, h)
    SFRDfunc = interpolate.splrep(np.flip(tr["lookbackTime"]*1e9), np.flip(SFRD), k=1)

    ccsn_types = ['IIP', 'II', 'Ib', 'Ic']
    ccsn = np.zeros(51)
    for j in ccsn_types:
         ccsn += SNe_rates[j]
    SNe_rates["ccsn"] = ccsn

    mass_per_bin = kea.rates.calculate_mass_per_bin(edges, SFRDfunc)
    event_rate = kea.rates.calculate_event_rates(edges, mass_per_bin, SNe_rates['ccsn'].values)

    event_rate_check = pickle.load(open(data_folder+"TEST_DATA/event_rates_no_Z.p", "rb"))
    assert np.allclose(event_rate, event_rate_check)


def test_calculate_2D_SFRD():
    SFRD = kea.rates.calculate_2D_SFRD(data_folder+"TEST_DATA/Millenium_Simulation.dat")
    SFRD_check = pickle.load(open(data_folder+"TEST_DATA/2D_SFRD.p", "rb"))
    assert np.allclose(SFRD, SFRD_check)
