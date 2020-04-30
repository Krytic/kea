"""

Test for the loading submodule

"""
import kea.load
import pandas as pd
from numpy import allclose
def test_load_BPASS():
    """
    Test to check the load BPASS function.
    Want a pandas DataFrame as an output
    """
    data_folder = "../data/"
    SNe_rates = kea.load.loadBPASS(data_folder+"bpass_v2.2.1_imf135_300/supernova-bin-imf135_300.z020.dat")
    SNe_check = pd.read_csv(data_folder+"TEST_DATA/loadBPASS_output.txt", index_col=0)
    assert allclose(SNe_rates, SNe_check)


def test_load_GW():
    """
    Test to check the loadGW function.
    Want a pandas DataFrame as an output.
    """
    data_folder = "../data/"
    GW_rates = kea.load.loadGW(data_folder+"GWrates/v2.2hobbs/gwmergerdata.z020.dat")
    GW_check = pd.read_csv(data_folder+"TEST_DATA/loadGW_output.txt", index_col=0)
    assert allclose(GW_rates, GW_check)
