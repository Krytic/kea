"""
A script to plot the rates of SNe & Compact mergers without taking metallicity
into account.


Several different plots are outputted:

* The Stellar Formation Rate

  1. Actual rate
  2. Interpolated rate

* The BPASS event rates per solar mass
* The Event rates per :math:`Gpc^3`

Author: Max Briel
"""

import pandas as pd
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
from kea.load import packnload
import kea.hist as kea
import kea.rates as kea

import argparse

def ratesNoMetallicity(data_folder, out_folder):

    h = 0.73
    now = 13.799

    # Read the cosmological simulation and the lookbacktime relation + GW events
    cosmod = pd.read_csv(data_folder+"cosmological_model.dat", comment="#")
    tr = pd.read_csv(data_folder+"timerel.dat", comment="#")
    data = pd.read_csv( data_folder+"GWrates/v2.2hobbs/gwmergerdata.z002.dat",
                        sep= "\s+",
                        names=["log_age", "BHBH", "BHNS", "NSNS", "age_yrs"],
                        engine="python")

    # Extract the star formation rate density
    SFRD = kea.getSFRD(cosmod, tr, 62.5, h)
    SFRDfunc = interpolate.splrep(np.flip(tr["lookbackTime"]*1e9), np.flip(SFRD), k=1)
    ynew = interpolate.splev(tr["lookbackTime"]*1e9, SFRDfunc, der=0)

    # plot the SFRD and its linear interpolation
    fig = plt.figure()
    plt.step(tr["lookbackTime"].values, SFRD, where="mid", label="MilliMil data")
    plt.plot(tr["lookbackTime"], ynew, label="Linear interpolation")
    plt.xlabel("lookback Time (Gyr)")
    plt.legend()
    plt.ylabel("SFRD (M$_\odot$/yr /Mpc$^3$)")
    plt.ylim(0)
    plt.tight_layout()
    plt.savefig(out_folder+"SFRD_MilliMil.pdf")


    # Get the different SNe rates
    SNe_rates = packnload(data_folder+"/bpass_v2.2.1_imf135_300/supernova-bin-imf135_300.z002.dat")

    # store the rates into histograms
    #types = ["IIP", "II", "Ib", "Ic", "Type_Ia", "LGRB", "PISNe", "BHBH", "BHNS", "NSNS"]
    types = ["ccsn", "Type_Ia", "LGRB", "PISNe", "BHBH", "BHNS", "NSNS"]
    rates = {i:kea.BPASS_hist() for i in types}

    log_bins = rates["ccsn"].getLogBins()
    rates["ccsn"].Fill(log_bins, SNe_rates[["IIP", "II", "Ib", "Ic"]].sum(axis=1))
    #rates["IIP"].Fill(log_bins, SNe_rates["IIP"].values)
    #rates["II"].Fill(log_bins, SNe_rates["II"].values)
    #rates["Ib"].Fill(log_bins, SNe_rates["Ib"].values)
    #rates["Ic"].Fill(log_bins, SNe_rates["Ic"].values)
    rates["Type_Ia"].Fill(log_bins,SNe_rates["Ia"].values)
    rates["LGRB"].Fill(log_bins,SNe_rates["LGRB"].values)
    rates["PISNe"].Fill(log_bins, SNe_rates["PISNe"].values)

    rates["BHBH"].Fill(log_bins, data["BHBH"])
    rates["BHNS"].Fill(log_bins, data["BHNS"])
    rates["NSNS"].Fill(log_bins, data["NSNS"])

    # normalise the rates to be in a Events/M/yr
    bin_widths = np.array([rates["ccsn"].getBinWidth(i)*1e9 for i in range(0, rates["ccsn"].getNBins())])
    for i in rates:
        rates[i] = rates[i]/1e6/bin_widths


    # plot the event rates per solar mass
    fig2 = plt.figure()
    for i in rates:
        rates[i].plotLog(label=i)

    plt.ylabel(r"Events/M$_\odot$/yr")
    plt.yscale("log")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_folder+"SNe_rates.pdf")


    s_rate = 1000
    events = kea.getEventRates(SFRDfunc, rates, s_rate, now) # return in Events/yr/Mpc
    fig3 = plt.figure()

    # Get the event rate against lookbacktime
    for i in events:
        events[i] = events[i]*1e9 #Get it it units Events/yr/Gpc
        events[i].plot(label=i)
    plt.ylabel("Events/yr/Gpc$^3$")
    plt.xlabel("Lookback Time (Gyr)")
    plt.yscale("log")
    plt.xlim(0,15)
    plt.ylim(10, 2e6)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_folder+"Event_rates_LB.pdf")

    # Get the event rate against redshift
    LBZrel = interpolate.splrep(np.flip(tr["lookbackTime"]), np.flip(tr["Z"]), k=2)
    end = 920
    Z = interpolate.splev(events["ccsn"].getBinEdges()[:end],LBZrel)


    fig4 = plt.figure()

    for i in events:
        plt.plot(Z, events[i].getValues()[:end], label=i)

    plt.yscale("log")
    plt.ylabel("Events/yr/Gpc$^3$")
    plt.xlabel("Redshift")
    plt.ylim(1e1, 2e6)
    plt.xlim(0, 6)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_folder+"Event_rates_Z.pdf")


parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                dest="data_folder",
                type=str,
                nargs=1,
                required=True,
                help="Folder containing the data files"
                )

parser.add_argument("-o",
                dest="output_folder",
                type=str,
                nargs=1,
                required=True,
                help="Output folder for plots"
                )
parser.set_defaults(func=ratesNoMetallicity)


if __name__ == "__main__":
    args = parser.parse_args()
    ratesNoMetallicity(args.data_folder, args.output_folder)
