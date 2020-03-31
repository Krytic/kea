#
# Functions to plot frequent figures
#
# Author: Max Briel
import scipy.interpolate
import matplotlib.pyplot as plt

def plotSFRD(SFRD, func, x):
    ynew = scipy.interpolate.splev(x*1e9, func, der=0)
    plt.step(x, SFRD, where="mid", label="MilliMil data")
    plt.plot(x, ynew, label="Spline interpolation")
    plt.xlabel("lookback Time (Gyr)")
    plt.legend()
    plt.ylabel("SFRD (M$_\odot$/yr/Mpc$^3$)")
    plt.ylim(0)
    plt.tight_layout()
    return None

def plotRates(histograms):
    for i in histograms:
        histograms[i].plotLog(label=i)

    plt.ylabel(r"Events/yr/M$_\odot$")
    plt.yscale("log")
    plt.legend()
    plt.tight_layout()
    return None

def plotEventRates_LB(histograms):
    for i in histograms:
        histograms[i].plot(label=i)
    plt.ylabel("Events/yr/Gpc$^3$")
    plt.xlabel("Lookback Time (Gyr)")
    plt.yscale("log")
    plt.xlim(0,15)
    plt.ylim(10, 2e6)
    plt.legend()
    plt.tight_layout()
    return None

def plotEventRates_Z(histograms, Z):
    for i in histograms:
        plt.plot(Z, histograms[i].getValues()[:len(Z)], label=i)

    plt.yscale("log")
    plt.ylabel("Events/yr/Gpc$^3$")
    plt.xlabel("Redshift")
    plt.ylim(1e1, 2e6)
    plt.xlim(0, 6)
    plt.legend()
    plt.tight_layout()
    return None
