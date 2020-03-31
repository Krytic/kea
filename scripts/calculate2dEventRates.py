"""
Script to calculate the 2D event rate (using pickles SFRD and BPASS)

Example:
python calculate2dEventRates.py -s "../../data/2dSFRD.p"
                                -t "../../data/timerel.dat"
                                -i "../../data/"
                                -o "../../data/2dEventRates.p"

"""


import kea.rates
import kea.load
import numpy as np
import argparse
from scipy import interpolate
import pickle
import pandas as pd

now = 13.799

def calculate2dEventRates(SFRD_file, time_file, data_folder, output_file):

    SFRD = kea.load.load2dSFRD(SFRD_file)
    rates = kea.load.loadAllRates(data_folder)
    tr = pd.read_csv(data_folder+"timerel.dat", comment="#")

    eventRates = {}
    for column in SFRD.columns:
        Z_SFR = SFRD[column]
        func = interpolate.splrep(np.flip(tr["lookbackTime"]*1e9), np.flip(Z_SFR.values), k=1)

        eventRates[column] = kea.rates.getEventRates(func, rates[column], 100, now)   #event/yr/Msun -> #events/yr/Mpc3

    # Normalise the event rates to be in #events/yr/Gpc$^3$
    for Z in eventRates:
        for i in eventRates[Z]:
            eventRates[Z][i] = eventRates[Z][i]*1e9

    pickle.dump(eventRates, open(output_file, "wb"))
    return None

parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-s",
                dest="SFRD_file",
                type=str,
                required=True,
                help="The file containing a 2D SFRD over time and metallicity"
                )
parser.add_argument("-t",
                    dest="time_file",
                    type=str,
                    required=True,
                    help="The file containing time information of the simulation")

parser.add_argument("-i",
                    dest="data_folder",
                    type=str,
                    required=True,
                    help="The folder containing the BPASS data")

parser.add_argument("-o",
                dest="output_folder",
                type=str,
                required=True,
                help="Output filename for the pickled Event Rates"
                )
parser.set_defaults(func=calculate2dEventRates)


if __name__ == "__main__":
    args = parser.parse_args()
    calculate2dEventRates(args.SFRD_file, args.time_file, args.data_folder, args.output_folder)
