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
from kea.constants import *
import numpy as np
import argparse
from scipy import interpolate
import pickle
import pandas as pd

def calculate_2D_event_rates(SFRD_file, time_file, data_folder, output_file, nr_bins):

    rates = kea.load.load_all_rates(data_folder)
    SFRD = pickle.load(open(SFRD_file, "rb"))
    tr = pd.read_csv(time_file, comment="#")
    edges = np.linspace(0,NOW, nr_bins+1)
    event_rates = pd.DataFrame(columns=rates.columns, index=edges[1:], dtype=np.float64)
    event_rates.index.name = "age_Gyr"
    for column in SFRD.columns:
        print(column)
        Z_SFR = SFRD[column]
        interpolated_SFRD = interpolate.splrep(np.flip(tr["lookbackTime"]*1e9), np.flip(Z_SFR.values), k=1)
        mass_per_bin = kea.rates.calculate_mass_per_bin(edges, interpolated_SFRD)

        for t in event_rates.columns.unique(level="Event Type"):
            event_rates[(t,column)] = kea.rates.calculate_event_rates(edges, mass_per_bin, rates[(t, column)].values)

    pickle.dump(event_rates, open(output_file, "wb"))
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

parser.add_argument("-N",
                    dest="nr_bins",
                    type=int,
                    required=True,
                    help="Number of bins used for final binning.")

parser.set_defaults(func=calculate_2D_event_rates)


if __name__ == "__main__":
    args = parser.parse_args()
    calculate_2D_event_rates(args.SFRD_file,
        args.time_file,
        args.data_folder,
        args.output_folder,
        args.nr_bins)
