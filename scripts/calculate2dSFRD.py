"""
A script to calculate the event rates of SNe & Compact mergers
while taking metallicities into account. Returns a Msun/yr/Mpc3

Example: python calculate2dSFRD.py
                                -i "../../data/large_galaxy.dat"
                                -t "../../data/timerel.dat"
                                -o "../../data/2dSFRD.p"
Author: Max Briel
"""
import pandas as pd
import numpy as np
import pickle
import sys
import kea.load
import argparse

def calculate2dSFRD(data_file, time_file, output_file):
    print("Going through all the galaxies...", end='')
    SFR = kea.load.calculate2dSFRD(data_file, time_file)
    pickle.dump(SFR, open(output_file, 'wb'))
    print("DONE")
    return None

parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                dest="data_file",
                type=str,
                required=True,
                help="The file containing the galaxy data"
                )
parser.add_argument("-t",
                    dest="time_file",
                    type=str,
                    required=True,
                    help="The file containing time information of the simulation")

parser.add_argument("-o",
                dest="output_folder",
                type=str,
                required=True,
                help="Output filename for the pickled SFR"
                )
parser.set_defaults(func=calculate2dSFRD)


if __name__ == "__main__":
    args = parser.parse_args()
    calculate2dSFRD(args.data_file, args.time_file, args.output_folder)
