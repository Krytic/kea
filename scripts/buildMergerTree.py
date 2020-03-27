"""
A script to build the merger tree history per galaxy for the inputted data.

Author: Max Briel

"""

import pandas as pd
import kea.mergerHistory as kea
import numpy as np
import pickle

import argparse

def buildMergerTree(data_folder, output_folder):
    cosmod = pd.read_csv(data_folder+"large_galaxy.dat", comment="#")
    root_nodes = kea.buildHistory(cosmod)
    pickle.dump(root_nodes, open(output_folder+"mergertrees.p", "wb"))

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
parser.set_defaults(func=buildMergerTree)

if __name__ == "__main__":
    args = parser.parse_args()
    buildMergerTree(args.data_folder, args.output_folder)
