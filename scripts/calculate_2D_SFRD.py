"""
A script to calculate the event rates of SNe & Compact mergers
while taking metallicities into account. Returns a Msun/yr/Mpc3

Example: python calculate2dSFRD.py
                                -i "../../data/large_galaxy.dat"
                                -o "../../data/2dSFRD.p"
Author: Max Briel
"""
import pickle
import kea.rates
import argparse

def calculate_2D_SFRD(data_file, output_file):
    print("Going through all the galaxies...", end='')
    SFRD = kea.rates.calculate_2D_SFRD(data_file)
    pickle.dump(SFRD, open(output_file, 'wb'))
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

parser.add_argument("-o",
                dest="output_folder",
                type=str,
                required=True,
                help="Output filename for the pickled SFR"
                )
parser.set_defaults(func=calculate_2D_SFRD)


if __name__ == "__main__":
    args = parser.parse_args()
    calculate_2D_SFRD(args.data_file, args.output_folder)
