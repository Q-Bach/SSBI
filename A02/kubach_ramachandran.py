"""
Authors: Zhoutao Zhang, Noel Kubach
"""

import argparse
import os
import sys
from Bio.PDB import PDBParser
import numpy as np

args_parser = argparse.ArgumentParser(
    prog="RamachandranPlot",
    description="Plots a Ramachandran diagram to a PDF file.")
args_parser.add_argument(
    "-i", "--input",
    dest="input",
    help="Paths to .pdb files",
    type=str,
    nargs='+',
    required=True
)
args_parser.add_argument(
    "-o", "--output",
    dest="output",
    help="Paths to output .pdf",
    type=str,
    required=True
)
args = args_parser.parse_args()


def main():
    # reading in all pdb files
    pdbs = []
    parser = PDBParser()
    for path in args.input:
        try:
            pdbs.append(parser.get_structure(os.path.basename(path), path))
        except:
            print("Error while reading PDB file: " + path, file=sys.stderr)
            exit(0)

def dihedral(a, b, c, d) -> (float, float):


print(args.input + args.output)
