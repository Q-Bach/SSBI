"""
Authors: Zhoutao Zhang, Noel Kubach
"""

import argparse
import os
import sys
from Bio.PDB import PDBParser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Setting up argument parser
args_parser = argparse.ArgumentParser(
    prog="RamachandranPlot",
    description="Plots a Ramachandran diagram to a output file.")
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
    help="Paths to output PDF.",
    type=str,
    required=True
)
args = args_parser.parse_args()


def main():
    # reading in all pdb files
    pdbs = read_pdbs(args.input)

    # extracting file names for titles of figures:
    titles = [os.path.splitext(os.path.basename(name))[0] for name in args.input]
    # Iterating over pdb files and generting figures for each of them:
    figures = []
    for i in range(len(pdbs)):
        pdb = pdbs[i]
        # calculate phi and psi angles for whole structure
        phi_psi = calc_phi_psi(pdb)
        # generating ramachandran plot
        figures.append(plot_ramachandran(phi_psi, args.output, titles[i]))

    # saving figures to PDF:
    pdf = PdfPages(args.output)
    for fig in figures:
        fig.savefig(pdf, format='pdf')
    pdf.close()


def read_pdbs(paths: []):
    """
    Reads in one or more PDB files and returns a list of
    parsed objects.
    :param paths: list of paths to .pdb files
    :return: list of parsed pdb files
    """
    pdbs = []
    parser = PDBParser()
    for path in paths:
        try:
            pdbs.append(parser.get_structure(os.path.basename(path), path))
        except:
            print("Error while reading PDB file: " + path, file=sys.stderr)
            exit(0)
    return pdbs


def calc_phi_psi(pdb) -> list:
    """
    Iterates over chains and residues of a pdb object.
    For each enclosed residue, the phi and psi angle is calculated and stored
    in a two dimensional array.
    :param pdb: pdb object
    :return: [[phi, psi], [phi, psi], ...]
    """
    phi_psi = []
    # iterating over all models in one pdb:
    for model in pdb.get_models():
        # iterating over chains in one model
        for chain in model.get_chains():

            # Initializing previous two residues with None
            aa1 = None
            aa2 = None

            # iterating over all amino acids in one chain
            for aa0 in chain.get_residues():

                # ensures, that aa0, aa1 and aa2 have all required atoms
                # and that aa1 is enclosed by aa2 and aa0 (not for first and last)
                if aa1 and aa2\
                        and "CA" in aa2.child_dict.keys()\
                        and "N" in aa1.child_dict.keys()\
                        and "CA" in aa1.child_dict.keys()\
                        and "C" in aa1.child_dict.keys()\
                        and "N" in aa0.child_dict.keys():

                    # accessing all coordinates of the atoms required to calculate phi and psi
                    # for amino acid 'aa1'
                    CA_prev = aa2.child_dict["CA"].get_coord()
                    N = aa1.child_dict["N"].get_coord()
                    CA = aa1.child_dict["CA"].get_coord()
                    C = aa1.child_dict["C"].get_coord()
                    N_next = aa0.child_dict["N"].get_coord()

                    # calculating phi and psi for aa1
                    phi = dihedral(CA_prev, N, CA, C)
                    psi = dihedral(N, CA, C, N_next)
                    phi_psi.append([phi, psi])

                # updating previous two amino acids
                aa2 = aa1
                aa1 = aa0
    return phi_psi


def dihedral(a, b, c, d):
    """
    Computes the dihedral angle from four vectors.
    The computed angle corresponds to a rotation along the bond between b and c.
    The theoretical background was studied using: "https://leimao.github.io/blog/Dihedral-Angles/".
    :param a: numpy array with three coordinates
    :param b: numpy array with three coordinates
    :param c: numpy array with three coordinates
    :param d: numpy array with three coordinates
    :return: angle for rotation of bond b-c
    """

    # calculating vectors in planes
    ab = b - a
    bc = c - b
    cd = d - c

    # calculating normal vectors n1 and n2 of planes:
    n1 = np.cross(ab, bc)
    n2 = np.cross(bc, cd)

    # cos(thetha) for planes ab bc and bc cd
    cos_theta = (np.dot(n1, n2))/(np.linalg.norm(n1)*np.linalg.norm(n2))

    # sin(thetha) for the same planes
    sin_theta = np.dot(np.cross(n1, n2), bc)/(np.linalg.norm(n1)*np.linalg.norm(n2)*np.linalg.norm(bc))

    # calculating theta
    theta = np.arctan2(sin_theta, cos_theta)

    return theta*180/np.pi


def plot_ramachandran(phi_psi: list, path: str, title= list) -> plt.figure:
    """
    Generates a ramachandran plot with mathplotlib from the given phi psi angles
    :param title: Title of the diagram
    :param phi_psi: list with phi psi angles
    :param path: path (with ending) to save plot
    """
    fig = plt.figure()
    # generating scatter plot with phi and psi
    plt.scatter([phi[0] for phi in phi_psi], [psi[1] for psi in phi_psi], marker='.')

    # setting plot title and axis labels
    plt.title(title)
    plt.xlabel("phi")
    plt.ylabel("psi")

    # setting axis dimensions and ticks
    plt.xlim([-180, 180])
    plt.ylim([-180, 180])
    plt.xticks([-180, -90, 0, 90, 180])
    plt.yticks([-180, -90, 0, 90, 180])

    # making plot squared
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')

    # saving plot to defined location
    #plt.savefig(path, dpi=300)
    return fig

main()
