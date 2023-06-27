"""
Authors: Zhoutao Zhang, Noel Kubach
"""

import argparse
import os
from Bio.PDB import PDBParser
from Bio.PDB import Residue
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


# Setting up argument parser
args_parser = argparse.ArgumentParser(
    prog="DSSP",
    description="Calculates hydrogen bond energies from the first record"
                "and the first chain of a pdb file based on the DSSP algorithm.")
args_parser.add_argument(
    "-i", "--input",
    dest="input",
    help="Path to .pdb file",
    type=str,
    nargs=1,
    required=True
)
args_parser.add_argument(
    "-o", "--output",
    dest="output",
    help="Paths to output folder. Files dssp_matrix.tsv as well as two figures will be created automatically.",
    nargs=1,
    type=str,
    required=True
)
args = args_parser.parse_args()


def main():
    # reading in first chain of first model of PDB file:
    parser = PDBParser()
    chain = parser.get_structure(os.path.basename(args.input[0]), args.input[0]).get_models().__next__().get_chains().__next__()

    # generating array with residues of the chain:
    residues = [residue for residue in chain.get_residues()]

    # Initializing matrix for bond energies with zeroes
    E = np.array([[float(0) for i in range(len(chain))] for j in range(len(chain))])

    # Iterating over all pairs of residues and calculating the bond energy
    for i in range(len(residues)):
        for j in range(len(residues)):
            E[i][j] = dssp(residues[i], residues[j])

    # saving output .tsv
    write_tsv(E, args.output[0])

    # generating plots
    create_plots(E, args.output[0])


def dssp(acceptor: Residue, donor: Residue) -> float:
    """
    Calculates the hydrogen bond energy between two residues of a protein
    according to the DSSP algorithm.
    :param acceptor: acceptor residue (with =O group)
    :param donor: donor residue (with NH group)
    :return: float with bond energy in kj/mol
    """
    # extracting atoms and positions:
    at_donor = residue_to_dict(donor)
    at_acceptor = residue_to_dict(acceptor)

    # checking if all necessary atoms are present:
    if "N" not in at_donor.keys() or\
            "H" not in at_donor.keys() or\
            "C" not in at_acceptor.keys() or\
            "O" not in at_acceptor.keys():
        return np.nan

    # calculating distances:
    rON = dist(at_acceptor["O"], at_donor["N"])
    rCH = dist(at_acceptor["C"], at_donor["H"])
    rOH = dist(at_acceptor["O"], at_donor["H"])
    rCN = dist(at_acceptor["C"], at_donor["N"])

    # calculating energy:
    # 0.084: partial charges
    # 332: factor from paper
    # 4.184: kcal -> kj
    # result in kj/mol (paper)
    return 0.084*(1/rON + 1/rCH - 1/rOH - 1/rCN)*332*4.184


def dist(point1: np.array, point2: np.array) -> float:
    """
    Calculates the distance between two points in three-dimensional space.
    :param point1: numpy array with coordinates of first point
    :param point2: numpy array with coordinates of first point
    :return: float with distance between points
    """
    squared_sum = float(0)
    for i in range(len(point1)):
        squared_sum += np.power((point1[i] - point2[i]), 2)
    return np.sqrt(squared_sum)


def residue_to_dict(residue: Residue) -> dict:
    """
    Extracts all atoms from a Bio.PDB.Residue object and
    returns a dictionary with the atom ids and the coordinates.
    :param residue: Biopython PDB file residue
    :return: dictionary with atom ids -> coordinates
    """
    atoms = {}
    for atom in residue.get_atoms():
        atoms[atom.get_id()] = np.array(atom.get_coord())
    return atoms


def write_tsv(E: np.array, path: str) -> None:
    np.savetxt(os.path.join(path, "dssp_matrix.tsv"), E, delimiter="\t")


def create_plots(E: np.array, path: str) -> None:
    # parameters:
    ticks = [i + 0.5 for i in range(0, len(E))]
    tick_labels = range(1, len(E)+1)
    tick_label_size = 7
    titlesize = 20

    # first plot with continuous color scheme
    plt.figure()
    plt.figure(figsize=(10, 10))

    im = sns.heatmap(E, vmax=-2.1, cmap="binary_r", linewidths=0.5, cbar=True, square=True, linecolor="lightgrey",
                     cbar_kws={"shrink": 0.6})
    im.xaxis.tick_top()
    im.set_title("Hydrogen Bond Energy", pad=40)
    im.title.set_size(titlesize)
    colorbar = im.collections[0].colorbar
    colorbar.ax.tick_params(labelsize=15)
    colorbar.outline.set_linewidth(2)
    colorbar.outline.set_edgecolor("black")

    plt.xticks(ticks=ticks, labels=tick_labels, fontsize=tick_label_size)
    plt.yticks(ticks=ticks, labels=tick_labels, fontsize=tick_label_size)

    plt.tight_layout()
    plt.savefig(os.path.join(path, "dssp_continuous.png"), dpi=300)

    # second plot with discrete values for below/above threshold
    E1 = np.array([[1 if cell <= -2.1 else 0 for cell in row] for row in E])
    colors = ((1.0, 1.0, 1.0), (0.0, 0.0, 0.0))
    colormap = matplotlib.colors.LinearSegmentedColormap.from_list("custom", colors, len(colors))

    # creating plot
    plt.figure()
    plt.figure(figsize=(10, 10))

    im = sns.heatmap(E1, cmap=colormap, linewidths=0.5, cbar=True, square=True, center=0.5, linecolor="lightgrey",
                     cbar_kws={"shrink": 0.2})
    im.set_title("Hydrogen Bond Energy below -2.1 kj/mol", pad=40)
    im.title.set_size(titlesize)
    im.xaxis.tick_top()

    plt.xticks(ticks=ticks, labels=tick_labels, fontsize=tick_label_size)
    plt.yticks(ticks=ticks, labels=tick_labels, fontsize=tick_label_size)

    # setting colorbar ticks
    colorbar = im.collections[0].colorbar
    colorbar.set_ticks([0.25, 0.75])
    colorbar.set_ticklabels(["Above\nThreshold\n(-2.1 kj/mol)", "Hydrogen\nBonds"])
    colorbar.ax.tick_params(labelsize=15)
    colorbar.outline.set_linewidth(2)
    colorbar.outline.set_edgecolor("black")

    plt.tight_layout()
    plt.savefig(os.path.join(path, "dssp_discrete.png"), dpi=300)


main()
