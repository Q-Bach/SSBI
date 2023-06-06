#!/usr/bin/env python
# Authors: Zhang Zhoutao, Noel Kubach

import Bio.PDB as PDB
from Bio import pairwise2
from Bio.PDB import Chain, Residue, Atom, Model
import numpy as np
import sys
import argparse


def parse_args():
    """
    Function to set up an argument parser for our script
    :return: parsed arguments
    """
    args_parser = argparse.ArgumentParser(
        prog="RMSD",
        description="Aligns two protein structures and calculates the RMSD.")
    args_parser.add_argument(
        "pdb1",
        help="Path to first .pdb file",
        nargs=1,
        type=str
    )
    args_parser.add_argument(
        "pdb2",
        help="Path to second .pdb file",
        nargs=1,
        type=str
    )
    args_parser.add_argument(
        "-a", "--align",
        help="Whether or not to align the structures before calculation of the RMSD",
        required=False,
        action="store_true"
    )
    args = args_parser.parse_args()
    return args


def main():
    # parsing args
    args = parse_args()

    # trying to read in PDBs:
    try:
        parser = PDB.PDBParser(QUIET=True)
        pdb1 = parser.get_structure("pdb1", args.pdb1[0]).get_models().__next__()
        pdb2 = parser.get_structure("pdb2", args.pdb2[0]).get_models().__next__()
    except:
        print("Error while reading in the .pdb files.", file=sys.stderr)
        exit(1)

    # extracting matching atoms
    atoms1, atoms2 = match_and_extract_atoms(pdb1, pdb2)

    # aligning structures
    if args.align:
        super_imposer = PDB.Superimposer()
        super_imposer.set_atoms(fixed=atoms1, moving=atoms2)
        super_imposer.apply(atoms2)

    # printing calculated values
    print("number of compared atoms: " + str(len(atoms1)))
    print("RMSD: " + str(rmsd(atoms1, atoms2)))


def chain_to_seq(chain: PDB.Chain) -> (str, [PDB.Residue]):
    """
    Extracts the amino acid sequence as one-letter code from a PDB chain object.
    :param chain: PDB chain object
    :return: string with sequence
    """
    # dictionary for conversion from three letter code to one letter code
    three_to_one = {
        'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
        'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    seq = ''
    residues = []
    for residue in chain.get_residues():
        three_letter = residue.resname
        if three_letter in three_to_one.keys():
            seq = seq + three_to_one[three_letter]
            residues.append(residue)
    return seq, residues


def match_and_extract_atoms(model1: PDB.Model, model2: PDB.Model) -> (list[PDB.Atom], list[PDB.Atom]):
    """
    Function to extract matching atoms from two PDB models.
    It performs sequence alignments for all pairs of chains
    and puts the atoms corresponding to matching residues in two lists.
    :param model1: first model
    :param model2: second model
    :return: two lists with matching atoms
    """
    # checking, whether both models have the same number of chains:
    if len(model1) != len(model2):
        print("The models have a different number of chains.", file=sys.stderr)
        exit(1)

    chains1 = [chain for chain in model1.get_chains()]
    chains2 = [chain for chain in model2.get_chains()]

    atoms1 = []
    atoms2 = []

    # iterating over the chains:
    for i in range(len(model1)):
        chain1 = chains1[i]
        chain2 = chains2[i]

        # extracting amino acid sequence and corresponding residues
        seq1, residues1 = chain_to_seq(chain1)
        seq2, residues2 = chain_to_seq(chain2)

        # aligning sequences:
        alignment = pairwise2.align.globalxx(seq1, seq2)

        # collecting atoms that match between both sequences
        # for the same amino acid all atoms are collected
        # for a match between different amino acids only CA is collected
        n = 0  # residue index in chain1
        m = 0  # residue index in chain2
        seq1 = alignment[0][0]
        seq2 = alignment[0][1]
        for j in range(len(seq1)):
            if seq1[j] == seq2[j]:  # match case

                # making dictionary from atoms of first residue
                atoms_dict = {}
                for atom in residues1[n]:
                    atoms_dict[atom.id] = atom
                # adding atom only if present in both residues
                for atom in residues2[m]:
                    if atom.id in atoms_dict.keys():
                        atoms1.append(atoms_dict[atom.id])
                        atoms2.append(atom)
                # increasing counter
                n += 1
                m += 1
            elif seq1[j] == "-":  # gap in seq1
                m += 1
            elif seq2[j] == "-":  # gap in seq2
                n += 1
            else:  # missmatch
                # increasing counter
                n += 1
                m += 1

    return atoms1, atoms2


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


def atom_distance(atom1: PDB.Atom, atom2: PDB.Atom) -> float:
    """
    Calculates the distance between two PDB.Atom objects
    :param atom1: first atom
    :param atom2: second atom
    :return: distance in Angstrom
    """
    pos1 = np.array(atom1.get_coord())
    pos2 = np.array(atom2.get_coord())
    return dist(pos1, pos2)


def rmsd(atoms1: list[PDB.Atom], atoms2: list[PDB.Atom]) -> float:
    """
    Calculates the RMSD given two lists of atoms.
    :param atoms1: list of atoms from first structure
    :param atoms2: list of atoms from second structure
    :return: RMSD (Angstrom)
    """
    # calculating squared pairwise distances
    distances = []
    for i in range(len(atoms1)):
        distances.append(atom_distance(atoms1[i], atoms2[i]) ** 2)

    # calculating rmsd
    return np.sqrt(np.sum(distances) / len(atoms1))


if __name__ == '__main__':
    main()
