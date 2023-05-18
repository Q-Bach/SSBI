"""
SSBI - Assignment 3
Authors: Zhoutao Zhang, Noel Kubach
"""

import argparse
from Bio.PDB import PDBParser
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Bio.PDBPaser will give a warning each time a chain ends.
# To avoid this, simply ignore all warnings
import warnings
warnings.filterwarnings("ignore")

# Read in path to pdb file or files
args_parser = argparse.ArgumentParser(
    prog = "Investigating secondary structure elements",
    description = "Calculate the abundance of secondary structures, the number of different amino acids in some secondary structures and the distance of N and O in helices"
)
args_parser.add_argument(
    "-i", "--input",
    dest = "input",
    help = "Path to .pdb files",
    type = str,
    required = True
)
args = args_parser.parse_args()


def abundance_ss(pdb_files):
    '''
    Function to solve assignment 2a
    How abundant are sheets, right-handed α-helices, and right-handed 310-helices in these proteins?
    ''' 
    # Initialize counter 
    Helix_index = [1,5]
    HELIX = {
        1 : "Right-handed alpha      ",
        5 : "Right-handed 310        "
    } # Two types of helices we interessted in
    HELIX_length = {
        1 : 0,
        5 : 0
    }
    SHEET_length = 0
    total_length = 0
    # Loop through pdb files
    for pdb_file in pdb_files:
        with open(pdb_file, "r") as f:
            lines = f.readlines();
            # all helix length of a pdb file
            for line in lines:
                # In pdb files, helices entry begins with HELIX
                if line[:5]=="HELIX":
                    # length is at col 72-76
                    length = int(line[71:76])
                    # type is at col 39-40
                    helix_type = int(line[38:40])
                    if helix_type in Helix_index:
                        HELIX_length[helix_type] += length
                    total_length += length       
            # all sheet length of a pdb file
            for line in lines:
                # In pdb files, sheets entry begins with SHEET
                if line[:5]=="SHEET":
                    # There's no length entry like in helix in sheet of pdb file
                    # But there are begin and end indeces, So we calcilate the length 
                    # using end-begin+1
                    length = int(line[33:37])-int(line[22:26])+1
                    SHEET_length += length
                    total_length += length
    # Calculate the relativ length
    for key in HELIX_length:
        HELIX_length[key] /= total_length
    SHEET_length /= total_length
    # Print the result
    print("Secondary structure     Abundance")
    for key in HELIX_length:
        print(HELIX[key], HELIX_length[key])
    print("SHEET                   ", SHEET_length)

def propensities_amino_acid(pdb_files):
    '''
    Function to solve assignment 2b
    What are the propensities of the different amino acids to form these structures?
    '''
    # We only interessed in 20 amino acids, ignore other residues like NEG and HOH. 
    amino_acids = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    # Initialize counter
    HELIX_aa_counter = {}
    SHEET_aa_counter = {}
    # Add entries for each amino acid in the counter
    for i in range(len(amino_acids)):
        HELIX_aa_counter[amino_acids[i]] = 0
        SHEET_aa_counter[amino_acids[i]] = 0
    # Loop through all files
    for pdb_file in pdb_files:
        pdbparser = PDBParser()
        # Read chains from pdb files
        structure = pdbparser.get_structure(os.path.splitext(os.path.basename(pdb_file)), pdb_file)
        chains = {chain.id: [residue.get_resname() for residue in chain.get_residues()] for chain in structure.get_chains()}
        with open(pdb_file, "r") as f:
            lines = f.readlines();
            for line in lines:
                if (line[:5]=="HELIX"):
                    if int(line[38:40])==1:
                        # Fine the locations of helices on the chain
                        chain = line[19] 
                        start = int(line[21:25])-1
                        end = int(line[33:37])
                        # count amino acids
                        for aa in chains[chain][start:end]:
                            if aa in amino_acids:
                                HELIX_aa_counter[aa] += 1
                if (line[:5]=="SHEET"):
                    # Find the locations of sheets on the chain
                    chain = line[21] 
                    start = int(line[22:26])-1
                    end = int(line[33:37])
                    # count amino acids
                    for aa in chains[chain][start:end]:
                        if aa in amino_acids:
                            SHEET_aa_counter[aa] += 1
    # Calculate the relative frequency
    total_aa_helix = 0
    for aa in amino_acids:
        total_aa_helix += HELIX_aa_counter[aa]
    for aa in amino_acids:
        HELIX_aa_counter[aa] /= total_aa_helix 
    total_aa_sheet = 0
    for aa in amino_acids:
        total_aa_sheet += SHEET_aa_counter[aa]
    for aa in amino_acids:
        SHEET_aa_counter[aa] /= total_aa_sheet
    print("Helices       Sheets") 
    for aa in amino_acids:
        print(aa+" & ", end="")
        print("{:.3f}".format(HELIX_aa_counter[aa]), end=" & ")
        print(aa+" & ", end="")
        print("{:.3f}".format(SHEET_aa_counter[aa]), end="\\\\ \n")
    
def distances(pdb_files):
    '''
    Function to solve assignment 2c
    For all proteins, compute the distances between backbone N atoms and the O atom
    of backbone C = O groups four residues earlier. Plot all these distances in a 
    histogram. Next, do the same for right-handed alpha-helices only, and plot them in a
    second histogram.
    '''
    # We use list to save distances
    distances = []
    distance_helices = []
    # loop through all files
    for pdb_file in pdb_files:
        pdbparser = PDBParser()
        structure = pdbparser.get_structure(os.path.splitext(os.path.basename(pdb_file)), pdb_file)
        chains = {chain.id:chain for chain in structure.get_chains()}
        Helices = []
        # Record the location of helices in a pdb file in  a string
        with open(pdb_file, "r") as f:
            lines = f.readlines();
            for line in lines:
                if (line[:5]=="HELIX"):
                    if int(line[38:40])==1:
                        chain = line[19] 
                        start = int(line[21:25])
                        end = int(line[33:37])
                        Helices.append((chain, start, end))
        # Find the coordinate of back bone N and O atom.
        # We have to use two dictionary to save the coordinates of
        # N and O atoms.
        # In o_coord_all, the keys are chain ids, in each entry, there 
        # will be a list saving all coordinates. It's easier for the 
        # calculation of distance of all amino acids. 
        # In o_coord, the keys are chain ids, in each entry, there 
        # will be a new dictionary to save the coordinates of N and O using 
        # residue ids as index, because the residue ids are not continuous and
        # do not always begins with 1.
        o_coord = {}
        o_coord_all = {}
        n_coord = {}
        n_coord_all = {}
        for chain_id in chains:
            n_coord_all[chain_id] = []
            n_coord[chain_id] = {}
            o_coord_all[chain_id] = []
            o_coord[chain_id] = {}
            for residue in chains[chain_id].get_residues():
                # Get the coordinate of  N and O from each residue
                atom_counter = 0
                n_index = -1
                o_index = -1
                for atom in residue.get_atoms():
                    if (atom.get_name()=="N" and n_index==-1):
                        # It seems that N atom does not always appear at
                        # the begining as the first atom of a residue 
                        # so we have to look for it 
                        n_index = atom_counter
                        n_coord[chain_id][residue.id[1]] = atom.get_coord()
                        n_coord_all[chain_id].append(atom.get_coord())
                    if (atom.get_name()=="O" and o_index==-1):
                        # O is also the same case here
                        o_index = atom_counter
                        o_coord[chain_id][residue.id[1]] = atom.get_coord()
                        o_coord_all[chain_id].append(atom.get_coord())
                # And there are residues without N and O, we ignore them 
                if n_index == -1:
                    n_coord[chain_id][residue.id[1]] = []
                    n_coord_all[chain_id].append([])
                if o_index == -1:
                    o_coord[chain_id][residue.id[1]] = []
                    o_coord_all[chain_id].append([])
        # Calculate the distances for all amino acids
        for chain_id in chains:
            chain_length = len(o_coord_all[chain_id])
            for i in range(chain_length):
                # Avoid index out of range
                if (i+4)>=chain_length:
                    break
                else:
                    o = o_coord_all[chain_id][i]
                    n = n_coord_all[chain_id][i+4]
                    # If the residue contains no o and n, we ignore it
                    if o==[] or n==[]:
                        continue
                    square_distance = 0
                    for j in range(len(o)):
                        square_distance += (n[j]-o[j])**2
                    distances.append(np.sqrt(square_distance))
        # Calculate the distance for helices
        for (chain_id, start, end) in Helices:
            # Avoid index out of range
            if start+4 >= end:
                continue
            else:
                try:
                    # To avoid missing index (like 1w2p chain A id 346 is missing), 
                    # amino acid other than 20 amino acids we know and other problems
                    # just ignor the errors
                    for i in range(start, end-4):
                        o = o_coord[chain_id][i]
                        n = n_coord[chain_id][i+4]
                        if o==[] or n==[]:
                            continue
                        square_distance = 0
                        for j in range(len(o)):
                            square_distance += (n[j]-o[j])**2
                        distance_helices.append(np.sqrt(square_distance))
                except:
                    continue
    # Plot the histogram
    plt.figure()
    plt.hist(distance_helices, bins=[i for i in np.arange(2,16,1)],density=True, label="Distance of only helices amino acids.", color='salmon', alpha=0.7)
    plt.hist(distances,bins=[i for i in np.arange(2,16,1)],density=True, label="Distance of all backbone amino acids.", color='palegreen', alpha=0.7)
    plt.xlabel("Distance (Angström) between backbone N atoms and the O atom \n of backbone C = O groups four residues earlier")
    plt.ylabel("Probability density")
    plt.legend()
    plt.show()


def main():
    # Get path of input file
    path = args.input
    # If it is a file or a folder
    if path[-4:]==".pdb":
        pdb_files = [path]
    else:
        pdb_files = []
        # Go through the folder to get all pdb files
        for file in os.listdir(path):
            # We only need pdb files, ignore other format
            if file.endswith(".pdb"):
                pdb_files.append(os.path.join(path, file))
    abundance_ss(pdb_files)
    propensities_amino_acid(pdb_files)
    distances(pdb_files)

if __name__ == '__main__':
    main()
        