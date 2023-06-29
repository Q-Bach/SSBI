from Bio import SeqIO
import argparse

# Table of monoisotopic masses from https://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html
aa_mass = {
    "A": 71.03711,
    "R": 156.10111,
    "N": 114.04293,
    "D": 115.02694,
    "C": 103.00919,
    "E": 129.04259,
    "Q": 128.05858,
    "G": 57.02146,
    "H": 137.05891,
    "I": 113.08406,
    "L": 113.08406,
    "K": 128.09496,
    "M": 131.04049,
    "F": 147.06841,
    "P": 97.05276,
    "S": 87.03203,
    "T": 101.04768,
    "W": 186.07931,
    "Y": 163.06333,
    "V": 99.06841
}

atom_mass = {
    "H": 1.00782,
    "O": 15.99491
}

def parse_args():
    """
    Function to set up an argument parser for our script
    :return: parsed arguments
    """
    args_parser = argparse.ArgumentParser(
        prog="Kubach_Noel_07_02",
        description="Generates the peptides expected from a digestion with trypsin"
                    "and calculates the masses of the peptides")
    args_parser.add_argument(
        "fasta",
        help="Path to the fasta containing the protein sequence.",
        nargs=1,
        type=str
    )
    args = args_parser.parse_args()
    return args

def b_y_mass(seq):
    """
    Function to calculate the mass of all possible b and y ions
    :return: masses of b-ions and y-ions
    """
    # Calculate the mass of all possible b and y ions
    b_masses = [] 
    y_masses = []
    for i in range(len(seq)):
        b_mass = 0
        y_mass = 0
        for j in seq[:i+1]:
            b_mass += aa_mass[j]
        for k in seq[-i-1:]:
            y_mass += aa_mass[k]
        b_mass += atom_mass["H"] # The weight in the table does not contain H+ and OH-, for b ions, we add H+ to it
        y_mass += 3 * atom_mass["H"] + atom_mass["O"] # Add 2 H+ and 1 OH-
        b_masses.append(b_mass)
        y_masses.append(y_mass)
    return b_masses, y_masses


def main():
    args = parse_args()

    with open(args.fasta[0], "r") as f:
        seq = SeqIO.read(f, format="fasta")

    print("Found sequence of length: " + str(len(seq)))

    b_mass, y_mass = b_y_mass(seq)

    # print the results
    print("Seq | No. |     B     |    m/z    ")
    print("----------------------------------")
    for i in range(len(seq)):
        print(f"{seq[i]:4}|{i+1:5}|{b_mass[i]:11.9}|{b_mass[i]/1:11.9}")
    
    print("Seq | No. |     Y     |    m/z    ")
    print("----------------------------------")
    for i in range(len(seq)):
        print(f"{seq[i]:4}|{len(seq)-i:5}|{y_mass[-i-1]:11.9}|{y_mass[-i-1]/1:11.9}")
    
if __name__ == '__main__':
    main()