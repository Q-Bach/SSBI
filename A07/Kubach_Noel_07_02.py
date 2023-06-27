from Bio import SeqIO
import argparse


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


def main():
    args = parse_args()

    with open(args.fasta[0], "r") as f:
        sequence = SeqIO.read(f, format="fasta")

    print("Found sequence of length: " + str(len(sequence)))

    # performing digestion
    peptides = trypsin(sequence)

    print("Number of fragments after tryptic digestion: " + str(len(peptides)))

    # calculating masses and printing report:
    # generating formatting string:
    longest_peptide = len(peptides[0])
    template = "{:<" + str(longest_peptide + 1) + "s} | {:<22.4f} | {:<22.4f}"
    print(str("{:<" + str(longest_peptide + 2) + "s} | {:<22s} | {:<22s}")
          .format("\nPeptides", "monoisotopic mass (Da)", "average mass (Da)"))  # Header
    print("-" * (longest_peptide + 44))
    for peptide in trypsin(str(sequence.seq)):
        print(template.format(peptide, monoisotopic_mass(peptide), average_mass(peptide)))


def trypsin(seq: str) -> []:
    """
    Function to predict the peptides resulting from a
    tryptic digestion of a protein.
    :param seq: input protein sequence
    :return: List of peptides
    """

    seq = seq.upper()
    peptides = []
    prev_cut = -1
    for i in range(len(seq)):
        aa = seq[i]
        aa_next = seq[i + 1] if i + 1 < len(seq) else ""
        if (aa == "R" or aa == "K") and aa_next != "P":
            peptides.append(seq[prev_cut + 1:i + 1])
            prev_cut = i
    if prev_cut + 1 != len(seq):
        peptides.append(seq[prev_cut + 1:])

    return sorted(peptides, key=lambda x: len(x), reverse=True)


def monoisotopic_mass(peptide: str) -> float:
    """
    Calculates the monoisotopic mass of a peptide.
    :param peptide: Sequence of the peptide
    :return: monoisotopic mass of the peptide
    """
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

    peptide = peptide.upper()
    mass = 15.99491 + 2 * 1.00782  # monoisotopic mass of water (C and N-terminus)
    # iterating over amino acids:
    for aa in peptide:
        if aa in aa_mass.keys():
            mass += aa_mass.get(aa)

    return mass


def average_mass(peptide: str) -> float:
    """
    Calculates the average mass of a peptide.
    :param peptide: Sequence of the peptide
    :return: average mass of the peptide
    """
    # Table of average masses from https://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html
    aa_mass = {
        "A": 71.0788,
        "R": 156.1875,
        "N": 114.1038,
        "D": 115.0886,
        "C": 103.1388,
        "E": 129.1155,
        "Q": 128.1307,
        "G": 57.0519,
        "H": 137.1411,
        "I": 113.1594,
        "L": 113.1594,
        "K": 128.1741,
        "M": 131.1926,
        "F": 147.1766,
        "P": 97.1167,
        "S": 87.0782,
        "T": 101.1051,
        "W": 186.2132,
        "Y": 163.1760,
        "V": 99.1326
    }

    peptide = peptide.upper()
    mass = 15.99977 + 2 * 1.00811  # average mass of water (C and N-terminus)
    # iterating over amino acids:
    for aa in peptide:
        if aa in aa_mass.keys():
            mass += aa_mass.get(aa)

    return mass


main()
