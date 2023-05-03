"""
Authors: Zhoutao Zhang, Noel Kubach
"""
import argparse
import os

# Configuring parser with argparse to read in arguments properly and provide a help page.
parser = argparse.ArgumentParser(
    prog="Nussinov RNA predictor",
    description="Predicts the secondary structure of RNA with the Nussinov algorithm.")
parser.add_argument(
    '-i', '--input',
    dest='path',
    help='Path to fastA file with RNA sequence.',
    type=str,
    required=True)
parser.add_argument(
    '--min-loop-length',
    help='Set the minimal allowed length for loops (default 3).',
    dest='min_loop_length',
    type=int,
    default=3,
    required=False)
parser.add_argument(
    '--score-GC',
    help='Provide a custom score for GC pairing (default 4).',
    dest='gc_score',
    type=int,
    default=4,
    required=False)
parser.add_argument(
    '--score-AU',
    help='Provide a custom score for AU pairing (default 3).',
    dest='au_score',
    type=int,
    default=3,
    required=False)
parser.add_argument(
    '--score-GU',
    help='Provide a custom score for GU pairing (default 1).',
    dest='gu_score',
    type=int,
    default=1,
    required=False)
args = parser.parse_args()


def main():
    """
    Main program loop, which initializes Nussinov Algorithm and prints results to the console.
    """
    # reading in fastA file
    fastas = read_fasta(args.path)

    for fasta in fastas:

        # running nussinov algorithm
        folding = nussinov(fasta)

        # printing run summary
        print('\nFilename: ' + os.path.basename(args.path))     # input values
        print('Sequence: ' + fasta.header)                      #
        print()                                                 #
        print('Min-loop-length: ' + str(args.min_loop_length))  #
        print('GC-score:        ' + str(args.gc_score))         #
        print('AU-score:        ' + str(args.au_score))         #
        print('GU-score:        ' + str(args.gu_score))         #
        print('Score:           ' + str(folding.score))         #
        print()
        print('bpseq:\n' + folding.as_bpseq())                  # bpseq format
        print('bracket-dot:\n' + folding.as_bracket_dot())      # bracket-dot format

        # printing DP matrix:
        print('\nDP matrix:\n   ', end='')
        for base in folding.seq:
            print('{:>3.3s}'.format(base), end='')              # First row with RNA sequence
        print()
        for row, base in zip(folding.matrix, folding.seq):
            print('{:>3.3s}'.format(base), end='')              # First column with RNA sequence
            for value in row:
                print('{: 3d}'.format(value), end='')           # Rows with scores
            print()


class Folding:
    """
    Class for folding object, to store the original RNA sequence, the DP matrix,
    the archived score and a dictionary with all base pairings.
    """

    seq = str               # String with RNA sequence
    matrix = [[]]           # DP matrix
    score = int             # folding score
    pairing = {}            # base pairing (1-based base indices)

    def as_bpseq(self) -> str:
        """
        Method to format the base pairing in bpseq format.
        :return: string with base pairing in bpseq format
        """
        bpseq = ""
        for base in self.pairing.keys():
            bpseq += '{:d} {:s} {:d}\n'.format( # Template
                base,                           # first index
                self.seq[base - 1],             # base of first index
                self.pairing.get(base)          # second index
            )
        return bpseq

    def as_bracket_dot(self) -> str:
        """
        Method to format the base pairing in bracket-dot format.
        :return: string with base pairing in bracket-dot format
        """
        text = self.seq + '\n'                      # first row with RNA-sequence
        for base in self.pairing.keys():
            if self.pairing.get(base) == 0:
                text += '.'                         # "." in case of the base not being paired
            else:
                if self.pairing.get(base) > base:
                    text += '('                     # "(" if base is paired with base with larger index
                else:
                    text += ')'                     # ")" if base is paired with base with larger index
        return text


class FastA:
    """
    Class to store a read in fastA file.
    Stores the sequence and the header of a fastA.
    """
    header = str    # header of fastA
    seq = str       # sequence of fastA file

    def __init__(self, header: str, seq: str):
        self.header = header
        self.seq = seq


def read_fasta(path: str) -> [FastA]:
    """
    Function to read in a fastA file.
    :param path: path to the file
    :return: list with all found sequences as Fasta objects
    """
    fastas = []                         # list to store Fasta objects
    with open(path, 'r') as f:
        text = f.read()                 # read entire file
    entries = text.split('>')           # split file by '>' character (assumes no '>' within header line)
    for entry in entries:
        if entry:                       # for non-empty entries
            lines = entry.split('\n')   # split lines of fastA entry
            # Adding new Fasta object to list. First line becomes header, other lines become sequence
            fastas.append(FastA(lines[0], ''.join(lines[1:])))
    return fastas


def nussinov(fasta: FastA) -> Folding:
    """
    Nussinov algorithm to compute an ideal folding of an RNA sequence
    :param fasta: input Fasta object with RNA sequence to fold
    :return: Folding object with the calculated folding
    """
    def score(a: str, b: str) -> int:
        """
        Scoring function for the Nussinov algorithm
        :param a: first base
        :param b: second base
        :return: score of the two bases being paired
        """
        # Turning input to uppercase
        a = a.upper()
        b = b.upper()

        # returning correct score
        if (a == 'G' and b == 'C') or (a == 'C' and b == 'G'):
            return args.gc_score
        elif (a == 'A' and b == 'U') or (a == 'U' and b == 'A'):
            return args.au_score
        elif (a == 'G' and b == 'U') or (a == 'U' and b == 'G'):
            return args.gu_score
        else:
            return 0

    def traceback(folding: Folding) -> Folding:
        """
        Function to initialize the recursive traceback through the DP matrix
        :param folding: Folding object with filled in DP matrix
        :return: Folding object with added base pairing
        """

        # Extracting important variables
        m = folding.matrix
        seq = folding.seq
        # Initializing dictionary for base pairing. Setting all connections to 0 (no connections)
        pairing = {x: 0 for x in range(1, len(seq) + 1)}

        def recursion(i, j):
            """
            Recursive functino to find a traceback in a filled in DP matrix.
            Adds the found base pairs to the pairing dictionary.
            The code was taken from the lecture (2D, slide 41).
            :param i: row index
            :param j: column index
            """
            if i < j:
                if m[i][j] == m[i + 1][j]:                                  # case 1
                    recursion(i + 1, j)
                elif m[i][j] == m[i][j - 1]:                                # case 2
                    recursion(i, j - 1)
                elif m[i][j] == m[i + 1][j - 1] + score(seq[i], seq[j]):    # case 3, i and j are paired
                    pairing[i + 1] = j + 1                                  # storing new pair
                    pairing[j + 1] = i + 1
                    recursion(i + 1, j - 1)
                else:                                                       # case 4
                    for k in range(i + 1, j - 1):
                        if m[i][j] == m[i][k] + m[k + 1][j]:
                            recursion(i, k)
                            recursion(k + 1, j)
                            break

        # Initializing traceback in upper right corner of DP matrix
        recursion(0, len(seq) - 1)
        folding.pairing = pairing
        return folding

    # extracting important variable
    seq = fasta.seq

    # Initializing matrix with zeroes
    m = [[0 for _ in range(len(seq))] for _ in range(len(seq))]

    # Loop to fill in matrix
    for x in range(args.min_loop_length + 1, len(seq)):  # diagonal index (starting index specifies min-loop-length)
        for i in range(0, len(seq) - x):  # row index
            j = i + x  # column index

            # calculating values of different cases
            case1 = m[i + 1][j]
            case2 = m[i][j - 1]
            case3 = m[i + 1][j - 1] + score(seq[i], seq[j])
            # case 4
            k = [x for x in range(i + 1, j - 1)]    # list for values of i<k<j
            if i + 2 >= j:
                case4 = 0
            else:
                row = [m[i][_k] for _k in k]        # row left of m[i][j]
                col = [m[_k + 1][j] for _k in k]    # col below m[i][j]
                row_col = [x + y for x, y in zip(row, col)]
                case4 = max(row_col)

            # assigning maximum value of all cases
            m[i][j] = max(case1, case2, case3, case4)

    # creating new Folding object
    folding = Folding()
    folding.seq = seq
    folding.matrix = m
    folding.score = m[0][len(seq) - 1]

    # performing traceback
    folding = traceback(folding)

    return folding


main()
