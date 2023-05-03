"""
Authors: Zhang, Kubach
"""

import argparse
import os

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
    folding = nussinov(read_fasta(args.path))
    print('Filename: ' + os.path.basename(args.path))
    print('Min-loop: ' + str(args.min_loop_length))
    print('GC: ' + str(args.gc_score))
    print('AU: ' + str(args.au_score))
    print('GU: ' + str(args.gu_score))
    print('Score: ' + str(folding.score))
    print()
    print('bpseq:\n' + folding.as_bpseq())
    print('bracket-dot:' + folding.as_bracket_dot())

    # printing DP matrix:
    print('\nDP matrix:\n   ', end='')
    for base in folding.seq:
        print('{:>3.3s}'.format(base), end='')
    print()
    for row, base in zip(folding.matrix, folding.seq):
        print('{:>3.3s}'.format(base), end='')
        for value in row:
            print('{: 3d}'.format(value), end='')
        print()


class Folding:
    matrix = list(list())
    pairing = {}
    score = int
    seq = str

    def as_bpseq(self) -> str:
        bpseq = ""
        for base in self.pairing.keys():
            bpseq += '{:d} {:s} {:d}\n'.format(
                base,
                self.seq[base - 1],
                self.pairing.get(base)
            )
        return bpseq

    def as_bracket_dot(self) -> str:
        text = self.seq + '\n'
        for base in self.pairing.keys():
            if self.pairing.get(base) == 0:
                text += '.'
            else:
                if self.pairing.get(base) > base:
                    text += '('
                else:
                    text += ')'
        return text


class FastA:
    header = str
    seq = str

# Todo handle more than one sequence in fasta
def read_fasta(path: str) -> FastA:
    fasta = FastA()
    with open(path, 'r') as f:
        fasta.header = f.readline().removeprefix('>')
        fasta.seq = f.read().replace('\n', '')
    return fasta


def score(a: str, b: str) -> int:
    a = a.upper()
    b = b.upper()
    if (a == 'G' and b == 'C') or (a == 'C' and b == 'G'):
        return args.gc_score
    elif (a == 'A' and b == 'U') or (a == 'U' and b == 'A'):
        return args.au_score
    elif (a == 'G' and b == 'U') or (a == 'U' and b == 'G'):
        return args.gu_score
    else:
        return 0


def nussinov(fasta: FastA) -> Folding:
    seq = fasta.seq

    # Initializing matrix with zeroes
    m = [[0 for _ in range(len(seq))] for _ in range(len(seq))]

    # Loop to fill in matrix
    for x in range(args.min_loop_length + 1, len(seq)):  # diagonal index
        for i in range(0, len(seq) - x):  # row index
            j = i + x  # column index

            # calculating values of different cases
            case1 = m[i + 1][j]
            case2 = m[i][j - 1]
            case3 = m[i + 1][j - 1] + score(seq[i], seq[j])
            # case 4
            k = [x for x in range(i + 1, j - 1)]  # list for values of i<k<j
            if i + 2 >= j:
                case4 = 0
            else:
                row = [m[i][_k] for _k in k]  # row left of m[i][j]
                col = [m[_k + 1][j] for _k in k]  # col below m[i][j]
                row_col = [x + y for x, y in zip(row, col)]
                case4 = max(row_col)

            # assigning maximum value
            m[i][j] = max(case1, case2, case3, case4)

    # Printing matrix for debugging
    #for row in m:
    #    print(row)

    # creating new Folding object
    folding = Folding()
    folding.seq = seq
    folding.matrix = m
    folding.score = m[0][len(seq) - 1]

    # performing traceback
    folding = traceback(folding)

    return folding


def traceback(folding: Folding) -> Folding:
    m = folding.matrix
    seq = folding.seq
    pairing = {x: 0 for x in range(1, len(seq) + 1)}  # setting all connections to 0 (no connections)

    def recursion(i, j):
        if i < j:
            if m[i][j] == m[i + 1][j]:  # case 1
                recursion(i + 1, j)
            elif m[i][j] == m[i][j - 1]:  # case 2
                recursion(i, j - 1)
            elif m[i][j] == m[i + 1][j - 1] + score(seq[i], seq[j]):  # case 3, i and j are paired
                pairing[i + 1] = j + 1  # storing new pair
                pairing[j + 1] = i + 1
                recursion(i + 1, j - 1)
            else:  # case 4
                for k in range(i + 1, j - 1):
                    if m[i][j] == m[i][k] + m[k + 1][j]:
                        recursion(i, k)
                        recursion(k + 1, j)
                        break

    recursion(0, len(seq) - 1)
    folding.pairing = pairing
    return folding


main()