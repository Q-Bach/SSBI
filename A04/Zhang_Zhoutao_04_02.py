import numpy as np
import argparse

# The Sequence can be inputted as an argument
args_parser = argparse.ArgumentParser(
    prog = "Chou-Fasman Algorithm to predict protein secondary structure."
)
args_parser.add_argument(
    "-s", "--sequence",
    #dest = "sequence",
    help = "",
    type = str,
    default="MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE"
)
args = args_parser.parse_args()

## Scores for Chou-Fasman
# helix p_alpha relative probabilities
p_a = {
    'E': 1.53, 'A': 1.45, 'L': 1.34, 'H': 1.24, 'M': 1.20, 'Q': 1.17, 'W': 1.14, 'V': 1.14, 'F': 1.12,  # builders
    'K': 1.07, 'I': 1.00, 'D': 0.98, 'T': 0.82, 'S': 0.79, 'R': 0.79, 'C': 0.77,  # indifferent
    'N': 0.73, 'Y': 0.61, 'P': 0.59, 'G': 0.53  # breakers
}

# helix builder/indifferent/breaker class weights
w_a = {aa: (1 if score > 1.1 else -1 if score < 0.75 else 0.5) for aa, score in p_a.items()}

# strand p_beta relative probabilities
p_b = {
    'M': 1.67, 'V': 1.65, 'I': 1.60, 'C': 1.30, 'Y': 1.29, 'F': 1.28, 'Q': 1.23,  'L': 1.22, 'T': 1.20, 'W': 1.19,  # builders
    'A': 0.97, 'R': 0.90, 'G': 0.81, 'D': 0.80,  # indifferent
    'K': 0.74, 'S': 0.72, 'H': 0.71, 'N': 0.65, 'P': 0.62, 'E': 0.26  # breakers
}

# strand builder/indifferent/breaker class weights
# changed to 2/0/-1 encoding for convenience, so that a window of 5 is a strand core if sum(window) >= 5.
w_b = {aa: (2 if score > 1 else -1 if score < 0.78 else 0) for aa, score in p_b.items()}


# Find protential Helices in the sequence
def find_helix_(seq):
    seq_length = len(seq)
    # Generate an array containing score for every aa in the sequence
    weights = []
    probability = []
    helix = ["-" for _ in range(seq_length)]
    for i in range (0, seq_length):
        weights.append(w_a[seq[i]])
        probability.append(p_a[seq[i]])
    # Find the helix core and extend it
    for i in range(0, seq_length-6):
        weight = np.sum(weights[i:i+6])
        # Find the helix core
        if weight >= 4:
            # Extend left bound
            left_bound = i
            while True:
                if left_bound-1<0 or left_bound+3>seq_length:
                    break
                # If the sum of probability is greater than 4, extend the left bound
                elif  (np.sum(probability[left_bound-1:left_bound+3]) >= 4):
                    left_bound -= 1
                else:
                    break
            # Extend right bound
            right_bound = i+6
            while True:
                if right_bound-3<0 or right_bound+1>seq_length:
                    break
                # If the sum of probability is greater than 4, extend the right bound
                elif (np.sum(probability[right_bound-3:right_bound+1])>=4):
                    right_bound += 1
                else:
                    break
            # Fill the protential helix list
            for j in range(left_bound, right_bound):
                helix[j] = "H"
    return helix

def find_strands_(seq):
    seq_length = len(seq)
    # Generate an array containing score for every aa in the sequence
    weights = []
    probability = []
    # List to store protential sheets
    sheet = ['-' for _ in range(seq_length)]
    for i in range (0, seq_length):
        weights.append(w_b[seq[i]])
        probability.append(p_b[seq[i]])
    # Find the sheet core and extend it
    for i in range(0, seq_length-5):
        weight = np.sum(weights[i:i+5])
        if weight >= 5:
            # extend left bond
            left_bound = i
            while True:
                if left_bound-1<0 or left_bound+3>seq_length:
                    break
                # If the sum of probability is greater than 4, extend the left bound
                elif  (np.sum(probability[left_bound-1:left_bound+3]) >= 4):
                    left_bound -= 1
                else:
                    break
            right_bound = i+5
            while True:
                if right_bound-3<0 or right_bound+1>seq_length:
                    break
                # If the sum of probability is greater than 4, extend the right bound
                elif (np.sum(probability[right_bound-3:right_bound+1])>=4):
                    right_bound += 1
                else:
                    break
            for j in range(left_bound, right_bound):
                sheet[j] = 'S'
    return sheet

def conflict_handling(he, st, seq):
    pred = []
    i = 0
    while i < len(seq):
        if(he[i]=='H' and st[i]=='S'):
            end = i+1
            if end < len(seq):
                while(he[end]=='H' and st[end]=='S'):
                    end += 1
                    if end >= len(seq):
                        break
            p_h = 0
            p_s = 0
            for j in range(i,end):
                p_h += p_a[seq[j]]
                p_s += p_b[seq[j]]
            if p_h > p_s:
                for _ in range(i, end):
                    pred.append('H') 
            else:
                for _ in range(i, end):
                    pred.append('S')
            i = end
        elif(he[i]=='H' and st[i]!='S'):
            pred.append('H')
            i += 1
        elif(st[i]=='S'and he[i]!='H'):
            pred.append('S')
            i += 1
        else:
            pred.append('-')
            i += 1
    return ''.join(pred)

def main():
    seq = args.sequence
    helix = find_helix_(seq)
    sheet = find_strands_(seq)
    pred = conflict_handling(helix, sheet, seq)
    print("Sequence:         " + seq)
    print("Helices:          " + ''.join(helix))
    print("Sheets:           " + ''.join(sheet))
    print("Final prediction: " + pred)


if __name__ == '__main__':
    main()