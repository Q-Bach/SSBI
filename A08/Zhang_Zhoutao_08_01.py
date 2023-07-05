import argparse

# Table of monoisotopic masses from https://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html
aa_mass_mono = {
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

def find_in_range(mass, error=0.055):
    """
    Function to find amino acids with certen mass within error range
    :return: list of possible amino acids
    """
    temp = []
    for aa in aa_mass_mono:
        if (aa_mass_mono[aa]-error <= mass <= aa_mass_mono[aa]+error):
            temp.append(aa)
    return temp 

def find_closest(mass):
    """
    Function to find an amino acid with  mass closest to given mass
    :return: one possible amino acid
    """
    temp = "x"
    min_delta = 65535
    for aa in aa_mass_mono:
        if(abs(aa_mass_mono[aa] - mass) < min_delta):
            temp = aa
            min_delta = abs(aa_mass_mono[aa] - mass)
    return temp

def decode(masses):
    """
    Function to decode b/y spectrum to amino acid sequences
    :return: list of pepide sequences
    """
    seq = [""]
    for mass in masses:
        aa = find_in_range(mass)
        if (len(aa)==0):
            aa = find_closest(mass)
            seq = [s + aa for s in seq]
        else:
            temp = seq.copy()
            seq = []
            for a in aa:
                seq.extend([s + a for s in temp])
    return seq

def reverse(seqs):
    """
    Function to reverse sequence generated by y spectrum
    :return: list of reversed sequence
    """
    if type(seqs)==str:
        rev_seq = ''
        for i in range(1,len(seqs)+1):
            rev_seq += seqs[-i]
        return rev_seq
    else:
        output = []
        for seq in seqs:
            output.append(reverse(seq))
        return output
    
def main():
    # set up an argument parser 
    args_parser = argparse.ArgumentParser(
        prog="Zhang_Zhoutao_08_02",
        description="Calculate peptid sequence based on b/y spectrum.")
    args_parser.add_argument(
        "-i",
        "--input",
        help="input b/y spectrum file path",
        nargs=1,
        type=str
    )
    args = args_parser.parse_args()

    # read spectrum
    path = args.input[0]
    spectrum1 = []
    spectrum2 = []
    with open("./materials/b_y_spectrum.txt", "r") as f:
        lines = f.readlines()
        for i in range(int(len(lines)/2)):
            spectrum1.append(float(lines[i*2].split()[0]))
            spectrum2.append(float(lines[i*2+1].split()[0]))
    
    # Determin which one is b spectrum
    min_1 = 65535
    min_2 = 65535
    aa_1 = spectrum1[0] - atom_mass["H"]
    aa_2 = spectrum2[0] - atom_mass["H"]
    for aa in aa_mass_mono:
        if(abs(aa_1-aa_mass_mono[aa])<min_1):
            min_1 = abs(aa_1-aa_mass_mono[aa])
        if(abs(aa_2-aa_mass_mono[aa])<min_2):
            min_2 = abs(aa_2-aa_mass_mono[aa])
    if(min_1<min_2):
        b_spectrum = spectrum1
        y_spectrum = spectrum2
    else:
        b_spectrum = spectrum2
        y_spectrum = spectrum1

    # Decode the sequence
    b_masses = [b_spectrum[0] - atom_mass["H"]]
    y_masses = [y_spectrum[0] - 3 * atom_mass["H"] - atom_mass["O"]]
    for i in range(len(b_spectrum)-1):
        b_masses.append(b_spectrum[i+1]-b_spectrum[i])
        y_masses.append(y_spectrum[i+1]-y_spectrum[i])
    b_seqs = decode(b_masses)
    y_seqs = reverse(decode(y_masses))

    # Search for the same sequences
    result = []
    for seq in b_seqs:
        if seq in y_seqs:
            result.append(seq)
    
    print("SSBI - Assignment 8")
    print("Noel Kubach, Zhoutao Zhang")
    print("---------------------------------------")
    print("b/y spectrum file: " + path)
    print("Sequence length: ", len(b_spectrum))
    print(len(result), "possible sequences found:")
    print("---------------------------------------")
    for i, seq in enumerate(result):
        print(seq+", ", end="")
        if ((i+1)%4==0):
            print()
    
if __name__ == "__main__":
    main()



    




