Noel Kubach, Zhoutao Zhang

# Assignment 04

## Task 1
### Required packages
The python script `Kubach_Noel_04_01` requires the packages
numpy, Biopython, mathplotlib and seaborn

### Usage
```
options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to .pdb file
  -o OUTPUT, --output OUTPUT
                        Paths to output folder. Files dssp_matrix.tsv as well as two figures will be created
                        automatically.
```

### Example
```shell
python Kubach_Noel_04_01.py -i resources/5jxv.pdb -o resources
```
This generates the files [`dssp_matrix.tsv`](resources/dssp_matrix.tsv),
[`dssp_continuous.png`](resources/dssp_continuous.png) and
[`dssp_discrete.png`](resources/dssp_discrete.png).

## Task 2

### Required packages
In out implementation we used `argparse` and `numpy`.
These packages need to be installed in order to run our program.

### Usage
```
usage: zhang_chou_fasman.py [-h] -s SEQUENCE

Predict the secondary structure of the input sequence

options:
  -h, --help            show this help message and exit
  -s, --sequence        input sequence can be specified here
                        if no sequence was given, the program will run on 5JXV
```
### Example
```shell
python zhang_chou_fasman.py -s MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE
```
