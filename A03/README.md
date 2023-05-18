Zhoutao Zhang, Noel Kubach

# Assignment 03

## Task 2

### Required packages
In out implementation we used `argparse`, `Bio` (biopython), `numpy` and `mathplotlib`.
These packages need to be installed in order to run our program.

### Usage
```
usage: zhang_assignment3 [-h] -i INPUT [INPUT ...]

Calculate how abundant are sheets, right-handed α-helices, and right-handed 310-helices in these proteins?

Calculate the propensities of the different amino acids to form these structures.

Plot histogram of the distance between backbone N atom and O atom.

options:
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        Paths to .pdb file(s)
```
### Example
```shell
python zhang_assignment3.py -i .\Supplementary
```
