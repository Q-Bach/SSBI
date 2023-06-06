Zhang Zhoutao, Noel Kubach

# Assignment 05
## Task 3
### Required packages
* Biopython
* Numpy
* Argparse

### Usage
```
usage: RMSD [-h] [-a] pdb1 pdb2

Aligns two protein structures and calculates the RMSD.

positional arguments:
  pdb1         Path to first .pdb file
  pdb2         Path to second .pdb file

options:
  -h, --help   show this help message and exit
  -a, --align  Whether or not to align the structures before calculation of the RMSD
```

### Example
RMSD between the structures 1M9K and 1M9M:
````shell
python Kubach_Noel_05_03.py resources/1m9k.pdb resources/1m9m.pdb
````
```
number of compared atoms: 6335
RMSD: 7.867950709899855
```

RMSD between the structures 1M9K and 1M9M with structural alignment:
````shell
python Kubach_Noel_05_03.py -a resources/1m9k.pdb resources/1m9m.pdb
````
```
number of compared atoms: 6335
RMSD: 0.6681580913016868
```