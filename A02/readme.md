Noel Kubach, Zhoutao Zhang

# Assignment 02

## Task 3

### Dokumentation 
```
Required packages:
  - NumPy
  - biopython
  
usage: 
kubach_ramachandran.py [-h] -i < PATH_TO_PDB_FILE_1 > < PATH_TO_PDB_FILE_2 > ... < PATH_TO_PDB_FILE_n >
-o < PATH_FOR_PDF_OUTPUT >

Compute the Ramachandran map for a given set of protrins.

options:
  -h, --help            show this help message and exit
  -i PATH, --input PATH
                        Path to pdb files with protein structure
  -o PATH, --output PATH
                        Path to save the generated Ramachandran maps.
```

### Usage example:
```bash
python kubach_ramachandran.py -i testfiles/test1.pdb testfiles/test2.pdb testfiles/test3.pdb -o output/
```
