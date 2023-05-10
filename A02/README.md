Zhoutao Zhang, Noel Kubach

# Assignment 02
## Task 3 - Ramachandran Maps

### Required packages
In out implementation we used `argparse`, `Bio` (biopython), `numpy` and `mathplotlib`.
These packages need to be installed in order to run our program.

### Usage
```
usage: RamachandranPlot [-h] -i INPUT [INPUT ...] -o OUTPUT

Plots a Ramachandran diagram to a output file.

options:
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        Paths to .pdb files
  -o OUTPUT, --output OUTPUT
                        Paths to output. Ending specifies file type (e.g. .pdf).
```
### Example
```shell
python kubach_ramachandran.py -i resources/1igt.pdb -o resources/1IGT.pdf
python kubach_ramachandran.py -i resources/5ire.pdb -o resources/5IRE.pdf
python kubach_ramachandran.py -i resources/2mwk.pdb -o resources/2MWK.pdf
python kubach_ramachandran.py -i resources/5ire.pdb resources/5ire.pdb resources/2mwk.pdb -o resources/all.pdf
```
The resulting pdfs can be found in the [resources](resources) folder.\
![1IGT.png](resources%2F1IGT.png)\
![5IRE.png](resources%2F5IRE.png)\
![2MWK.png](resources%2F2MWK.png)\
![all.png](resources%2Fall.png)