Zhang Zhoutao, Noel Kubach

# Assignment 10

## Task 1 – Protein-Protein-Interaction Networks
### Required packages
* Argparse
* numpy
* tqdm
* pandas
* matplotlib
* powerlaw
* sys
* warnings

### Usage
````
usage: Zhang_Zhoutao_10_01.py [-h] -i <filepath>

options:
  -i, --input <file path>   input file path 
  -h, --help                show this help message and exit
````

### Example
````shell
python Zhang_Zhoutao_10_1.py -i human_interactome.csv
Reading file: 100%|████████████████████████████████████████████████████████| 342353/342353 [00:02<00:00, 152598.07it/s]
Number of nodes:  21557
Number of edges:  342353
Calculating best minimal value for power law fit
Searching for connected components: 100%|█████████████████████████████████████| 21557/21557 [00:00<00:00, 23334.30it/s]
Number of connected components:  26
The largest connected components contains  21521  nodes
Calculating distance: 100%|██████████████████████████████████████████████████████████| 100/100 [14:52<00:00,  8.93s/it]
````



