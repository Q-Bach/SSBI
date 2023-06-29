Zhang Zhoutao, Noel Kubach

# Assignment 07

## Task 2 - _in silico_ digestion with trypsin
### Required packages
* Biopython
* Argparse

### Usage
````
usage: Kubach_Noel_07_02 [-h] fasta

Generates the peptides expected from a digestion with trypsinand calculates the masses of the peptides

positional arguments:
  fasta       Path to the fasta containing the protein sequence.

options:
  -h, --help  show this help message and exit
````

### Example
Digestion of the Alcohol dehydrogenase P07327:

````shell
python Kubach_Noel_07_02.py P07327.fasta
Found sequence of length: 375
Number of fragments after tryptic digestion: 38

Peptides                                           | monoisotopic mass (Da) | average mass (Da)     
---------------------------------------------------------------------------------------------
MVAVGICGTDDHVVSGTMVTPLPVILGHEAAGIVESVGEGVTTVKPGDK  | 4841.4652              | 4844.5893             
LDTMMASLLCCHEACGTSVIVGVPPDSQNLSMNPMLLLTGR          | 4317.0468              | 4320.1269             
KPIHHFLGISTFSQYTVVDENAVAK                          | 2800.4548              | 2802.1830             
VTPGSTCAVFGLGGVGLSAIMGCK                           | 2224.1054              | 2225.6652             
VCLIGCGFSTGYGSAVNVAK                               | 1944.9437              | 1946.2674             
NDVSNPQGTLQDGTSR                                   | 1687.7761              | 1688.7288             
EMTDGGVDFSFEVIGR                                   | 1757.7930              | 1758.9211             
FSLDALITHVLPFEK                                    | 1728.9450              | 1730.0375             
KPFSIEEVEVAPPK                                     | 1568.8449              | 1569.8190             
ELGATECINPQDYK                                     | 1579.7188              | 1580.7303             
INEGFDLLHSGK                                       | 1328.6724              | 1329.4759             
VIPLAIPQCGK                                        | 1137.6580              | 1138.4345             
NPESNYCLK                                          | 1066.4753              | 1067.1823             
IDAASPLEK                                          | 942.5022               | 943.0655              
AAVLWELK                                           | 928.5382               | 929.1278              
IIAVDINK                                           | 884.5331               | 885.0721              
KPIQEVLK                                           | 953.5909               | 954.1785              
GAILGGFK                                           | 761.4435               | 761.9200              
LVADFMAK                                           | 893.4680               | 894.0975              
MSTAGK                                             | 593.2843               | 593.6967              
AAGAAR                                             | 515.2816               | 515.5706              
AHEVR                                              | 610.3187               | 610.6715              
ECVPK                                              | 574.2785               | 574.6937              
TILMF                                              | 623.3352               | 623.8091              
FTCR                                               | 525.2369               | 525.6240              
VIK                                                | 358.2580               | 358.4821              
ICK                                                | 362.1988               | 362.4883              
FAK                                                | 364.2110               | 364.4455              
TWK                                                | 433.2325               | 433.5084              
SIR                                                | 374.2278               | 374.4411              
CK                                                 | 249.1147               | 249.3289              
IK                                                 | 259.1896               | 259.3495              
CR                                                 | 277.1209               | 277.3423              
DK                                                 | 261.1324               | 261.2787              
AK                                                 | 217.1426               | 217.2689              
SK                                                 | 233.1375               | 233.2683              
R                                                  | 174.1117               | 174.2035              
K                                                  | 146.1055               | 146.1901              

````

## Task 4

### Required packages
* Biopython
* Argparse

### Usage
````
usage: Zhang_Zhoutao_07_04.py [-h] -s <sequence>

Calculate the masses of all b- and y-ions and their m/z ratio

options:
  -s, --sequence <sequence> input sequence 
  -h, --help                show this help message and exit
````

### Example
````shell
python Zhang_Zhoutao_07_04.py -s HFEEDMGRK
Input sequence is: HFEEDMGRK
Length:  9
Seq | No. |     B     |    m/z
----------------------------------
H   |    1|  138.06673|  138.06673
F   |    2|  285.13514|  285.13514
E   |    3|  414.17773|  414.17773
E   |    4|  543.22032|  543.22032
D   |    5|  658.24726|  658.24726
M   |    6|  789.28775|  789.28775
G   |    7|  846.30921|  846.30921
R   |    8| 1002.41032| 1002.41032
K   |    9| 1130.50528| 1130.50528
Seq | No. |     Y     |    m/z
----------------------------------
H   |    9| 1148.51583| 1148.51583
F   |    8| 1011.45692| 1011.45692
E   |    7|  864.38851|  864.38851
E   |    6|  735.34592|  735.34592
D   |    5|  606.30333|  606.30333
M   |    4|  491.27639|  491.27639
G   |    3|   360.2359|   360.2359
R   |    2|  303.21444|  303.21444
K   |    1|  147.11333|  147.11333
````