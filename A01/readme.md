Zhoutao Zhang, Noel Kubach

# Assignment 01
## Task 1

### Documentation:
```
usage: Nussinov RNA predictor [-h] -i PATH [--min-loop-length MIN_LOOP_LENGTH] [--score-GC GC_SCORE] [--score-AU AU_SCORE] [--score-GU GU_SCORE]

Predicts the secondary structure of RNA with the Nussinov algorithm.

options:
  -h, --help            show this help message and exit
  -i PATH, --input PATH
                        Path to fastA file with RNA sequence.
  --min-loop-length MIN_LOOP_LENGTH
                        Set the minimal allowed length for loops (default 3).
  --score-GC GC_SCORE   Provide a custom score for GC pairing (default 4).
  --score-AU AU_SCORE   Provide a custom score for AU pairing (default 3).
  --score-GU GU_SCORE   Provide a custom score for GU pairing (default 1).

```

### Usage example:
```bash
python kubach_nussinov_predictor.py -i testfiles/test.fasta --min-loop-length 5 --score-GC 10 --score-AU 5 --score-GU 2
```
### Output:
```
Filename: test.fasta
Sequence: RNA test sequence

Min-loop-length: 5
GC-score:        10        
AU-score:        5
GU-score:        2
Score:           57        

bpseq:
1 G 14
2 U 0
3 U 12
4 C 11
5 A 0
6 U 0
7 A 0
8 A 0
9 G 0
10 A 0
11 G 4
12 G 3
13 U 0
14 C 1
15 A 30
16 A 29
17 C 28
18 A 27
19 G 0
20 C 26
21 A 0
22 A 0
23 C 0
24 G 0
25 G 0
26 G 20
27 U 18
28 G 17
29 U 16
30 U 15

bracket-dot:
GUUCAUAAGAGGUCAACAGCAACGGGUGUU
(.((......)).)((((.(.....)))))

DP matrix:
     G  U  U  C  A  U  A  A  G  A  G  G  U  C  A  A  C  A  G  C  A  A  C  G  G  G  U  G  U  U
  G  0  0  0  0  0  0  0  5  5  7 10 12 14 22 22 22 30 30 30 35 35 35 45 45 45 45 45 47 52 57
  U  0  0  0  0  0  0  0  5  5  7 10 12 12 12 17 20 20 20 25 25 30 35 35 35 37 39 39 45 47 47
  U  0  0  0  0  0  0  0  0  2  5 10 12 12 12 15 15 15 15 25 25 30 30 30 35 37 37 37 42 45 47
  C  0  0  0  0  0  0  0  0  0  0 10 10 10 10 10 10 10 15 25 25 25 25 25 35 35 37 37 42 42 45
  A  0  0  0  0  0  0  0  0  0  0  0  2  7  7 10 10 10 15 15 20 25 25 25 25 27 27 32 35 40 45
  U  0  0  0  0  0  0  0  0  0  0  0  2  5  5 10 10 10 15 15 20 25 25 25 25 25 27 27 35 40 40
  A  0  0  0  0  0  0  0  0  0  0  0  0  5  5  5  5 10 10 10 20 20 20 22 22 25 25 27 30 35 40
  A  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 10 10 10 20 20 20 22 22 22 22 27 27 32 35
  G  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 10 10 10 20 20 20 22 22 22 22 27 27 30 35
  A  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 10 10 10 12 12 12 22 22 22 22 27 27 30 35
  G  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 10 10 10 12 12 12 22 22 22 22 25 25 30 35
  G  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  2 12 12 12 15 15 20 22 24 25 30 35
  U  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  2  2  5  5  5 10 20 22 22 25 30 35
  C  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 10 20 20 20 25 30 35
  A  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 10 10 10 15 25 30 35
  A  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 10 10 10 15 25 30 30
  C  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 10 10 10 15 25 25 25
  A  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 10 15 15 20 22
  G  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 10 12 15 17 17
  C  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 10 10 15 15 15
  A  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  5  5  5 10
  A  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  5  5
  C  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  2
  G  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  2
  G  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  G  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  U  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  G  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  U  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  U  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0

```