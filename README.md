# Barcode Generator

Barcode Generator is a Python script for generating and analyzing unique DNA barcodes. This script automates the process of designing DNA barcodes that are likely biologically inert and  meet specific criteria, including G/C content, homopolymer content, the minimum free energy of secondary structures and the level of similarity to other known sequences.

This script uses  multiple Python packages, sequence analysis tools, and parallel computing to process large numbers of sequences efficiently. The script's functionality can be divided into several steps, including generating sequences, analyzing secondary structure, conducting BLAST searches, and ranking the best barcodes.

## Dependencies

The script uses Python 3.8 and requires the following libraries:

- Biopython (BLAST)
- Pandas (data handling)
- NumPy (computations)
- RNAfold (secondary structure prediction)
- tqdm (progress bars)
- multiprocessing (parallel processing)
- itertools (generating combinations)

## Installation

`conda install -c conda-forge bioconda pandas numpy biopython`

`pip install tqdm`

`python -m pip install viennarna`


## Workflow

### Sequence Generation

The script begins by generating a specified number of random DNA sequences of a set length and with a G/C content between 40-60% and without  homopolymeric repeats >4. Only sequences likely to be biologically inert (i.e. lacking common restriction enzyme sites and common promoter and RBS consensus sequences) are accepted. It also saves these sequences to a fasta file.

### Secondary Structure Analysis

Next, the script uses RNAfold to analyze the secondary structure of each sequence, calculating the maximum minimum free energy (MFE) and average MFE for each sequence. Sequences with a max MFE greater than a specified threshold are removed. The remaining sequences are sorted by average MFE, and the top 10% are saved to a fasta file.

### BLAST Searches

The top sequences are then run through a BLAST search against the NCBI 'nt' database. The output of the BLAST search, which includes various alignment statistics, is saved to a csv file. This file is then cleaned to remove any misformatted  entries.

### Barcode Ranking

Next, the script finds the best barcodes based on minimal similarity to existing sequences, as determined by the  BLAST search results, and minimal predicted secondary structure. The results are filtered using a custom function, and the unique barcodes are selected. Information for each unique barcode is collected into a DataFrame, which is saved to a csv file.

### Barcode Pair Ranking

Finally, the script finds the best pairs of barcodes. Ranking each pair by various criteria including max MFE, average MFE, alignment length, percent identity, and query coverage.  The pairs are then sorted by weighted average rank and intra pair orthogonality, and the results are saved to a csv file.

## Running the Script

This script is intended to be run as a standalone Python script. The user may adjust the number of sequences to generate, the length of the sequences, and other parameters as needed. To execute the script, navigate to the directory containing the script and run:

`python generate_barcodes.py`

## Output

1. `full_barcode_set.fasta`: All the generated sequences.
2. `MFE_table.csv`: Secondary structure analysis results.
3. `2ndary_structure_optimized_barcodes.fasta`: The top 10% of sequences based on secondary structure.
4. `blast_alignments.out` and `blast_results.csv`: Results of the BLAST search.
5. `FINAL_barcodes.csv`: The best individual barcodes.
6. `ranked_barcode_pairs.csv`: The best pairs of barcodes.
