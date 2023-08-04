#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 10:50:46 2023

@author: petevoor
"""

import random
from Bio.SeqUtils import MeltingTemp as mt
import RNA
import pandas as pd
from tqdm import tqdm
import os
import subprocess
import re
from itertools import combinations

#==============================================================================
# INPUTS
#==============================================================================

# Set the size of the possible barcode population
num_sequences = 1000
# Set barcode sequence length
length = 300
# Define the window size for the secondary structure search
window_size = 50

alignment_length_filter = 50
query_coverage_filter = 30
percent_identity_filter = 95

#==============================================================================
# COMMON VARIABLES
#==============================================================================

# Get the directory where this script is located
script_directory = os.path.dirname(os.path.abspath(__file__))

common_enzymes = {
    "BbsI-HF": "GAAGAC",
    "BfuAI": "ACCTGC",
    "BsaI-HFv2": "GGTCTC",
    "BsmBI-v2": "CGTCTC",
    "PaqCI": "CACCTGC",
    "SapI": "GCTCTTC",
    "EcoRI": "GAATTC",
    "HindIII": "AAGCTT",
    "XhoI": "CTCGAG",
    "BamHI": "GGATCC",
    "KpnI": "GGTACC",
    "SalI": "GTCGAC",
    "SmaI": "CCCGGG",
    "XbaI": "TCTAGA",
    "NotI": "GCGGCCGC",
    "NcoI": "CCATGG",
    "NheI": "GCTAGC",
    "SacI": "GAGCTC",
    "AvaI": "CYCGRG",
    "AflII": "CTTAAG",
    "BglII": "AGATCT",
    "ScaI": "AGTACT",
    "MluI": "ACGCGT",
}

# Define the consensus sequences for the different promoter elements
regulatory_elements = {
    'sigma70': {'-10': 'TATAAT', '-35': 'TTGACA'},
    'RpoS': {'-10': 'CTATACT', '-35': 'TTGACT'},
    'RpoN': {'-12': 'TGGCACG', '-24': 'GG.{4}TGC'},
    'RpoH': {'-10': 'CCTAT[AT]TATA[TA]T', '-35': 'CTTGAA'},
    'Shine-Dalgarno': {'RBS': 'AGGAGG'} # Added the Shine-Dalgarno sequence
}

#==============================================================================
# FUNCTIONS
#==============================================================================

def has_homopolymer(seq, length=5):
    for i in range(len(seq) - length):
        if seq[i:i+length] == length * seq[i]:
            return True
    return False

def melting_temp_in_range(seq, window=20, min_temp=50, max_temp=70):
    for i in range(len(seq) - window + 1):
        sub_seq = seq[i:i+window]
        temp = mt.Tm_Wallace(sub_seq)
        if temp < min_temp or temp > max_temp:
            return False
    return True

def find_regulatory_elements(sequence,regulatory_elements):
    for element_type, elements in regulatory_elements.items():
        for element, consensus in elements.items():
            if re.search(consensus, sequence): # Use regular expression search
                return False

    return True

def generate_random_sequence(length, gc_content_min=40, gc_content_max=60):
    while True:
        sequence = ''.join(random.choice('ATGC') for _ in range(length))
        gc_content = (sequence.count('G') + sequence.count('C')) / length * 100
        repeats = has_homopolymer(sequence, length=5)
        found_enzymes = {key: value for key, value in common_enzymes.items() if value in sequence}
        # Check if the conditions are met
        if not found_enzymes and not find_regulatory_elements(sequence, regulatory_elements) \
                and gc_content_min <= gc_content <= gc_content_max \
                and repeats == False \
                and melting_temp_in_range(sequence):
            return sequence

def generate_sequences(num_sequences, length):
    return [generate_random_sequence(length) for _ in tqdm(range(num_sequences), desc="Generating Sequences")]

def find_secondary_structure(sequence, window_size):
    mfe_values = []
    max_mfe = float('inf') # Positive infinity ensures that any real number will be lower
    # Check all 50 bp windows of the sequence
    for i in range(len(sequence) - window_size + 1):
        subsequence = sequence[i:i+window_size]
        (structure, mfe) = RNA.fold(subsequence)
        mfe_values.append(mfe)
        if mfe < max_mfe:
            max_mfe = mfe
    return max_mfe, sum(mfe_values) / len(mfe_values)

def clean_blast_results(blast_output_path):
    # Define the correct column names
    columns = [
        "Query Accession Version", "Percent Identity", "Alignment Length", 
        "E-value", "Bit Score", "Query Coverage", "Subject Length",
        "Subject Accession Version", "Subject Title", "Location", "Type"
    ]

    # Read the BLAST results
    blast_results = pd.read_csv(blast_output_path, sep=',', header=None, names=columns)

    # Find rows where 'Percent Identity' cannot be converted to a number
    non_numeric_rows = pd.to_numeric(blast_results['Percent Identity'], errors='coerce').isna()

    # Get Query Accessions associated with non-numeric rows
    accessions_to_remove = blast_results.loc[non_numeric_rows, 'Query Accession Version'].unique()

    # Remove rows associated with these Query Accessions
    filtered_blast_results = blast_results.loc[~blast_results['Query Accession Version'].isin(accessions_to_remove)]

    # Save the filtered DataFrame back to the original file
    filtered_blast_results.to_csv(blast_output_path, sep=',', index=False)

    
def hamming_distance(seq1, seq2):
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

# Define the BLAST results filter conditions
def filter_func(group):
    if len(group) <= 3:
        if all(group["Alignment Length"] < alignment_length_filter) and all(group["Query Coverage"] < query_coverage_filter) and all(group["Percent Identity"] < percent_identity_filter):
            return True
    return False


#==============================================================================
# GENERATE SEQUENCES
#==============================================================================


sequences = generate_sequences(num_sequences, length)

# Print the original sequences
file_path = os.path.join(script_directory, "full_barcode_set.fasta")
print(f"Writing the full set of generated barcodes to {file_path}")
with open(file_path, "w") as f:
    for sequence in sequences:
        query_number = sequences.index(sequence) + 1
        f.write(f">barcode_{query_number}\n{sequence}\n")
        
#==============================================================================
# ANALYZE SECONDARY STRUCTURE
#==============================================================================

# Lists to hold the data for the DataFrame
sequence_names = []
whole_sequences = []
max_mfe_values = []
avg_mfe_values = []

# Analyzing secondary structure for each generated sequence
for idx, sequence in tqdm(enumerate(sequences), total=len(sequences), desc="Analyzing Secondary Structure"):
    max_mfe, avg_mfe = find_secondary_structure(sequence, window_size)
    sequence_names.append(f"barcode_{idx + 1}")
    whole_sequences.append(sequence)
    max_mfe_values.append(max_mfe)
    avg_mfe_values.append(avg_mfe)

# Create the DataFrame
df = pd.DataFrame({
    'sequence_name': sequence_names,
    'whole_sequence': whole_sequences,
    'max_MFE': max_mfe_values,
    'avg_MFE': avg_mfe_values
})

# Remove rows with max MFE greater than -15
df = df[df['max_MFE'] >= -15]

# Sort the DataFrame by average MFE from least negative to most negative
df = df.sort_values(by='avg_MFE', ascending=False)

# Write the secondary structure screening results to a csv file.
csv_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "MFE_table.csv")
print(f"Writing results of secondary structure screening to {csv_file_path}")
df.to_csv(csv_file_path, index=False)

# Create a file path for the top 10% of sequences
top_10_file_path = os.path.join(script_directory, "2ndary_structure_optimized_barcodes.fasta")

top_n = int(len(df) * 0.1) # Calculates the top 10%
top_sequences = df.nlargest(top_n, 'avg_MFE')
top_10_names = top_sequences['sequence_name'].tolist() # Use top_sequences DataFrame
top_10_sequences = top_sequences['whole_sequence'].tolist()

print(f"Writing barcodes with minimal secondary structure to {top_10_file_path}")
with open(top_10_file_path, "w") as f:
    for name, sequence in zip(top_10_names, top_10_sequences):
        f.write(f">{name}\n{sequence}\n")
        
#==============================================================================
# BLAST SEARCH
#==============================================================================
        
# BLAST command line
blast_command = [
    "blastn",
    "-db", "nt",
    "-query", top_10_file_path,
    "-word_size", "11",
    "-gapopen", "6",
    "-gapextend", "4",
    "-max_target_seqs", "5",
    "-outfmt", "10 qaccver pident length evalue bitscore qcovs slen saccver stitle",
    "-remote"
]

# Output file path
blast_output_path = os.path.join(script_directory, "blast_results.csv")

# Execute the BLAST command
print(f"Starting BLAST search with {len(top_10_sequences)}. This may take some time...")
with open(blast_output_path, "w") as blast_output:
    subprocess.run(blast_command, stdout=blast_output)

print(f"BLAST results written to {blast_output_path}")

clean_blast_results('blast_results.csv')

#==============================================================================
# FIND BEST BARCODES
#==============================================================================

# Assume that blast_results.csv has been cleaned, and we'll read it into a DataFrame
blast_results = pd.read_csv('blast_results.csv')

# Filter the results
filtered_barcodes = blast_results.groupby('Query Accession Version').filter(filter_func)

# Get the unique barcodes
unique_barcodes = filtered_barcodes['Query Accession Version'].unique()

# Initialize the final DataFrame
FINAL_barcodes = pd.DataFrame(columns=[
    "barcode_name", "barcode_sequence", "max_MFE", "average_MFE",
    "Query Accession Version", "Percent Identity", "Alignment Length", 
    "E-value", "Bit Score", "Query Coverage", "Subject Length", 
    "Subject Accession Version", "Subject Title", "Location", "Type"
])

# Add the required information to the final DataFrame
for barcode_name in unique_barcodes:
    barcode_info = filtered_barcodes[filtered_barcodes['Query Accession Version'] == barcode_name].iloc[0]
    sequence_info = df[df['sequence_name'] == barcode_name].iloc[0]

    new_row = {
        "barcode_name": barcode_name,
        "barcode_sequence": sequence_info['whole_sequence'],
        "max_MFE": sequence_info['max_MFE'],
        "average_MFE": sequence_info['avg_MFE'],
        "Query Accession Version": barcode_info["Query Accession Version"],
        "Percent Identity": barcode_info["Percent Identity"],
        "Alignment Length": barcode_info["Alignment Length"],
        "E-value": barcode_info["E-value"],
        "Bit Score": barcode_info["Bit Score"],
        "Query Coverage": barcode_info["Query Coverage"],
        "Subject Length": barcode_info["Subject Length"],
        "Subject Accession Version": barcode_info["Subject Accession Version"],
        "Subject Title": barcode_info["Subject Title"],
        "Location": barcode_info["Location"],
        "Type": barcode_info["Type"]
    }
    new_row_df = pd.DataFrame([new_row])
    FINAL_barcodes = pd.concat([FINAL_barcodes, new_row_df], ignore_index=True)

# Optionally, save the final DataFrame to a CSV file
final_file_path = os.path.join(script_directory, "FINAL_barcodes.csv")
FINAL_barcodes.to_csv(final_file_path, index=False)

print(f"Final barcodes written to {final_file_path}")

# Create all possible pairs of barcodes
barcode_pairs = list(combinations(FINAL_barcodes['barcode_name'], 2))

# Store the ranked pairs
ranked_pairs = []

for pair in barcode_pairs:
    barcode1, barcode2 = pair
    barcode1_data = FINAL_barcodes.loc[FINAL_barcodes['barcode_name'] == barcode1].iloc[0]
    barcode2_data = FINAL_barcodes.loc[FINAL_barcodes['barcode_name'] == barcode2].iloc[0]

    # Calculate Hamming distance
    distance = hamming_distance(barcode1_data['barcode_sequence'], barcode2_data['barcode_sequence'])

    # Define the ranking criteria
    rank_criteria = [
        ('max_MFE', 'asc'),
        ('average_MFE', 'asc'),
        ('Alignment Length', 'asc'),
        ('Percent Identity', 'asc'),
        ('Query Coverage', 'asc')
    ]

    # Calculate ranks for each criteria
    ranks = []
    for column, order in rank_criteria:
        sorted_df = FINAL_barcodes.sort_values(column, ascending=(order == 'asc'))
        rank1 = sorted_df.index.get_loc(sorted_df[sorted_df['barcode_name'] == barcode1].index[0])
        rank2 = sorted_df.index.get_loc(sorted_df[sorted_df['barcode_name'] == barcode2].index[0])
        ranks.append((rank1 + rank2) / 2)

    # Average the ranks and store the results
    avg_rank = sum(ranks) / len(ranks)
    ranked_pairs.append((pair, avg_rank, distance))

# Sort the pairs by Hamming distance, then by average rank
ranked_pairs.sort(key=lambda x: (-x[2], x[1]))

# Create a DataFrame for the results
ranked_pairs_df = pd.DataFrame(ranked_pairs, columns=['barcode_pair', 'average_rank', 'hamming_distance'])

# Save ranked_pairs_df as a csv file
ranked_pairs_df.to_csv('ranked_barcode_pairs.csv', index=False)







