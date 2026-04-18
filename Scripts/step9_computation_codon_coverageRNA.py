# --------------------------------------------------------
# --- CodonPausingKit                                  ---
# --- Copyright (c) 2025-2026 Aude Trinquier           ---
# --- All rights reserved.                             ---
# --- PolyForm Noncommercial License 1.0.0.            ---
# --- step9_computation_codon_coverageRNA.py           ---
# --------------------------------------------------------

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pyBigWig
import pandas as pd
import os
import glob
import re
from collections import defaultdict


# Note: this step 9 is the equivalent of step 7, but this time for RNA (instead of RIBO)


# %% CODONS COVERAGE COMPUTATION (new)
# ------------------------------------

# Utility function for 9.1
def get_coverage_per_codon(bigwig_file, bed_file, output_file):
    """
    Compute average coverage per codon based on BigWig data and codon positions in BED file,
    and count how many times a coverage occurrence is found for each codon.
    
    Parameters:
    - bigwig_file: Path to the input BigWig file with per-nt coverage data.
    - bed_file: Path to the BED file with codon positions (3-nt intervals).
    - output_file: Path to save the output with codon coverage results and the occurrence counter.
    """
    # Open the BigWig file
    bw = pyBigWig.open(bigwig_file)
    
    # Load the BED file as a DataFrame
    bed_data = pd.read_csv(bed_file, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand', 'Gene_name'])
    
    # Prepare the output file
    with open(output_file, 'w') as out_file:
        out_file.write('Chromosome\tCodon_Start\tCodon_End\tCodon_Seq\tAverage_Coverage\tOccurrences\tGene_name\n')
        
        # Loop through each codon in the BED file
        for index, row in bed_data.iterrows():
            chrom = row['Chromosome']
            codon_start = row['Start']
            codon_end = row['End']
            codon_seq = row['Name']
            Gene_name = row['Gene_name']
            
            # Query the coverage for the 3-nt codon range
            coverage_values = bw.values(chrom, codon_start, codon_end)
            
            # Initialize occurrence counter
            occurrences = 0
            
            # Calculate the average coverage for the codon
            valid_coverage = [val for val in coverage_values if val is not None and val > 0]
            if valid_coverage:
                avg_coverage = sum(valid_coverage) / len(valid_coverage)
                occurrences = len(valid_coverage)  # Increment occurrences based on valid coverage points
            else:
                avg_coverage = 0  # Handle case where coverage is missing
            
            # Write the result to the output file
            out_file.write(f'{chrom}\t{codon_start}\t{codon_end}\t{codon_seq}\t{avg_coverage:.3f}\t{occurrences}\t{Gene_name}\n')
    
    # Close the BigWig file
    bw.close()


# Definition of substep 9.1
def get_coverage_per_codon_for_multiple_bigwigs_in_folder(input_folder, bed_file, output_dir):
    """
    Process all BigWig files in a folder, applying the same BED file, and output results per BigWig file.
    
    Parameters:
    - input_folder: Path to the folder containing BigWig files.
    - bed_file: Path to the BED file with codon positions.
    - output_dir: Directory to save output files, one per BigWig file.
    """
    print("\n Getting coverage per codon for all the bigwig files (RNA)...")
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all BigWig files in the input folder
    bigwig_files = glob.glob(os.path.join(input_folder, "*.bw"))
    
    if not bigwig_files:
        print(f"  No BigWig files found in the folder: {input_folder}")
        return
    
    # Loop through each BigWig file and process it
    for bigwig_file in bigwig_files:
        # Extract the sample name from the BigWig file for output naming
        sample_name = os.path.basename(bigwig_file).replace('.bw', '')
        
        # Define the output file name
        output_file = os.path.join(output_dir, f"{sample_name}_codon_RNAcoverage.txt")
        
        # Process the current BigWig file
        print(f"  Processing {bigwig_file}...")
        get_coverage_per_codon(bigwig_file, bed_file, output_file)
        print(f"  Finished processing {bigwig_file}, output saved to {output_file}")



# %% CODON RNA COVERAGE NORMALIZATION WITH RNA COVERAGE IN ORF (new)
# ------------------------------------------------------------------

# Definition of substep 9.2
def codon_coverage_normalization(input_folder):
    print("\n Performing codon coverage normalization (RNA)...")
    # Find all input files
    input_files = glob.glob(os.path.join(input_folder, "*_codon_RNAcoverage.txt"))

    # Process each file
    for file_path in input_files:
        # Read the input file
        df = pd.read_csv(file_path, sep="\t")

        # Calculate mean and standard deviation of Average_Coverage per Gene_name
        gene_stats = df.groupby("Gene_name")["Average_Coverage"].agg(
            GeneRNACoverage='mean',
            GRC_SD='std'
        ).reset_index()

        # Merge back to original data
        df = pd.merge(df, gene_stats, on="Gene_name", how="left")

        # Calculate Normalized_Codon_Cov safely
        def safe_division(row):
            if row["GeneRNACoverage"] == 0 or pd.isna(row["GeneRNACoverage"]):
                return 0
            else:
                return row["Average_Coverage"] / row["GeneRNACoverage"]

        df["Normalized_Codon_Cov"] = df.apply(safe_division, axis=1)

        # Fill NaN GEC_SD with 0 (in case only one codon was present for a gene)
        df["GRC_SD"] = df["GRC_SD"].fillna(0)

        # Prepare output filename
        filename = os.path.basename(file_path)
        output_filename = filename.replace("_codon_RNAcoverage.txt", "_codon_RNAcoverage_with_normalization.txt")
        output_path = os.path.join(input_folder, output_filename)

        # Save the output
        df.to_csv(output_path, sep="\t", index=False)

        print(f"Processed {filename} -> {output_filename}")

    print("All files processed.")



# %% ACROSS GENOME CODON RNA COVERAGE COMPUTATION (normalized)
# ------------------------------------------------------------

# Utility function for 9.3
def compute_average_coverage(input_file, output_file):
    codon_coverage = defaultdict(lambda: [0, 0])  # Dict to store [total coverage, total occurrences]
    
    with open(input_file, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            columns = line.strip().split('\t')  # Split by tab
            codon_seq = columns[3]  # Codon sequence (4th column)
            average_Normcoverage = float(columns[9])  # Average norm coverage (10th column)
            occurrences = int(columns[5])  # Occurrences (6th column)
            
            # Update total coverage and occurrences for each codon
            codon_coverage[codon_seq][0] += average_Normcoverage * occurrences
            codon_coverage[codon_seq][1] += occurrences
    
    # Compute the average coverage per codon sequence
    codon_avg_coverage = {}
    for codon, (total_coverage, total_occurrences) in codon_coverage.items():
        if total_occurrences > 0:
            codon_avg_coverage[codon] = total_coverage / total_occurrences
        else:
            codon_avg_coverage[codon] = 0.0

    # Write the results to the output file
    with open(output_file, 'w') as output:
        output.write('Codon_Seq\tAverage_NormCoverage\n')
        for codon, avg_coverage in codon_avg_coverage.items():
            output.write(f'{codon}\t{avg_coverage:.3f}\n')


# Definition of substep 9.3
def compute_average_coverage_for_multiple_files(input_folder, output_folder):
    print("\n Computing average coverage for multiple files (RNA)...")
    # Check if output folder exists, if not, create it
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Iterate through all .txt files in the input folder
    for filename in os.listdir(input_folder):
        if filename.endswith('normalization.txt'):
            input_file = os.path.join(input_folder, filename)
            output_file = os.path.join(output_folder, f'average_{filename}')
            print(f'Processing {filename}...')
            compute_average_coverage(input_file, output_file)
            print(f'Output written to {output_file}')



# %% Added November 2024: accross the board normalization
# Note: Not sure if should be used for RNA, since I think it wasn't used for RIBO
# --> Update with Aude: it should have been used for BOTH!
# -------------------------------------------------------

# # Utility function for 9.4
def normalize_codon_coverage(input_file, output_folder):
    # Read the input file into a DataFrame
    df = pd.read_csv(input_file, sep='\t')

    # Calculate the total coverage
    total_coverage = df['Average_NormCoverage'].sum()

    # Create a new DataFrame with Codon_Seq and Coverage_Percentage
    df_normalized = pd.DataFrame({
        'Codon_Seq': df['Codon_Seq'],
        'Coverage_Percentage': (df['Average_NormCoverage'] / total_coverage) * 100
    })

    # Generate output filename
    base = os.path.basename(input_file)
    output_file = os.path.join(output_folder, f"{os.path.splitext(base)[0]}_normalized.txt")

    # Save the new DataFrame to a file
    df_normalized.to_csv(output_file, sep='\t', index=False)
    print(f"  Normalized file saved as: {output_file}")


# Definition of substep 9.4
def normalize_codon_coverage_for_multiple_files(input_folder, output_folder):
    print("\n Across the board normalization of codon coverage (RNA)...")
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Iterate over all files ending with "_codon_RNAcoverage_with_normalization.txt" in the input folder
    for file in os.listdir(input_folder):
        if file.endswith("_codon_RNAcoverage_with_normalization.txt"):
            input_file = os.path.join(input_folder, file)
            normalize_codon_coverage(input_file, output_folder)



# %% SANITY CHECK FILE LENGTH
# Applied to the last step
# ---------------------------

# Definition of substep 9.5
def count_lines_in_txt_files(folder_path):
    # List all files in the folder
    files = [f for f in os.listdir(folder_path) if f.endswith('.txt')]
    
    # Iterate over each .txt file and count lines
    for file in files:
        file_path = os.path.join(folder_path, file)
        with open(file_path, 'r') as f:
            line_count = sum(1 for _ in f)
        print(f"{file}: {line_count} lines")



def main():
    if len(sys.argv) != 6:
        print("Usage: step9_computation_codon_coverageRNA.py <RNAbigwig_cleaned_from_step6> <input_bed_codons_to_analyze> <RNAcodon_computations_folder_name> <RNAresults_codon_coverage_folder_name> <RNAnormalized_codon_coverage_folder_name>")
        sys.exit(1)

    RNAbigwig_cleaned_from_step6_folder_name = sys.argv[1]  # input needed for step 9.1
    input_bed_codons_to_analyze = sys.argv[2]               # input needed for step 9.1
    RNAcodon_computations_folder_name = sys.argv[3]         # (output) folder for step 9.1, which will be also needed for step 9.2 (as input)
    RNAresults_codon_coverage_folder_name = sys.argv[4]     # (output) folder for step 9.3, which will be also needed for step 9.4 (as input)
    RNAnormalized_codon_coverage_folder_name = sys.argv[5]  # (output) folder for step 9.4, which will be also needed for the sanity check of step 9.5 (as input)

    # Calling substep 9.1) Getting the coverage per codon for multiple bigwig files
    get_coverage_per_codon_for_multiple_bigwigs_in_folder(RNAbigwig_cleaned_from_step6_folder_name, input_bed_codons_to_analyze, RNAcodon_computations_folder_name)

    # Calling substep 9.2) Normalizing the codon coverage
    codon_coverage_normalization(RNAcodon_computations_folder_name)

    # Calling substep 9.3) Computing average coverage
    compute_average_coverage_for_multiple_files(RNAcodon_computations_folder_name, RNAresults_codon_coverage_folder_name)

    # Calling sustep 9.4) Normalizing codon coverage "across the board"
    # I'm not sure that I should keep this step, since I think it was skipped for RIBO
    normalize_codon_coverage_for_multiple_files(RNAresults_codon_coverage_folder_name, RNAnormalized_codon_coverage_folder_name)

    # Calling substep 9.5) Doing a little sanity check
    count_lines_in_txt_files(RNAnormalized_codon_coverage_folder_name) # Npw printing the 62 for the real end, for RNA too


if __name__ == "__main__":
    main()
