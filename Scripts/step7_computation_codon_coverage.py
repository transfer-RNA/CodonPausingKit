#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# --------------------------------------------------------
# --- CodonPausingKit                                  ---
# --- Copyright (c) 2025-2026 Aude Trinquier           ---
# --- All rights reserved.                             ---
# --- PolyForm Noncommercial License 1.0.0.            ---
# --- step7_computation_codon_coverage.py              ---
# --------------------------------------------------------

import sys
import pyBigWig
import pandas as pd
import os
import glob
from collections import defaultdict


# %% CODONS COVERAGE COMPUTATION
# -------------------------------

# Utility function for 7.1
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


# Definition of substep 7.1
def get_coverage_per_codon_for_multiple_bigwigs_in_folder(input_folder, bed_file, output_dir):
    """
    Process all BigWig files in a folder, applying the same BED file, and output results per BigWig file.
    
    Parameters:
    - input_folder: Path to the folder containing BigWig files.
    - bed_file: Path to the BED file with codon positions.
    - output_dir: Directory to save output files, one per BigWig file.
    """
    print("\n Getting coverage per codon for all the bigwig files...")
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
        output_file = os.path.join(output_dir, f"{sample_name}_codon_coverage.txt")
        
        # Process the current BigWig file
        print(f"  Processing {bigwig_file}...")
        get_coverage_per_codon(bigwig_file, bed_file, output_file)
        print(f"  -> Finished processing {bigwig_file} (for coverage per codon), output saved to {output_file}")



# %% CODON COVERAGE NORMALIZATION WITH GENE ELONGATION COVERAGE
# --------------------------------------------------------------

# Definition of substep 7.2
def codon_coverage_normalization(input_folder):
    print("\n Performing codon coverage normalization...")
    # Find all input files
    input_files = glob.glob(os.path.join(input_folder, "*_codon_coverage.txt"))

    # Process each file
    for file_path in input_files:
        # Read the input file
        df = pd.read_csv(file_path, sep="\t")

        # Calculate mean and standard deviation of Average_Coverage per Gene_name
        gene_stats = df.groupby("Gene_name")["Average_Coverage"].agg(
            GeneElongationCoverage='mean',
            GEC_SD='std'
        ).reset_index()

        # Merge back to original data
        df = pd.merge(df, gene_stats, on="Gene_name", how="left")

        # Calculate Normalized_Codon_Cov safely
        def safe_division(row):
            if row["GeneElongationCoverage"] == 0 or pd.isna(row["GeneElongationCoverage"]):
                return 0
            else:
                return row["Average_Coverage"] / row["GeneElongationCoverage"]

        df["Normalized_Codon_Cov"] = df.apply(safe_division, axis=1)

        # Fill NaN GEC_SD with 0 (in case only one codon was present for a gene)
        df["GEC_SD"] = df["GEC_SD"].fillna(0)

        # Prepare output filename
        filename = os.path.basename(file_path)
        output_filename = filename.replace("_codon_coverage.txt", "_codon_coverage_with_normalization.txt")
        output_path = os.path.join(input_folder, output_filename)

        # Save the output
        df.to_csv(output_path, sep="\t", index=False)

        print(f"  Normalized {filename} -> {output_filename}")

    print("All files normalized.")



# %% ACROSS GENOME CODON COVERAGE COMPUTATION (normalized)
# --------------------------------------------------------

# Utility function for 7.3
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


# Definition of substep 7.3
def compute_average_coverage_for_multiple_files(input_folder, output_folder):
    print("\n Computing average coverage for multiple files...")
    # Check if output folder exists, if not, create it
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Iterate through all .txt files in the input folder
    for filename in os.listdir(input_folder):
        if filename.endswith('normalization.txt'):
            input_file = os.path.join(input_folder, filename)
            output_file = os.path.join(output_folder, f'average_{filename}')
            print(f'  Processing {filename} for average coverage')
            compute_average_coverage(input_file, output_file)
            print(f'  Output written to {output_file}')

    print(f'Processing complete. All output files are in {output_folder}')



# %% Added November 2024: across the board normalization
# -------------------------------------------------------

# Utility function for 7.4
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
    print(f"  Normalized (across the board) file saved as: {output_file}")


# Definition of substep 7.4
def normalize_codon_coverage_for_multiple_files(input_folder, output_folder):
    print("\n Across the board normalization of codon coverage...")
    # Check if output folder exists, if not, create it
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Iterate over all files ending with "_codon_coverage_with_normalization.txt" in the input folder
    for file in os.listdir(input_folder):
        if file.endswith("_codon_coverage_with_normalization.txt"):
            input_file = os.path.join(input_folder, file)
            print(f"  Normalizing codon coverage of: {input_file}")
            normalize_codon_coverage(input_file, output_folder)



# %% SANITY CHECK FILE LENGTH
# ---------------------------

# Definition of substep 7.5
def count_lines_in_txt_files(folder_path):
    print("\n Small sanity check: printing file length...")
    # List all files in the folder
    files = [f for f in os.listdir(folder_path) if f.endswith('.txt')]
    
    # Iterate over each .txt file and count lines
    for file in files:
        file_path = os.path.join(folder_path, file)
        with open(file_path, 'r') as f:
            line_count = sum(1 for _ in f)
        print(f"  {file}: {line_count} lines")


def main():
    if len(sys.argv) != 6:
        print("Usage: step7_computation_codon_coverage.py <bigwig_cleaned_from_step6> <input_bed_codons_to_analyze> <codon_computations_folder_name> <results_codon_coverage_folder_name> <normalized_codon_coverage_folder_name>")
        sys.exit(1)

    bigwig_cleaned_from_step6_folder_name = sys.argv[1]  # input needed for step 7.1
    input_bed_codons_to_analyze = sys.argv[2]            # input needed for step 7.1
    codon_computations_folder_name = sys.argv[3]         # (output) folder for step 7.1, which will be also needed for step 7.2 (as input)
    results_codon_coverage_folder_name = sys.argv[4]     # (output) folder for step 7.3, which will be also needed for step 7.4 (as input)
    normalized_codon_coverage_folder_name = sys.argv[5]  # (output) folder for step 7.4, which will be also needed for the sanity check of step 7.5 (as input)

    # Calling substep 7.1) Getting the coverage per codon for multiple bigwig files
    get_coverage_per_codon_for_multiple_bigwigs_in_folder(bigwig_cleaned_from_step6_folder_name, input_bed_codons_to_analyze, codon_computations_folder_name)

    # Calling substep 7.2) Normalizing the codon coverage
    codon_coverage_normalization(codon_computations_folder_name)

    # Calling substep 7.3) Computing average coverage
    compute_average_coverage_for_multiple_files(codon_computations_folder_name, results_codon_coverage_folder_name)

    # Calling sustep 7.4) Normalizing codon coverage "across the board"
    # I'm commented this step fow now since I think it's not used in the end
    normalize_codon_coverage_for_multiple_files(results_codon_coverage_folder_name, normalized_codon_coverage_folder_name)

    # Calling substep 7.5) Doing a little sanity check
    count_lines_in_txt_files(normalized_codon_coverage_folder_name) # Updated too to print the 62 lines on the last one, just like for RNA side


if __name__ == "__main__":
    main()
