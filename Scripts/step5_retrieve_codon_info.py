#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# --------------------------------------------------------
# --- CodonPausingKit                                  ---
# --- Copyright (c) 2025-2026 Aude Trinquier           ---
# --- All rights reserved.                             ---
# --- PolyForm Noncommercial License 1.0.0.            ---
# --- step5_retrieve_codon_info.py                     ---
# --------------------------------------------------------

import sys
import pandas as pd


# %% RETRIEVE CODON LIST FOR GENES TO ANALYZE
# -------------------------------------------

# Function to extract codon information for each gene
# Utility function to get_all_codon_info() defined right after
def get_codon_info(gene_name, nucleotide_seq, start, stop, strand, genome):
    codon_info = []
    codon_length = 3  # A codon is always 3 nucleotides long
    num_codons = len(nucleotide_seq) // codon_length
    
    for i in range(num_codons):
        codon_seq = nucleotide_seq[i * codon_length:(i + 1) * codon_length]
        codon_rank = i + 1
        
        if strand == "+":
            codon_start = start + i * codon_length
            codon_end = codon_start + codon_length - 1
        elif strand == "-":
            codon_end = stop - i * codon_length
            codon_start = codon_end - codon_length + 1

        codon_info.append({
            "Gene_name": gene_name,
            "Codon_rank": codon_rank,
            "Codon_seq": codon_seq,
            "Codon_start": codon_start,
            "Codon_end": codon_end,
            "Genome": genome  # Include the genome/chromosome information
        })
        
    return codon_info


# Definition of substep 5.1: getting all the codon info from the genes to analyze
def get_all_codon_info(input_csv_genes_to_analyze, output_codon_info_of_genes_to_analyze):
    print("\n Getting all the codon info...")
    # Load the CSV file into a pandas DataFrame
    df = pd.read_csv(input_csv_genes_to_analyze)

    # Loop over each filtered gene and extract codon information
    all_codon_info = []
    for _, row in df.iterrows():
        gene_name = row['Gene_name']
        nucleotide_seq = row['Nucleotide_seq']
        # Not useful it seems
        # aminoacid_seq = row['Aminoacid_seq']
        start = row['Start']  # Assumed to be gene start position
        stop = row['Stop']    # Assumed to be gene stop position
        strand = row['Strand']  # "+" or "-"
        genome = row['Genome']  # The chromosome/genome information
        
        codon_info = get_codon_info(gene_name, nucleotide_seq, start, stop, strand, genome)
        all_codon_info.extend(codon_info)

    # Convert all the codon info to a pandas DataFrame
    codon_df = pd.DataFrame(all_codon_info)

    # Save the result to a new CSV file
    codon_df.to_csv(output_codon_info_of_genes_to_analyze, index=False)

    print(f"  Codon information saved to: {output_codon_info_of_genes_to_analyze} with genome information")


# %% REMOVE FIRST TWO AND LAST TWO CODONS FOR EACH GENE
def remove_first2_and_last2_codons_of_each_gene(input_codon_info_of_genes_to_analyze, output_codons_to_analyze):
    # Load the codon information file generated in the previous step
    codon_df = pd.read_csv(input_codon_info_of_genes_to_analyze)

    # Function to remove the first 10 and last 5 codons for each gene
    def remove_codons(df):
        codons_to_analyze = []
        grouped = df.groupby("Gene_name")
        
        for gene_name, group in grouped:
            total_codons = len(group)
            # Ensure there are enough codons to remove
            if total_codons > 4:
                # Remove first 2 and last 2 codons
                filtered_group = group.iloc[2:total_codons-2]
                codons_to_analyze.append(filtered_group)
        
        return pd.concat(codons_to_analyze, ignore_index=True)

    # Apply the filtering function to the codon DataFrame
    codons_to_analyze_df = remove_codons(codon_df)

    # Save the result to a new CSV file
    codons_to_analyze_df.to_csv(output_codons_to_analyze, index=False)

    print(f"  Codon analysis file created (without first 2 and last 2 codons of each gene): {output_codons_to_analyze}")


# %% CONVERT CODON LIST TO BED FILE
def convert_to_bed(input_csv_file, output_bed_file):
    """
    Convert a CSV file containing codon information to BED format,
    retaining Gene_name as an additional column.
    
    Parameters:
    - input_csv_file: Path to the input CSV file.
    - output_bed_file: Path to save the output BED file.
    """
    # Load the codon data from CSV
    codon_data = pd.read_csv(input_csv_file)
    
    # Create a new DataFrame for the BED format
    bed_data = pd.DataFrame()
    
    # Populate the BED fields
    bed_data['Chromosome'] = codon_data['Genome']
    bed_data['Start'] = codon_data['Codon_start'] - 1  # BED uses 0-based start
    bed_data['End'] = codon_data['Codon_end']          # BED uses 1-based end
    bed_data['Name'] = codon_data['Codon_seq']         # Codon sequence as BED name
    bed_data['Score'] = 0                              # Placeholder score
    bed_data['Strand'] = '+'                           # Placeholder strand
    bed_data['Gene_name'] = codon_data['Gene_name']    # Additional column for gene name
    
    # Save to a BED file with tab-delimited format, no header or index
    bed_data.to_csv(output_bed_file, sep='\t', header=False, index=False)
    print(f"  Converted to BED file: {output_bed_file}")


def main():
    if len(sys.argv) != 5:
        print("Usage: step5_retrieve_codon_info.py <input_csv_genes_to_analyze> <output_codon_info_of_genes_to_analyze> <output_codons_to_analyze> <output_bed_codons_to_analyze>")
        sys.exit(1)

    input_csv_genes_to_analyze = sys.argv[1]            # input needed for step 5.1
    output_codon_info_of_genes_to_analyze = sys.argv[2] # output of step 5.1
    output_codons_to_analyze = sys.argv[3]              # output of step 5.2
    output_bed_codons_to_analyze = sys.argv[4]          # output of step 5.3

    # Calling substep 5.1) Getting all the codon info from the genes to analyze
    get_all_codon_info(input_csv_genes_to_analyze, output_codon_info_of_genes_to_analyze)

    input_codon_info_of_genes_to_analyze = output_codon_info_of_genes_to_analyze
    # Calling substep 5.2) Removing first and last 2 codons of each gene
    remove_first2_and_last2_codons_of_each_gene(input_codon_info_of_genes_to_analyze, output_codons_to_analyze)

    input_codons_to_analyze = output_codons_to_analyze
    # Calling substep 5.3) Conversion to bed file
    convert_to_bed(input_codons_to_analyze, output_bed_codons_to_analyze)


if __name__ == "__main__":
    main()
