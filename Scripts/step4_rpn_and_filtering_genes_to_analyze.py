# --------------------------------------------------------
# --- CodonPausingKit                                  ---
# --- Copyright (c) 2025-2026 Aude Trinquier           ---
# --- All rights reserved.                             ---
# --- PolyForm Noncommercial License 1.0.0.            ---
# --- step4_rpn_and_filtering_genes_to_analyze.py      ---
# --------------------------------------------------------

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import pandas as pd
from collections import defaultdict
import csv
import gffutils


# Definition of substep 4.1)
def make_identifiers_unique(input_csv, output_csv, report_file):
    """
    Making Gene_name identifier a unique identifier
    """
    print("\n Making unique identifiers...")
    # Load CSV
    df = pd.read_csv(input_csv)
    if 'Gene_name' not in df.columns:
        raise ValueError("  Error: 'Gene_name' column not found.")

    # Track gene name occurrences
    name_count = defaultdict(int)
    modifications = []

    unique_names = []
    for name in df['Gene_name']:
        name_count[name] += 1
        if name_count[name] == 1:
            unique_names.append(name)
        else:
            new_name = f"{name}_{name_count[name]-1}"
            unique_names.append(new_name)
            modifications.append(f"{name} -> {new_name}")

    # Update the DataFrame
    df['Gene_name'] = unique_names

    # Save new CSV
    df.to_csv(output_csv, index=False)

    # Save modifications report
    with open(report_file, 'w') as f:
        f.write("Modified Gene_name entries:\n")
        for line in modifications:
            f.write(line + "\n")

    print(f"  Output CSV written to: {output_csv}")
    print(f"  Modifications saved to: {report_file}")


# Definition of substep 4.2: Reads Per Nucleotide
def calculate_rpn(input_csv_with_unique_ids, input_csv_mapped_reads, output_rpn_file):
    """
    Making RPN calculations
    """
    print("\n Making RPN calculations...")
    # Read the data from the CSV files into DataFrames
    df = pd.read_csv(input_csv_with_unique_ids)
    mapped_reads_df = pd.read_csv(input_csv_mapped_reads)

    # Get the first 17 columns
    first_17_columns = df.columns[:17]

    # Filter columns that match the pattern "RIBO-xxx-x_rpkm"
    ribo_columns = [col for col in df.columns if col.startswith("RIBO-") and col.endswith("_rpkm")]

    # Keep only the selected columns
    selected_columns = list(first_17_columns) + ribo_columns
    df_short = df[selected_columns].copy()  # Explicitly make a copy to avoid SettingWithCopyWarning

    # Read the total mapped reads for each sample
    total_mapped_reads = mapped_reads_df.set_index('Sample')['Total_Mapped_Reads'].to_dict()

    # Perform the reads per nucleotide (RPN) calculation and append new columns
    for ribo_column in ribo_columns:
        sample = ribo_column.split('_')[0]  # Extract the sample name
        if sample in total_mapped_reads:
            df_short.loc[:, f'{sample}_rpn'] = (df_short[ribo_column] * total_mapped_reads[sample]) / 1e9

    # Save the modified DataFrame to a new CSV file
    df_short.to_csv(output_rpn_file, index=False)

    print(f"  File '{output_rpn_file}' created successfully with RPN calculations.")


# Definition of substep 4.3 : exporting the overlapping cds
def export_overlapping_cds(gff_input, output_csv):
    """
    Identify CDS features that overlap in a GFF file and export them to CSV.

    Parameters
    ----------
    gff_input : str
        Path to the input GFF file.
    output_csv : str
        Path to the output CSV file.

    Notes
    -----
    The output file contains two columns:
    - CDS: the name of the CDS feature
    - Overlapping CDS: comma-separated names of overlapping CDS features

    The function uses the 'gene' attribute when available, otherwise falls
    back to 'ID', then to the internal feature id.
    """
    print("\n Computing the overlapping cds...")

    # Helper function to get CDS gene attribute or fallback to ID if gene is not available
    def get_feature_name(feature):
        if "gene" in feature.attributes:
            return feature.attributes["gene"][0] # Use gene if present
        elif "ID" in feature.attributes:
            return feature.attributes["ID"][0] # Fallback to ID
        else:
            return feature.id # Fallback to feature id

    print(f"  Loading GFF file: {gff_input}")

    db = gffutils.create_db(
        gff_input,
        dbfn=":memory:",
        force=True,
        keep_order=True,
        merge_strategy="merge",
        sort_attribute_values=True
    )

    overlapping_cds = []

    for feature in db.features_of_type("CDS"):
        overlaps = db.region(
            region=(feature.seqid, feature.start, feature.end),
            completely_within=False,
            featuretype="CDS"
        )

        overlap_names = [
            get_feature_name(other_feature)
            for other_feature in overlaps
            if other_feature.id != feature.id
        ]

        if overlap_names:
            overlapping_cds.append((get_feature_name(feature), overlap_names))

    with open(output_csv, mode="w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["CDS", "Overlapping CDS"])
        for cds_name, overlap_names in overlapping_cds:
            writer.writerow([cds_name, ", ".join(overlap_names)])

    print(f"  Overlapping CDS written to: {output_csv}")


# Definition of substep 4.4: Filtering out CDS that overlap
def filter_overlapping_cds(input_rpn_csv, input_overlapping_cds_csv, output_filtered_csv):
    """
    Remove genes from the RPN table whose Gene_name appears in the list of
    overlapping CDS.

    Parameters
    ----------
    input_rpn_csv : str
        Path to the input CSV file containing RPN calculations.
    input_overlapping_cds_csv : str
        Path to the CSV file listing overlapping CDS.
    output_filtered_csv : str
        Path to the output CSV file after filtering.
    """
    print("\n Filtering overlapping cds...")

    overview_df = pd.read_csv(input_rpn_csv)
    overlapping_df = pd.read_csv(input_overlapping_cds_csv)

    if "Gene_name" not in overview_df.columns:
        raise ValueError("  Error: 'Gene_name' column not found in RPN file.")

    if "CDS" not in overlapping_df.columns:
        raise ValueError("  Error: 'CDS' column not found in overlapping CDS file.")

    # Extract the list of CDS names from the 'overlapping_cds.csv' file
    cds_to_filter = overlapping_df['CDS'].tolist()

    # Remove entries that are present in the overlapping CDS
    filtered_df = overview_df[~overview_df['Gene_name'].isin(cds_to_filter)].copy() # .copy() is new for calculating the difference

    n_removed = len(overview_df) - len(filtered_df)
    print(f"  Removed {n_removed} rows matching overlapping CDS.")
    
    # Save the remaining (non-overlapping) entries to a new CSV file
    filtered_df.to_csv(output_filtered_csv, index=False)
    print(f"  Filtered overview written to: {output_filtered_csv}")


# Definition of substep 4.5
# %% FILTERING OUT LOW RPN GENES
def filtering_low_rpn_genes(input_filtered_csv, output_filtered_low_rpn, output_filtered_out_low_rpn):
    print("\n Filtering low rpn genes...")

    # Read the previous output file into a DataFrame
    df_rpn = pd.read_csv(input_filtered_csv)

    # Get all RPN columns (those ending with '_rpn')
    rpn_columns = [col for col in df_rpn.columns if col.endswith('_rpn')]

    # New, to be more consistent in the future
    rpn_threshold = 0.33

    # Filter out rows where any RPN value is under 0.33
    filtered_out = df_rpn[(df_rpn[rpn_columns] < rpn_threshold).any(axis=1)].copy()
    df_filtered = df_rpn[(df_rpn[rpn_columns] >= rpn_threshold).all(axis=1)].copy()

    # For the filtered-out genes, find the samples that didn't meet the criteria
    # NOTE: Now using "< rpn_threshold" instead of the original leftover "<= 0.5", to be consistent with how we've filtered just above
    filtered_out.loc[:, 'Failed_Samples'] = filtered_out[rpn_columns].apply(lambda row: ', '.join([col for col in rpn_columns if row[col] < rpn_threshold]), axis=1)

    # Save the filtered-out genes to a file
    filtered_out.to_csv(output_filtered_out_low_rpn, columns=['Gene_name', 'Failed_Samples'], index=False)

    # Save the filtered DataFrame to a new file
    df_filtered.to_csv(output_filtered_low_rpn, index=False)

    print(f"  File '{output_filtered_low_rpn}' created successfully with filtered genes (no low rpn).")
    print(f"  File '{output_filtered_out_low_rpn}' created with the samples that didn't meet the rpn threshold criteria.")


# Definition of substep 4.6: more filtering
# %% FILTERING OUT GENES WITH an empty Aminoacid_seq OR where Codon_count is not an integer
def filtering_empty_aminoacid_or_codon_count_thats_not_integer(input_filtered_low_rpn, output_fully_filtered, output_empty_aminoacid_seq, output_codon_count_not_integer):
    print("\n Filtering empty aminoacid or codon count that's not integer...")

    # Load the CSV file from previous step into a pandas DataFrame
    df = pd.read_csv(input_filtered_low_rpn)

    # Selecting rows with an empty Aminoacid_seq and save them
    # (it uses the operator | on a panda boolean series to simulate a logical 'or')
    empty_aminoacid_seq = df[df['Aminoacid_seq'].isnull() | df['Aminoacid_seq'].str.strip().eq('')]
    empty_aminoacid_seq.to_csv(output_empty_aminoacid_seq, index=False)
    # Only keeps rows with non-empty Aminoacid_seq
    # (it uses the operator & on a panda boolean series to simulate a logical 'and')
    df_filtered_1 = df[df['Aminoacid_seq'].notnull() & df['Aminoacid_seq'].str.strip().ne('')]

    # Selecting rows where Codon_count is not an integer and save them
    non_integer_codon_count = df_filtered_1[df_filtered_1['Codon_count'] != df_filtered_1['Codon_count'].astype(int)]
    non_integer_codon_count.to_csv(output_codon_count_not_integer, index=False)
    # Only keeps rows where Codon_count is an integer
    df_filtered_2 = df_filtered_1[df_filtered_1['Codon_count'] == df_filtered_1['Codon_count'].astype(int)]

    # Save the final retained rows into a new CSV file
    df_filtered_2.to_csv(output_fully_filtered, index=False)

    print("  Filtered out rows saved to", output_empty_aminoacid_seq, "and to", output_codon_count_not_integer)
    print("  Retained rows after filtering saved to", output_fully_filtered)


# %% SANITY CHECK ON Gene_name DUPLICATES
def sanity_check_no_gene_name_duplicated(input_filtered_csv_to_check):
    print("\n Sanity check 1: verifying that there is now no Gene_name entry duplicated...")

    # Load the file
    df = pd.read_csv(input_filtered_csv_to_check)

    # Check if 'Gene_name' column exists
    if 'Gene_name' not in df.columns:
        raise ValueError("  Error: 'Gene_name' column not found in the file.")

    # Find duplicates
    duplicates = df[df['Gene_name'].duplicated(keep=False)]

    if duplicates.empty:
        print("  -> No duplicate Gene_name entries found.")
    else:
        print(f"  -> Found {duplicates['Gene_name'].nunique()} duplicate Gene_name entries.\n")
        print("  Here are the duplicates:")
        print(duplicates.sort_values(by='Gene_name'))


def main():
    if len(sys.argv) != 14:
        print("Usage: step4_rpn_and_filtering_genes_to_analyze.py <input_csv_intermediate> <input_mapped_reads> <input_gff> <output_csv_with_unique_ids> <output_report> <output_rpn> <output_overlapping_cds> <output_filtered_overlap> <output_filtered_low_rpn> <output_filtered_out_low_rpn> <output_fully_filtered> <output_empty_aminoacid_seq> <output_codon_count_not_integer>")
        sys.exit(1)

    csv_intermediate = sys.argv[1]           # input needed for step 4.1
    mapped_reads_input = sys.argv[2]         # input needed for step 4.2
    gff_input = sys.argv[3]                  # input needed for step 4.3
    output_csv_with_unique_ids = sys.argv[4] # output of step 4.1
    output_report_file = sys.argv[5]         # output of step 4.1
    output_rpn = sys.argv[6]                 # output of step 4.2
    output_overlapping_cds = sys.argv[7]     # output of step 4.3
    output_filtered_overlap = sys.argv[8]    # output of step 4.4
    output_filtered_low_rpn = sys.argv[9]      # output (main) of step 4.5
    output_filtered_out_low_rpn = sys.argv[10] # output (secondary) of step 4.5
    output_fully_filtered = sys.argv[11]          # output (main) of step 4.6
    output_empty_aminoacid_seq = sys.argv[12]     # output (secondary) of step 4.6
    output_codon_count_not_integer = sys.argv[13] # output (secondary) of step 4.6

    # Calling substep 4.1) Making identifiers unique
    # This intermediate steps consumes the csv_intermediate produced in mini step 3, and produces output_csv_with_unique_ids and output_report_file
    make_identifiers_unique(csv_intermediate, output_csv_with_unique_ids, output_report_file)

    # Calling subset 4.2) to compute Reads Per Nucleotide
    input_csv_with_unique_ids = output_csv_with_unique_ids
    calculate_rpn(input_csv_with_unique_ids, mapped_reads_input, output_rpn) # Note (14 Mar 2026): I changed the first argument from csv_intermediate since we need to operate on the one with unique ids I think

    # Calling subset 4.3) to identify CDS features that overlap in a GFF file and export them to CSV
    export_overlapping_cds(gff_input, output_overlapping_cds)

    # Calling substep 4.4) to remove overlapping CDS from the RPN table
    input_rpn = output_rpn
    filter_overlapping_cds(input_rpn, output_overlapping_cds, output_filtered_overlap)

    # Verify that the output of step 4.4 has no Gene_name duplication
    sanity_check_no_gene_name_duplicated(output_filtered_overlap)

    input_filtered_overlap = output_filtered_overlap
    # Calling substep 4.5) to filter the low rpn
    filtering_low_rpn_genes(input_filtered_overlap, output_filtered_low_rpn, output_filtered_out_low_rpn)
    
    input_filtered_low_rpn = output_filtered_low_rpn
    # Calling substep 4.6) to filter out genes with an empty aminoacid sequence or with a codon count that's not integer
    filtering_empty_aminoacid_or_codon_count_thats_not_integer(input_filtered_low_rpn, output_fully_filtered, output_empty_aminoacid_seq, output_codon_count_not_integer)

    # Verify the output of step 4.6 has no Gene_name duplication
    sanity_check_no_gene_name_duplicated(output_fully_filtered)


if __name__ == "__main__":
    main()
