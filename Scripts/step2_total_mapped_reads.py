#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# --------------------------------------------------------
# --- CodonPausingKit                                  ---
# --- Copyright (c) 2025-2026 Aude Trinquier           ---
# --- All rights reserved.                             ---
# --- PolyForm Noncommercial License 1.0.0.            ---
# --- step2_total_mapped_reads.py                      ---
# --------------------------------------------------------

import sys
import os
import pandas as pd


def calculate_total_mapped_reads(excel_input, csv_output):
    # Load the Excel file
    df = pd.read_excel(excel_input)

    # Get the path of the directory for the Intermediate files from the csv_output filename
    intermediate_dir = os.path.dirname(csv_output)
    # Construct the intermediate CSV filename by using the Excel filename, but placed in Intermediates
    excel_base = os.path.splitext(os.path.basename(excel_input))[0]
    csv_intermediate = os.path.join(intermediate_dir, excel_base + ".csv")

    # Save the DataFrame as a CSV file
    df.to_csv(csv_intermediate, index=False)
    print(f"Intermediate CSV has been saved to {csv_intermediate}")

    # Filter for CDS and pseudogene features (assuming these features are indicated in a column named 'Feature')
    cds_pseudogene_df = df[df['Class'].isin(['CDS', 'pseudogene'])]

    # Extract RIBO samples 
    ribo_samples = [col for col in cds_pseudogene_df.columns if 'RIBO' in col]

    # Sum the counts for each RIBO sample
    total_mapped_reads = cds_pseudogene_df[ribo_samples].sum()

    # Create a simple output table
    output_df = pd.DataFrame({
        'Sample': total_mapped_reads.index,
        'Total_Mapped_Reads': total_mapped_reads.values
    })

    # Save the output table to a new CSV file
    output_df.to_csv(csv_output, index=False)

    print(f"Total mapped reads have been calculated and saved to {csv_output}")


def main():
    if len(sys.argv) != 3:
        print("Usage: step2_total_mapped_reads.py <excel_input> <csv_output>")
        sys.exit(1)

    excel_input, csv_output = sys.argv[1], sys.argv[2]
    calculate_total_mapped_reads(excel_input, csv_output)


if __name__ == "__main__":
    main()
