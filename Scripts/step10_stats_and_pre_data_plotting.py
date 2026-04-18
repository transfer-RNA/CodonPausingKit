#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# --------------------------------------------------------
# --- CodonPausingKit                                  ---
# --- Copyright (c) 2025-2026 Aude Trinquier           ---
# --- All rights reserved.                             ---
# --- PolyForm Noncommercial License 1.0.0.            ---
# --- step10_stats_and_pre_data_plotting.py            ---
# --------------------------------------------------------

import sys
import os
import pandas as pd
import glob
import re
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt


# %% RENAME NORMALIZED FILES TO SIMPLIFIED SAMPLE NAMES
# -----------------------------------------------------

# Definition of substep 10.1
def rename_normalized_filed(input_output_folder):
    print("\n Renaming the (copied) normalized files for the stats and pre data plotting...")
    # Loop through all .txt files in the folder
    for filename in os.listdir(input_output_folder):
        if filename.endswith(".txt"):
            # Full path to the file
            full_path = os.path.join(input_output_folder, filename)

            # Check for average_Asite_RIBO pattern
            match_asite = re.match(r"average_(Asite_RIBO-[^\.]+)\.min.*\.txt", filename)
            if match_asite:
                new_name = match_asite.group(1) + ".txt"
                os.rename(full_path, os.path.join(input_output_folder, new_name))
                print(f"  Renamed to: {new_name}")
                continue  # Skip to next file

            # Check for average_RNA pattern
            match_rna = re.match(r"average_(RNA-[^\.]+)\.min.*\.txt", filename)
            if match_rna:
                new_name = match_rna.group(1) + ".txt"
                os.rename(full_path, os.path.join(input_output_folder, new_name))
                print(f"  Renamed to: {new_name}")



# %% AVERAGE PAIRED COVERAGE FILES (new for RNA and RIBO) (updated)
# -----------------------------------------------------------------

# Definition of substep 10.2
def average_paired_files_and_calculate_standard_deviation(input_output_folder):
    print("\n Computing the average of paired files and calculating the standard deviation...")
    # Create an output folder for averaged files
    output_folder = os.path.join(input_output_folder, "Averaged")
    os.makedirs(output_folder, exist_ok=True)

    # Get all .txt files
    file_list = glob.glob(os.path.join(input_output_folder, "*.txt"))

    # Group files by sample name (remove final -1 or -2 only)
    file_dict = {}
    for file_path in file_list:
        filename = os.path.basename(file_path)
        match = re.match(r"(.*)-([12])\.txt$", filename)
        if match:
            sample_prefix = match.group(1)  # Example: Asite_RIBO-Dt1 or RNA-Dt1
            replicate = match.group(2)      # "1" or "2"
            file_dict.setdefault(sample_prefix, {})[replicate] = file_path

    # Process each pair
    for sample_prefix, files in file_dict.items():
        if "1" in files and "2" in files:
            print(f"\n  Processing sample: {sample_prefix}")

            # Load both replicate files
            df1 = pd.read_csv(files["1"], sep="\t")
            df2 = pd.read_csv(files["2"], sep="\t")

            # Merge on Codon_Seq
            merged = pd.merge(df1, df2, on="Codon_Seq", suffixes=("_1", "_2"))

            # Check if Coverage_Percentage columns exist
            if "Coverage_Percentage_1" not in merged.columns or "Coverage_Percentage_2" not in merged.columns:
                print(f"  Skipping {sample_prefix}: Missing Coverage_Percentage columns.")
                continue

            # Compute mean and standard deviation of Coverage_Percentage
            merged["Coverage_Percentage_Mean"] = merged[["Coverage_Percentage_1", "Coverage_Percentage_2"]].mean(axis=1)
            merged["Coverage_Percentage_SD"] = merged[["Coverage_Percentage_1", "Coverage_Percentage_2"]].std(axis=1)

            # Keep Codon_Seq, Coverage_Percentage_Mean, and Coverage_Percentage_SD
            output_df = merged[["Codon_Seq", "Coverage_Percentage_Mean", "Coverage_Percentage_SD"]]

            # Save to new file
            output_filename = os.path.join(output_folder, f"{sample_prefix}.txt")
            output_df.to_csv(output_filename, sep="\t", index=False)
            print(f"  Averaged file written: {output_filename}")
        else:
            print(f"  Skipping incomplete pair for: {sample_prefix}")



# %% AA sorted files generation
# -----------------------------

# Definition of substep 10.3
def sort_in_amino_acid_order(input_output_folder):
    print("\n Sorting in amino acid order (and replacing Ts by Us)...")
    # Input/Output folder containing the averaged .txt files
    input_output_folder = f"{input_output_folder}/Averaged"
    # Desired codon order
    desired_order = [
        "GCA","GCC","GCG","GCU","AGA","AGG","CGA","CGC","CGG","CGU",
        "AAC","AAU","GAC","GAU","UGC","UGU","CAA","CAG","GAA","GAG",
        "GGA","GGC","GGG","GGU","CAC","CAU","AUA","AUC","AUU","CUA",
        "CUC","CUG","CUU","UUA","UUG","AAA","AAG","AUG","UUC","UUU",
        "CCA","CCC","CCG","CCU","AGC","AGU","UCA","UCC","UCG","UCU",
        "ACA","ACC","ACG","ACU","UGG","UAC","UAU","GUA","GUC","GUG",
        "GUU"
    ]

    # Create ranking dictionary for sorting
    order_dict = {codon: i for i, codon in enumerate(desired_order)}

    # Process files
    for file in os.listdir(input_output_folder):
        if file.endswith(".txt"):
            file_path = os.path.join(input_output_folder, file)

            # Read file (tab-delimited)
            df = pd.read_csv(file_path, sep="\t")

            # Replace T with U in Codon_Seq
            df["Codon_Seq"] = df["Codon_Seq"].str.replace("T", "U")

            # Add sorting key
            df["sort_key"] = df["Codon_Seq"].map(order_dict)

            # Check for unexpected codons
            if df["sort_key"].isnull().any():
                missing = df[df["sort_key"].isnull()]["Codon_Seq"].unique()
                print(f"Warning: {file} contains unexpected codons: {missing}")

            # Sort according to custom order
            df = df.sort_values("sort_key")

            # Drop helper column
            df = df.drop(columns=["sort_key"])

            # Output filename
            output_file = file.replace(".txt", "_AAsorted.txt")
            output_path = os.path.join(input_output_folder, output_file)

            # Save file
            df.to_csv(output_path, sep="\t", index=False)

            print(f"  Processed: {file} → {output_file}")

    print("All files processed.")



# %% CALCULATE RIBO/RNA COVERAGE RATIOS AND PROPAGATED STANDARD DEVIATION
# -----------------------------------------------------------------------

# Definition of substep 10.4
def calculate_ribo_rna_ratios_and_propagate_standard_deviation(input_output_folder):
    print("\n Calculating RIBO/RNA ratio and propagating the standard deviation of a quotient...")
    # Input/Output folder containing the averaged .txt files
    input_output_folder = f"{input_output_folder}/Averaged"

    # Output folder for the ratio files
    output_folder = os.path.join(input_output_folder, "Ratios_with_SD")
    os.makedirs(output_folder, exist_ok=True)

    # Get all _AAsorted.txt files
    file_list = glob.glob(os.path.join(input_output_folder, "*_AAsorted.txt"))

    # Organize files by sample (Dt1, Dt2, Wt1, Wt2, etc.)
    samples = {}

    for file_path in file_list:
        filename = os.path.basename(file_path)
        match = re.match(r"(Asite_RIBO|RNA)-(.*)\_AAsorted.txt$", filename)
        if match:
            file_type = match.group(1)  # Asite_RIBO or RNA
            sample_name = match.group(2)  # Dt1, Dt2, Wt1, Wt2
            samples.setdefault(sample_name, {})[file_type] = file_path

    # Compute ratios and propagated SD
    for sample_name, file_paths in samples.items():
        if "Asite_RIBO" in file_paths and "RNA" in file_paths:
            print(f"  Processing sample: {sample_name}")

            # Load the RIBO and RNA data
            ribo_df = pd.read_csv(file_paths["Asite_RIBO"], sep="\t")
            rna_df = pd.read_csv(file_paths["RNA"], sep="\t")

            # Merge based on Codon_Seq
            merged = pd.merge(ribo_df, rna_df, on="Codon_Seq", suffixes=("_RIBO", "_RNA"))

            # Calculate the ratio
            merged["RIBO_to_RNA_Coverage_Ratio"] = merged["Coverage_Percentage_Mean_RIBO"] / merged["Coverage_Percentage_Mean_RNA"]

            # Calculate the propagated standard deviation
            merged["RIBO_to_RNA_Coverage_Ratio_SD"] = merged["RIBO_to_RNA_Coverage_Ratio"] * (
                ( (merged["Coverage_Percentage_SD_RIBO"] / merged["Coverage_Percentage_Mean_RIBO"])**2 +
                (merged["Coverage_Percentage_SD_RNA"] / merged["Coverage_Percentage_Mean_RNA"])**2
                ) ** 0.5
            )

            # (Optional) Mask ratios and SDs where RNA coverage is very low (e.g., < 1%) to avoid noisy divisions
            threshold = 1
            mask = merged["Coverage_Percentage_Mean_RNA"] < threshold
            merged.loc[mask, ["RIBO_to_RNA_Coverage_Ratio", "RIBO_to_RNA_Coverage_Ratio_SD"]] = float('nan')

            # Create output file
            output_df = merged[["Codon_Seq", "RIBO_to_RNA_Coverage_Ratio", "RIBO_to_RNA_Coverage_Ratio_SD"]]
            output_filename = os.path.join(output_folder, f"Ratio_{sample_name}.txt")
            output_df.to_csv(output_filename, sep="\t", index=False)
            print(f"  Ratio + SD file written: {output_filename}")

        else:
            print(f"  Skipping sample {sample_name}: missing RIBO or RNA file")



# %% Z-score with quick plot
# --------------------------

# Definition of substep 10.5
def calculate_z_score_and_generate_quick_plot(input_output_folder):
    print("\n Calculating z-score and generating a quick plot...")
    # Input/Output folder containing the averaged .txt files
    input_output_folder = f"{input_output_folder}/Averaged/Ratios_with_SD"

    # Define input file names
    files = {
        "t1": ("Ratio_Dt1.txt", "Ratio_Wt1.txt"),
        "t2": ("Ratio_Dt2.txt", "Ratio_Wt2.txt")
    }

    # Function to plot Z-scores
    def plot_z_scores(df, timepoint, output_folder):
        plt.figure(figsize=(12, 6))
        df_sorted = df.sort_values("Z", ascending=False)
        plt.bar(df_sorted['Codon_Seq'], df_sorted['Z'])
        plt.axhline(0, color='gray', linestyle='--')
        plt.xticks(rotation=90)
        plt.ylabel("Z-score")
        plt.title(f"Codon Z-scores - {timepoint.upper()}")
        plt.tight_layout()

        # Save plot
        plot_file = os.path.join(output_folder, f"Zscore_plot_{timepoint}.png")
        plt.savefig(plot_file, dpi=300)
        plt.close()
        print(f"  Plot saved to: {plot_file}")

    # Process each timepoint
    for timepoint, (d_file, wt_file) in files.items():
        # Load data
        d_path = os.path.join(input_output_folder, d_file)
        wt_path = os.path.join(input_output_folder, wt_file)

        df_d = pd.read_csv(d_path, sep="\t")
        df_wt = pd.read_csv(wt_path, sep="\t")

        # Rename columns
        df_d.columns = ['Codon_Seq', 'Ratio_D', 'SD_D']
        df_wt.columns = ['Codon_Seq', 'Ratio_WT', 'SD_WT']

        # Merge
        merged = pd.merge(df_d, df_wt, on='Codon_Seq')

        # Z-score calculation
        merged['Z'] = (merged['Ratio_D'] - merged['Ratio_WT']) / np.sqrt(merged['SD_D']**2 + merged['SD_WT']**2)
        merged['p_value'] = 2 * (1 - norm.cdf(np.abs(merged['Z'])))

        # BH correction
        merged = merged.sort_values('p_value').reset_index(drop=True)
        m = len(merged)
        merged['rank'] = np.arange(1, m + 1)
        merged['adj_p'] = merged['p_value'] * m / merged['rank']
        merged['adj_p'] = np.minimum(merged['adj_p'], 1.0)

        # Save data output
        output_file = os.path.join(input_output_folder, f"Codon_Zscore_Results_{timepoint}.csv")
        merged.to_csv(output_file, sep="\t", index=False)
        print(f"  Results saved to: {output_file}")

        # Plot Z-scores
        plot_z_scores(merged, timepoint, input_output_folder)


def main():
    if len(sys.argv) != 2:
        print("Usage: step10_stats_and_pre_data_plotting.py <Asite_pausing_calculations_folder_name>")
        sys.exit(1)

    Asite_pausing_calculations_folder_name = sys.argv[1]

    # Calling substep 10.1 that renames the files of interest for easier later manipulations
    rename_normalized_filed(Asite_pausing_calculations_folder_name)

    # Calling substep 10.2 that averages paired files and calculates standard deviation
    average_paired_files_and_calculate_standard_deviation(Asite_pausing_calculations_folder_name)

    # Calling substep 10.3 to order each file, in amino acid order (with changes T -> U to reflect mRNA code)
    sort_in_amino_acid_order(Asite_pausing_calculations_folder_name)

    # Calling substep 10.4 to calculate RIBO/RNA ratios and to propagate standard deviation
    calculate_ribo_rna_ratios_and_propagate_standard_deviation(Asite_pausing_calculations_folder_name)

    # Calling substep 10.5 to calculate z-score and to generate a quick plot
    calculate_z_score_and_generate_quick_plot(Asite_pausing_calculations_folder_name)


if __name__ == "__main__":
    main()
