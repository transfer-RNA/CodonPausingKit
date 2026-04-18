#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# --------------------------------------------------------
# --- CodonPausingKit                                  ---
# --- Copyright (c) 2025-2026 Aude Trinquier           ---
# --- All rights reserved.                             ---
# --- PolyForm Noncommercial License 1.0.0.            ---
# --- tools.py                                         ---
# --------------------------------------------------------

# This file implements some common functionality between step 6 (on the RIBO side)
# and step 8 (on the RNA side), to avoid duplication of code.

import os
import subprocess
import csv
import glob


# %% CONVERTING BIGWIG INTO WIGGLE
# Used by substeps 6.2 and 8.1
# --------------------------------

def convert_bigwig_to_wig(bigwig_file, output_directory):
    # Define the output file path
    output_wig = os.path.join(output_directory, os.path.basename(bigwig_file).replace(".bw", ".wig"))

    # Construct the command to convert BigWig to Wig
    command = ["bigWigToWig", bigwig_file, output_wig]

    try:
        # Run the conversion command
        subprocess.run(command, check=True)
        print(f"  Converted {bigwig_file} to {output_wig}")
    except subprocess.CalledProcessError as e:
        print(f"  Error converting {bigwig_file}: {e}")


def batch_convert_bigwig_to_wig(bigwig_files, output_directory):
    for bigwig_file in bigwig_files:
        convert_bigwig_to_wig(bigwig_file, output_directory)


def discover_bigwig_files(bigwig_folder_name):
    """
    Discover all BigWig files in the given folder.
    """
    bigwig_files = sorted(glob.glob(os.path.join(bigwig_folder_name, "*.bw")))

    if not bigwig_files:
        raise ValueError(f"  Error: No .bw files found in {bigwig_folder_name}")

    return bigwig_files


def conversion_from_bigwig_to_wig(bigwig_folder_name, wig_folder_name):
    """
    Convert from bigwig to wiggle
    
    Parameters:
    - bigwig_folder_name: Path to the folder where the input files will be automatically dicovered
    - wig_folder_name: Path to the folder that will hold the produced Wig files.
    """
    print("\n Conversion from bigwig to wig...")
    # Ensure the output directory exists
    if not os.path.exists(wig_folder_name):
        os.makedirs(wig_folder_name)

    # Read the list of BigWig files from the CSV file
    bigwig_files = discover_bigwig_files(bigwig_folder_name)

    # Batch convert all BigWig files to Wig
    batch_convert_bigwig_to_wig(bigwig_files, wig_folder_name)



# %% GETTING ABSOLUTE COVERAGE IN THE REVERSE FILES
# Used by substeps 6.4 and 8.2
# -------------------------------------------------------------

# Function to adjust negative values in reverse wig files
def adjust_reverse_coverage(wig_file):
    """
    Convert negative coverage values to positive in one reverse wig file.
    The file is overwritten in place.
    """    
    with open(wig_file, 'r') as infile:
        lines = infile.readlines()
    
    adjusted_lines = []
    for line in lines:
        if line.startswith('track') or line.startswith('variableStep'):
            # Preserve the headers as they are
            adjusted_lines.append(line)
        else:
            # Process lines with position and coverage value (no end position)
            pos, value = line.strip().split()
            value = abs(float(value))  # Convert coverage value to positive
            adjusted_lines.append(f"{pos}\t{value}\n")
    
    # Overwrite the original file with adjusted values
    with open(wig_file, 'w') as outfile:
        outfile.writelines(adjusted_lines)


# Function to process all reverse wig files in a folder
def process_reverse_wig_files(shifted_wig_folder_name, pattern: str):
    """
    Process all reverse fiveprime wig files in a folder and convert their
    coverage values to absolute values.
    """
    print("\n Adjusting reverse wiggle coverage to absolute values...")

    processed_count = 0
    for filename in os.listdir(shifted_wig_folder_name):
        if filename.endswith(pattern):
            wig_file = os.path.join(shifted_wig_folder_name, filename)
            print(f"  Processing {wig_file}")
            adjust_reverse_coverage(wig_file)
            processed_count += 1

    print(f"  Processed {processed_count} reverse wig file(s).")


# %% TRANSFORMING WIGGLE FILES TO BEDGRAPH (going from 1-based to 0-based coordinates)
# Used by substeps 6.5 and 8.3
# ------------------------------------------------------------------------------------

def wiggle_to_bedgraph(wiggle_folder_name, bedgraph_folder_name):
    print("\n Wiggle to bedgraph...")
    if not os.path.exists(bedgraph_folder_name):
        os.makedirs(bedgraph_folder_name)

    for filename in os.listdir(wiggle_folder_name):
        if filename.endswith(".wig"):
            input_file = os.path.join(wiggle_folder_name, filename)
            output_file = os.path.join(bedgraph_folder_name, filename.replace(".wig", ".bedgraph"))

            with open(input_file, 'r') as wig_file, open(output_file, 'w') as bedgraph_file:
                chrom = None
                span = 1  # Default span

                for line in wig_file:
                    line = line.strip()

                    # Parse 'variableStep' header to get the chromosome and span
                    if line.startswith("variableStep"):
                        chrom = line.split("chrom=")[1].split()[0]
                        if "span=" in line:
                            span = int(line.split("span=")[1])
                        
                    # Process coverage data
                    elif line and chrom is not None:
                        position, coverage = line.split()
                        position = int(position)

                        # Convert 1-based wiggle position to 0-based BED format
                        bedgraph_file.write(f"{chrom}\t{position-1}\t{position-1+span}\t{coverage}\n")

            print(f"  Converted {filename} to BED graph format and saved to {output_file}")



# %% MERGING FORWARD AND REVERSE BEDGRAPH FILES FOR EACH SAMPLE
# Used by substeps 6.6 and 8.4
# -------------------------------------------------------------

# Function to read a bed graph file
def read_bedgraph(bedgraph_file):
    """
    Reads a bed graph file, skipping any 'track' or 'browser' lines.
    """
    data = []
    with open(bedgraph_file, 'r') as f:
        for line in f:
            if line.startswith('track') or line.startswith('browser'):
                continue  # Skip 'track' or 'browser' header lines
            data.append(line.strip())
    return data


# Function to concatenate multiple bed graph files
def concatenate_bedgraphs(bedgraph_files, output_bedgraph):
    """
    Concatenates multiple bed graph files (removing 'track' or 'browser' lines) 
    into a single bed graph file.
    """
    with open(output_bedgraph, 'w') as out_bed:
        for bedgraph_file in bedgraph_files:
            data = read_bedgraph(bedgraph_file)
            out_bed.write('\n'.join(data) + '\n')


# Function to find forward/reverse bed graph file pairs in the directory
def find_bedgraph_pairs(directory, pattern_forward:str, pattern_reverse:str):
    """
    Finds pairs of bed graph files in the given directory where one ends with '.forward.fiveprime.bedgraph'
    and the other ends with '.reverse.fiveprime.bedgraph'.
    Returns a list of tuples with the forward and reverse bed graph file paths.
    """
    forward_files = {}
    reverse_files = {}
    
    # Scan directory for bed graph files
    for file in os.listdir(directory):
        if file.endswith(pattern_forward):
            base_name = file.replace(pattern_forward, '')
            forward_files[base_name] = os.path.join(directory, file)
        elif file.endswith(pattern_reverse):
            base_name = file.replace(pattern_reverse, '')
            reverse_files[base_name] = os.path.join(directory, file)
    
    # Find matching pairs
    bedgraph_pairs = []
    for base_name in forward_files:
        if base_name in reverse_files:
            bedgraph_pairs.append((forward_files[base_name], reverse_files[base_name]))
    
    return bedgraph_pairs


# Function to process bed graph files, concatenate them, and save to an output directory
def concatenate_all_bedgraph_pairs(bedgraphs_folder_name, merged_bedgraphs_folder_name, pattern_forward:str, pattern_reverse:str, suffix_for_concatenated_bedgraphs:str):
    """
    Processes bed graph file pairs (forward and reverse), removes headers, concatenates them, 
    and saves to the specified output directory. Also displays which pairs are being concatenated.
    """
    print("\n Merging bedgraphs...")
    # Ensure the output directory exist
    if not os.path.exists(merged_bedgraphs_folder_name):
        os.makedirs(merged_bedgraphs_folder_name)
    
    bedgraph_pairs = find_bedgraph_pairs(bedgraphs_folder_name, pattern_forward, pattern_reverse)
    
    if not bedgraph_pairs:
        print("  No bed graph file pairs found.")
        return
    
    for forward_bed, reverse_bed in bedgraph_pairs:
        base_name = os.path.basename(forward_bed).replace(pattern_forward, '')
        
        # Define output paths
        output_bedgraph = os.path.join(merged_bedgraphs_folder_name, f'{base_name}{suffix_for_concatenated_bedgraphs}')
        
        print(f"  Concatenating files: {forward_bed} and {reverse_bed}")
        
        # Concatenate bed graph files and save to 'temp' folder
        concatenate_bedgraphs([forward_bed, reverse_bed], output_bedgraph)
        
        print(f"  Saved concatenated bed graph to: {output_bedgraph}")


# %% SORTING ALL CONCATENATED BEDGRAPH FILES
# Used by substeps 6.7 and 8.5
# ------------------------------------------

def sort_bedgraph_file(input_bed, output_bed):
    """
    Sorts a bed graph file by chromosomal position.
    """
    data = []
    
    with open(input_bed, 'r') as infile:
        for line in infile:
            line = line.strip()
            if line and not line.startswith('track') and not line.startswith('browser'):
                try:
                    # Parse position and coverage
                    parts = line.split()
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    coverage = float(parts[3])
                    data.append((chrom, start, end, coverage))
                except ValueError:
                    print(f"  Skipping line with invalid data: {line}")
        
    # Sort data by chromosomal position
    data.sort(key=lambda x: (x[0], x[1]))  # Sort by chromosome and then by start position

    # Write sorted data to output bed graph file
    with open(output_bed, 'w') as outfile:
        for chrom, start, end, coverage in data:
            outfile.write(f"{chrom}\t{start}\t{end}\t{coverage}\n")


def sort_all_bed_files(merged_bedgraphs_folder_name, sorted_merged_bedgraphs_folder_name, pattern_concatenated:str, pattern_sorted:str):
    """
    Processes all bed graph files in the input directory, sorts them by position, and saves to the output directory.
    """
    print("\n Sorting all the merged bedfiles...")
    # Ensure the output directory exists
    os.makedirs(sorted_merged_bedgraphs_folder_name, exist_ok=True)
    
    # Find all files ending with .concatenated.bedgraph in the input directory
    bed_files = glob.glob(os.path.join(merged_bedgraphs_folder_name, f"*{pattern_concatenated}"))
    
    for bed_file in bed_files:
        # Get the base name of the file and create the output file path
        base_name = os.path.basename(bed_file)
        output_bed = os.path.join(sorted_merged_bedgraphs_folder_name, base_name.replace(f"{pattern_concatenated}", f"{pattern_sorted}"))
        
        print(f"  Processing {bed_file}...")
        sort_bedgraph_file(bed_file, output_bed)
        print(f"  -> Sorted file saved to {output_bed}")



# %% REMOVE OVERLAPS IN SORTED BEDGRAPHS
# Used by substeps 6.8 and 8.6
# ---------------------------------------

def read_bedgraph_return_tuples(file_path):
    """Reads bed graph file and returns a list of tuples (chromosome, start, end, coverage)."""
    bedgraph_data = []
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            chromosome = row[0]
            start = int(row[1])
            end = int(row[2])
            coverage = float(row[3])
            bedgraph_data.append((chromosome, start, end, coverage))
    return bedgraph_data


def find_overlaps(bedgraph_data):
    """Identifies overlapping coordinates and returns non-overlapping data and overlaps."""
    non_overlapping_data = []
    overlapping_data = []
    
    # Sort data by chromosome and start coordinate for easy comparison
    bedgraph_data.sort(key=lambda x: (x[0], x[1]))
    
    prev_entry = bedgraph_data[0]
    for i in range(1, len(bedgraph_data)):
        curr_entry = bedgraph_data[i]
        
        # If the chromosome matches and there is an overlap in coordinates
        if (curr_entry[0] == prev_entry[0] and 
            curr_entry[1] < prev_entry[2]):  # Overlap condition

            overlapping_data.append(prev_entry)
            overlapping_data.append(curr_entry)
        else:
            non_overlapping_data.append(prev_entry)
        
        # The "current" becomes the "previous" in preparation for the next iteration
        prev_entry = curr_entry
    
    # Append the last entry
    non_overlapping_data.append(prev_entry)
    
    return non_overlapping_data, overlapping_data


def write_report(overlapping_data, report_path):
    """Writes the overlapping coordinates to a text report."""
    with open(report_path, 'w') as report_file:
        for entry in overlapping_data:
            report_file.write(f'{entry[0]}\t{entry[1]}\t{entry[2]}\t{entry[3]}\n')


def write_clean_bedgraph(non_overlapping_data, output_bedgraph_path):
    """Writes the cleaned bed graph data to a new file."""
    with open(output_bedgraph_path, 'w') as output_file:
        for entry in non_overlapping_data:
            output_file.write(f'{entry[0]}\t{entry[1]}\t{entry[2]}\t{entry[3]}\n')


def find_bedgraph_files(directory):
    """Finds all .bedgraph files in the specified directory."""
    return [f for f in os.listdir(directory) if f.endswith('.bedgraph')]


def clean_bedgraph_by_removing_overlaps(input_folder, output_folder, suffix_overlap_report):
    print("\n Cleaning the bedfiles by removing the overlaps...")
    # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Find all bedgraph files in the input folder
    bedgraph_files = find_bedgraph_files(input_folder)

    if not bedgraph_files:
        print(f"  No bedgraph files found in the directory: {input_folder}")
    else:
        for bedgraph_file in bedgraph_files:
            # Construct the input file path
            input_bedgraph_file = os.path.join(input_folder, bedgraph_file)
            
            # Generate the output file names by appending '_cleaned' and suffix_overlap_report
            file_base_name = os.path.splitext(bedgraph_file)[0]
            output_bedgraph_file = os.path.join(output_folder, f'{file_base_name}_cleaned.bedgraph')
            overlap_report_file = os.path.join(output_folder, f'{file_base_name}_{suffix_overlap_report}')

            # Read the bedgraph file
            bedgraph_data = read_bedgraph_return_tuples(input_bedgraph_file)

            # Find non-overlapping data and overlaps
            non_overlapping_data, overlapping_data = find_overlaps(bedgraph_data)

            # For testing my hypothesis
            # has_common = not set(non_overlapping_data).isdisjoint(overlapping_data)
            # common_elements = list(set(non_overlapping_data).intersection(overlapping_data))
            # print("-----------------------------------")
            # print("HAS COMMON ELEMENTS = ", has_common)
            # print("COMMON ELEMENTS:", common_elements)
            # input("Press Enter to continue...")

            # Write the cleaned bedgraph data and the report for overlapping coordinates
            write_clean_bedgraph(non_overlapping_data, output_bedgraph_file)
            write_report(overlapping_data, overlap_report_file)

            print(f"  Processing complete for {bedgraph_file}.")
            print(f"  Cleaned bedgraph saved to {output_bedgraph_file}.")
            print(f"  Overlapping coordinates report saved to {overlap_report_file}.")


# %% TRANSFORMING SORTED/CLEANED BEDGRAPH FILES IN BIGWIG
# Used by substeps 6.9 and 8.7
# -------------------------------------------------------

def bedgraph_to_bigwig(bedgraph_file, chrom_sizes, output_bw):
    """
    Convert a bedGraph (.bedgraph) file to bigwig (.bw) format using the UCSC bedGraphToBigWig tool.
    
    :param bedgraph_file: Path to the input bedgraph file (.bedgraph)
    :param chrom_sizes: Path to the genome chromosome sizes file
    :param output_bw: Path to the output bigwig file (.bw)
    """
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_bw)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Construct the command for bedGraphToBigWig
    command = ['bedGraphToBigWig', bedgraph_file, chrom_sizes, output_bw]
    
    try:
        # Run the command using subprocess
        subprocess.run(command, check=True)
        print(f"  Successfully converted {bedgraph_file} to {output_bw}")
    except subprocess.CalledProcessError as e:
        print(f"  Error occurred during conversion: {e}")
    except FileNotFoundError:
        print("  bedGraphToBigWig tool not found. Make sure it's installed and available in your PATH.")


def convert_all_bedgraphs_to_bigwigs(input_folder, chrom_sizes, output_folder, pattern_sorted_cleaned_suffix):
    """
    Convert all bedgraph (.bedgraph) files in a folder to bigwig (.bw) format, retaining sample names.
    
    :param input_folder: Path to the folder containing bedgraph files
    :param chrom_sizes: Path to the genome chromosome sizes file
    :param output_folder: Path to the folder where output bigwig files will be saved
    """
    print("\n Converting all bedgraphs to bigwigs...")
    # Iterate over all .bedgraph files in the input folder
    for file_name in os.listdir(input_folder):
        if file_name.endswith(pattern_sorted_cleaned_suffix):
            bedgraph_file = os.path.join(input_folder, file_name)
            # Extract the sample name (removing the pattern_sorted_cleaned_suffix)
            sample_name = file_name.replace(pattern_sorted_cleaned_suffix, "")
            output_bw = os.path.join(output_folder, f"{sample_name}.bw")
            print(f"  Converting {bedgraph_file} to {output_bw}...")
            # Call the conversion function
            bedgraph_to_bigwig(bedgraph_file, chrom_sizes, output_bw)
