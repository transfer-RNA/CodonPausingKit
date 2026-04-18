# --------------------------------------------------------
# --- CodonPausingKit                                  ---
# --- Copyright (c) 2025-2026 Aude Trinquier           ---
# --- All rights reserved.                             ---
# --- PolyForm Noncommercial License 1.0.0.            ---
# --- step6_bigwig_processing.py                       ---
# --------------------------------------------------------

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import subprocess
import csv
import pyBigWig
import glob

import tools


# %% LOOKING AT FIRST 3 LINES OF A GIVEN BIGWIG FILE
# --------------------------------------------------

# Definition of substep 6.1)
def print_header_and_first_3_lines(bigwig_file):
    print("\n Showing info of just a bigwig file (sanity test)...")
    # Open the BigWig file
    bw = pyBigWig.open(bigwig_file)
    
    # Print header information
    print("  BigWig Header Information:")
    print(f"  Chromosomes: {bw.chroms()}")
    print(f"  Total Number of Entries: {bw.header()}")

    # Get all the chromosomes in the BigWig file
    chroms = bw.chroms()

    # Initialize a counter to track the number of lines printed
    line_count = 0

    print("  Printing the first 3 lines of the bigwig file...")
    # Loop through chromosomes to extract data
    for chrom in chroms:
        chrom_length = chroms[chrom]
        # Retrieve all entries for the current chromosome
        entries = bw.intervals(chrom, 0, chrom_length)
        # Print the first 3 lines (position, coverage)
        for entry in entries:
            print(f"  {chrom}\t{entry[0]}\t{entry[1]}\t{entry[2]}")
            line_count += 1
            if line_count == 3:
                break
        if line_count == 3:
            break
    # Close the BigWig file
    bw.close()


# %% CONVERTING BIGWIG INTO WIGGLE
# --------------------------------

    # For step 6.2:
    # Now using tools


# %% SHIFTING BIGWIG TO GET A,P,E site coordinates (easier to do the shifting than on bedgraph?)
# -------------------------------------------------------------

# Utility function for 6.3
# Load the universal shift values from the CSV
def load_shift_values(csv_file):
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        # TODO: Fix that, the loop is useless
        for row in reader:
            return {
                'A': int(row['A_shift']),
                'P': int(row['P_shift']),
                'E': int(row['E_shift'])
            }


# Note: load_wiggle_files has been removed since I now construct the list of wiggle files on the fly,
# to avoid having to rely on a csv that the user has to create consistently with the files that he has.


# Utility function for 6.3
# Process each Wiggle file using the same shift values for A, P, E
def process_wiggle(wig_file, shifts, output_dir):
    print("Processing wiggle...")
    # Determine shift direction based on file name
    if wig_file.endswith(".forward.fiveprime.wig"):
        shift_direction = 1  # Positive shift
    elif wig_file.endswith(".reverse.fiveprime.wig"):
        shift_direction = -1  # Negative shift
    else:
        print(f"  Unrecognized file format for {wig_file}. Skipping...")
        return
    
    # Extract the sample name from the file name
    sample_name = os.path.splitext(os.path.basename(wig_file))[0]
    
    # Prepare output files
    output_files = {
        'A': os.path.join(output_dir, f"Asite_{sample_name}.wig"),
        'P': os.path.join(output_dir, f"Psite_{sample_name}.wig"),
        'E': os.path.join(output_dir, f"Esite_{sample_name}.wig")
    }
    
    # Create a wiggle file for each site (A, P, E)
    for site, shift in shifts.items():
        with open(output_files[site], 'w') as out_wig:
            # Write header for wiggle file
            out_wig.write(f"track type=wiggle_0 name={sample_name}_{site} description=\"Shifted {site} site\"\n")
            
            # Open and read the original wiggle file
            with open(wig_file, 'r') as infile:
                chrom = ""
                for line in infile:
                    line = line.strip()
                    if line.startswith('variableStep'):
                        # Write the variableStep line directly to the new file
                        chrom = line
                        out_wig.write(f"{chrom}\n")
                    elif line == '' or line.startswith('track'):
                        continue  # Skip empty or track lines
                    else:
                        # Parse position and value
                        pos, value = line.split()
                        pos = int(pos) + (shift * shift_direction)
                        value = float(value)
                        
                        # Write each shifted position and value to the new wiggle file
                        out_wig.write(f"{pos}\t{value}\n")


# Definition of substep 6.3
# Main function to process all Wiggle files
def process_all_wiggle_files(wig_folder_name, shift_csv, shifted_wig_folder_name):
    print("\n Processing all wiggles...")
    # Ensure the output directory exists
    if not os.path.exists(shifted_wig_folder_name):
        os.makedirs(shifted_wig_folder_name)

    # Load universal shift values
    shifts = load_shift_values(shift_csv)
    
    # Load list of Wiggle files
    wiggle_files = sorted(glob.glob(os.path.join(wig_folder_name, "*.fiveprime.wig")))

    # Process each Wiggle file
    for wig_file in wiggle_files:
        print("  Starting to process the wiggle file ", wig_file)
        process_wiggle(wig_file, shifts, shifted_wig_folder_name)



# %% GETTING ABSOLUTE COVERAGE IN THE REVERSE FILES
# -------------------------------------------------------------

    # For step 6.4
    # Now using tools


# %% TRANSFORMING WIGGLE FILES TO BEDGRAPH (going from 1-based to 0-based coordinates)
# ------------------------------------------------------------------------------------

    # For step 6.5
    # Now using tools


# %% MERGING FORWARD AND REVERSE BEDGRAPH FILES FOR EACH SAMPLE
# -------------------------------------------------------------

    # For step 6.6
    # Now using tools


# %% SORTING ALL CONCATENATED BEDGRAPH FILES
# ------------------------------------------

    # For step 6.7
    # Now using tools


# %% REMOVE OVERLAPS IN SORTED BEDGRAPHS
# ---------------------------------------

    # For step 6.8
    # Now using tools

# %% TRANSFORMING SORTED/CLEANED BEDGRAPH FILES IN BIGWIG
# -------------------------------------------------------

    # Blah=


def main():
    if len(sys.argv) != 12:
        print("Usage: step6_bigwig_processing.py <just_a_bigwig_file> <bigwig_folder_name> <wig_folder_name> <shift_csv> <shifted_wig_folder_name> <shifted_bedgraph_folder_name> <merged_bedgraph_folder_name> <sorted_merged_bedgraphs_folder_name> <cleaned_bedgraphs_no_overlap_folder_name> <chrom_sizes_txt_file> <bigwig_output_folder_name>")
        sys.exit(1)

    just_a_bigwig_file = sys.argv[1]                       # input needed for step 6.1
    bigwig_folder_name = sys.argv[2]                       # name of folder needed for step 6.2
    wig_folder_name = sys.argv[3]                          # name of folder needed for step 6.2 (in which wig files will be produced) AND for step 6.3 (where they will be found)
    shift_csv = sys.argv[4]                                # input needed for step 6.3
    shifted_wig_folder_name = sys.argv[5]                  # name of folder needed for step 6.3 (in which shifted wig files will be produced) AND for step 6.4 AND 6.5 (where they will be found)
    shifted_bedgraphs_folder_name = sys.argv[6]            # name of folder needed for step 6.5 (in which shifted bedgraph files will be produced) AND for step 6.6 (where they will be found)
    merged_bedgraphs_folder_name = sys.argv[7]             # name of folder needed for step 6.6 (in which merged bedgraph files will be produced) AND for step 6.7 (where they will be found) 
    sorted_merged_bedgraphs_folder_name = sys.argv[8]      # name of folder needed for step 6.7 (in which sorted merged bedgraphs files will be produced) AND for step 6.8 (where they will be found)
    cleaned_bedgraphs_no_overlap_folder_name = sys.argv[9] # name of folder needed for step 6.8 (in which the cleaned bedgraphs with no overlaps will be produced) AND for step 6.9 (where they will be found)
    chrom_sizes_txt_file = sys.argv[10]                    # input needed for step 6.9 (text file containing chromozomes sizes, obtained long ago in step 1)
    bigwig_output_folder_name = sys.argv[11]               # name of folder needed for step 6.9 (in which the resulting bigbig files are produced)

    # Calling substep 6.1) Extracting the first 3 lines
    # TODO: This step will probably be removed since it happens on just one bigwig file and I think it was for a test
    print_header_and_first_3_lines(just_a_bigwig_file)

    # Calling substep 6.2) Converting bigwig to wig
    tools.conversion_from_bigwig_to_wig(bigwig_folder_name, wig_folder_name)

    # Calling substep 6.3) Processing all the Wiggle files
    process_all_wiggle_files(wig_folder_name, shift_csv, shifted_wig_folder_name)

    # Calling substep 6.4) Converting (in place) reverse shifted wig values to absolute values
    tools.process_reverse_wig_files(shifted_wig_folder_name, ".reverse.fiveprime.wig")

    # Calling substep 6.5) Converting wiggle files to bedgraphs
    tools.wiggle_to_bedgraph(shifted_wig_folder_name, shifted_bedgraphs_folder_name)

    # Calling substep 6.6) Concatenating all bedgraphs pairs
    tools.concatenate_all_bedgraph_pairs(shifted_bedgraphs_folder_name, merged_bedgraphs_folder_name, ".forward.fiveprime.bedgraph", ".reverse.fiveprime.bedgraph", ".concatenated.bedgraph")

    # Calling substep 6.7) Sorting all bedgraphs
    tools.sort_all_bed_files(merged_bedgraphs_folder_name, sorted_merged_bedgraphs_folder_name, ".concatenated.bedgraph", ".sorted.bedgraph")

    # Calling substep 6.8) Cleaning bedgraphs by removing overlaps
    tools.clean_bedgraph_by_removing_overlaps(sorted_merged_bedgraphs_folder_name, cleaned_bedgraphs_no_overlap_folder_name, "overlap_report.txt")

    # Calling substep 6.9) Converting back all the produced bedgraphs to bigwigs
    tools.convert_all_bedgraphs_to_bigwigs(cleaned_bedgraphs_no_overlap_folder_name, chrom_sizes_txt_file, bigwig_output_folder_name, ".sorted_cleaned.bedgraph")



if __name__ == "__main__":
    main()