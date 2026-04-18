#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# --------------------------------------------------------
# --- CodonPausingKit                                  ---
# --- Copyright (c) 2025-2026 Aude Trinquier           ---
# --- All rights reserved.                             ---
# --- PolyForm Noncommercial License 1.0.0.            ---
# --- step8_getting_rna_coverage.py                    ---
# --------------------------------------------------------

import sys
import os
import subprocess
import csv

import tools # Using some common functionality shared between RIBO and RNA sides


# Note: this step 8 is the equivalent of step 6, but this time for RNA (instead of RIBO)


# %% CONVERTING BIGWIG INTO WIGGLE
# --------------------------------

    # For step 8.1:
    # Using the common functionality provided in tools.py


# %% GETTING ABSOLUTE COVERAGE IN THE RNA REVERSE FILES
# -----------------------------------------------------

    # For step 8.2:
    # Using the common functionality provided in tools.py


# %% TRANSFORMING WIGGLE FILES TO BEDGRAPH (going from 1-based to 0-based coordinates)
# -----------------------------------------------------

    # For step 8.3:
    # Using the common functionality provided in tools.py


# %% MERGING FORWARD AND REVERSE RNA BEDGRAPH FILES FOR EACH SAMPLE
# -----------------------------------------------------------------

    # For step 8.4:
    # Using the common functionality provided in tools.py


# %% SORTING ALL CONCATENATED RNA BEDGRAPH FILES
# ----------------------------------------------

    # For step 8.5:
    # Using the common functionality provided in tools.py


# %% REMOVE OVERLAPS IN SORTED RNA BEDGRAPHS
# ------------------------------------------

    # For step 8.6:
    # Using the common functionality provided in tools.py


# %% TRANSFORMING SORTED/CLEANED RNA BEDGRAPH FILES IN BIGWIG
# -----------------------------------------------------------

    # For step 8.7:
    # Using the common functionality provided in tools.py


def main():
    if len(sys.argv) != 9:
        print("Usage: step8_getting_rna_cobverage.py <RNAbigwig_folder_name> <RNAwig_folder_name> <RNAbedgraphs_folder_name> <RNAmerged_bedgraphs_folder_name> <RNAsorted_merged_bedgraphs_folder_name> <RNAcleaned_bedgraphs_no_overlap_folder_name> <chrom_sizes_txt_file> <RNAbigwig_output_folder_name>")
        sys.exit(1)

    RNAbigwig_folder_name = sys.argv[1]                      # name of folder needed for step 8.1
    RNAwig_folder_name = sys.argv[2]                         # name of folder needed for step 8.2 (in which wig files will be produced) AND for step 8.3 (where they will be found)
    RNAbedgraphs_folder_name = sys.argv[3]                   # name of folder needed for step 8.3 (in which bedgraphs files will be produced) AND for step 8.4
    RNAmerged_bedgraphs_folder_name = sys.argv[4]            # name of folder needed for step 8.4 (in which merged bedgraph files will be produced) AND for step 8.5 (where they will be found) 
    RNAsorted_merged_bedgraphs_folder_name = sys.argv[5]     # name of folder needed for step 8.5 (in which sorted bedgraphs files will be produces) AND for step 8.6 (where they will be found)
    RNAcleaned_bedgraphs_no_overlap_folder_name = sys.argv[6]# name of folder needed for step 8.6 (in which cleaned bedgraphs files will be produced) AND for step 8.7 (where they will be found)
    chrom_sizes_txt_file = sys.argv[7]                       # input needed for step 8.7 (text file containing chromozomes sizes, obtained long ago in step 1)
    RNAbigwig_output_folder_name = sys.argv[8]               # name of folder needed for step 8.7 (in which the resulting bigbig files are produced)

    # Calling substep 8.1) Converting bigwig to wig
    tools.conversion_from_bigwig_to_wig(RNAbigwig_folder_name, RNAwig_folder_name)

    # Calling substep 8.2) Converting (in place) reverse shifted wig values to absolute values
    tools.process_reverse_wig_files(RNAwig_folder_name, ".reverse.global.wig")

    # Calling substep 8.3) Converting wiggle files to bedgraphs
    tools.wiggle_to_bedgraph(RNAwig_folder_name, RNAbedgraphs_folder_name)

    # Calling substep 8.4) Concatenating all bedgraphs pairs
    tools.concatenate_all_bedgraph_pairs(RNAbedgraphs_folder_name, RNAmerged_bedgraphs_folder_name, ".forward.global.bedgraph", ".reverse.global.bedgraph", ".concatenatedRNA.bedgraph")

    # Calling substep 8.5) Sorting all bedgraphs
    tools.sort_all_bed_files(RNAmerged_bedgraphs_folder_name, RNAsorted_merged_bedgraphs_folder_name, ".concatenatedRNA.bedgraph", ".sortedRNA.bedgraph")

    # Calling substep 8.6) Cleaning bedgraphs by removing overlaps
    tools.clean_bedgraph_by_removing_overlaps(RNAsorted_merged_bedgraphs_folder_name, RNAcleaned_bedgraphs_no_overlap_folder_name, "overlap_RNAreport.txt")

    # Calling substep 8.7) Converting back all the produced bedgraphs to bigwigs
    tools.convert_all_bedgraphs_to_bigwigs(RNAcleaned_bedgraphs_no_overlap_folder_name, chrom_sizes_txt_file, RNAbigwig_output_folder_name, ".sortedRNA_cleaned.bedgraph")
 

if __name__ == "__main__":
    main()
