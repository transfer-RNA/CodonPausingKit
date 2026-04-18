# --------------------------------------------------------
# --- CodonPausingKit                                  ---
# --- Copyright (c) 2025-2026 Aude Trinquier           ---
# --- All rights reserved.                             ---
# --- PolyForm Noncommercial License 1.0.0.            ---
# --- step1_fasta_chrom_sizes.py                       ---
# --------------------------------------------------------

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
from Bio import SeqIO
import csv


def calculate_chromosome_sizes(fasta_file, output_file):
    with open(fasta_file, "r") as handle:
        sizes = []
        for record in SeqIO.parse(handle, "fasta"):
            chromosome_name = record.id
            # Remove the NZ_ prefix if present to match chromosome naming in HRIBO files
            if chromosome_name.startswith("NZ_"):
                chromosome_name = chromosome_name[3:]
            chromosome_size = len(record.seq)
            sizes.append([chromosome_name, chromosome_size])

    with open(output_file, "w") as out_handle:
        out_handle.write("Chromosome,Size\n")
        for chrom, size in sizes:
            out_handle.write(f"{chrom},{size}\n")
    print(f"Chromosome sizes (.csv) have been saved to {output_file}")


def remove_header_from_chrom_sizes_csv(file_path, output_file):
    with open(file_path, 'r', encoding='utf-8') as infile, open(output_file, 'w', encoding='utf-8') as outfile:
        csv_reader = csv.reader(infile)
        first_row = next(csv_reader)
        if 'chrom' in first_row[0].lower() or 'name' in first_row[0].lower():
            print("Header detected and removed.")
        else:
            outfile.write("\t".join(first_row) + "\n")
        for row in csv_reader:
            outfile.write("\t".join(row) + "\n")
    print(f"Chromosomes sizes (.txt) have been saved to {output_file}")


def main():
    if len(sys.argv) != 3:
        print("Usage: step1_fasta_chrom_sizes.py <fasta_file> <txt_output>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    txt_output = sys.argv[2]
    # Construct intermediate CSV name by replacing .txt with .csv
    base, ext = os.path.splitext(txt_output)
    csv_intermediate = base + ".csv"

    # Will produce the intermediate csv file
    calculate_chromosome_sizes(fasta_file, csv_intermediate)
    # Will consume the intermediate csv file and produce its txt equivalent
    remove_header_from_chrom_sizes_csv(csv_intermediate, txt_output)


if __name__ == "__main__":
    main()
