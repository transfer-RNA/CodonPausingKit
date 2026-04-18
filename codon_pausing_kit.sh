#!/bin/bash

# ---------------------------------------------------
# --- CodonPausingKit                             ---
# --- Copyright (c) 2025-2026 Aude Trinquier      ---
# --- All rights reserved.                        ---
# --- PolyForm Noncommercial License 1.0.0.       ---
# --- codon_pausing_kit.sh                        ---
# ---------------------------------------------------

# This is the main component of the project
# Execute this wrapper bash script by doing
# ./codon_pausing_kit.sh

echo "Starting Codon Pausing Kit..."

# Bash strict mode (exit on error, treat unset vars as errors and safe pipelines)
set -euo pipefail

# Trap errors and print a helpful message
trap 'echo "❌ Error at line $LINENO: command exited with status $?"' ERR

# Detect Conda base path
CONDA_BASE="$(conda info --base)"
CONDA_SH="$CONDA_BASE/etc/profile.d/conda.sh"
if [ ! -f "$CONDA_SH" ]; then
    echo "Error: conda.sh not found at $CONDA_SH"
    exit 1
fi
source "$CONDA_SH"
# Just in case we were already in a conda environment
conda deactivate

# Get the main directory where this bash script is located (do NOT edit)
# and makes it the current directory
BASE_DIR="$(cd "$(dirname "$0")" && pwd)"

# All the python scripts are in this folder (do NOT edit)
SCRIPTS_DIR="$BASE_DIR/Scripts"

# All the inputs will be fetched from this directory (do NOT edit, instead place your files there)
INPUT_DIR="$BASE_DIR/Inputs"

# Create run folder with timestamp (do NOT edit)
TIMESTAMP="$(date +'%Y%m%d-%H%M%S')"
RUN_DIR="$BASE_DIR/Run-$TIMESTAMP"
# Inside of this timestamped run folder, there will be a subfolder for hosting the intermediate artifacts produced
INTERMEDIATE_DIR="$RUN_DIR/Intermediates"
mkdir -p "$INTERMEDIATE_DIR"

# Redirect all stdout and stderr to 'output_run.log' in the $RUN_DIR (in addition to printing on the screen)
LOG_FILE="$RUN_DIR/output_run.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "----------------------------------------------"
echo "--- Welcome to CodonPausingKit             ---"
echo "--- Copyright (c) 2025-2026 Aude Trinquier ---"
echo "--- All rights reserved.                   ---"
echo "--- PolyForm Noncommercial License 1.0.0.  ---"
echo "----------------------------------------------"
echo ""

# Create the final output folder
# Note that some files of interest can also be found in the INTERMEDIATE_DIR folder
OUTPUT_DIR="$RUN_DIR/Outputs"
mkdir -p "$OUTPUT_DIR"

### --- First Step ------------------
### ---------------------------------
echo "1st step: chromosome sizes"
echo "--------------------------"

# Define the name of the intermediate file that will contain the chromosome sizes
# This is the ouput of the first analysis (no need to customize it)
PATH_OUTPUT_CHROM_SIZES_TXT_FILE="$INTERMEDIATE_DIR/chromosome_sizes.txt"   # output of step 1, which will be needed much later for step 6.9

# Verifying the presence of a unique .fasta file in the Inputs folder
# It can been obtained from the NCBI genome
FASTA_FILES=("$INPUT_DIR"/*.fasta)
if [ ! -e "${FASTA_FILES[0]}" ]; then
    echo "❌ Error: No FASTA file found in $INPUT_DIR"
    exit 1
elif [ "${#FASTA_FILES[@]}" -ne 1 ]; then
    echo "❌ Error: Multiple FASTA files found in $INPUT_DIR. Please keep only one."
    exit 1
fi
PATH_INPUT_FASTA_FILE="${FASTA_FILES[0]}"
echo -e "FASTA file successfully found: $PATH_INPUT_FASTA_FILE \n"

echo "- Activating conda environment CodonPausingKit-steps1-to-3-and-10_condaEnv..."
conda activate CodonPausingKit-steps1-to-3-and-10_condaEnv
echo "- Running python script step1_fasta_chrom_sizes.py..."
python -u "$SCRIPTS_DIR/step1_fasta_chrom_sizes.py" "$PATH_INPUT_FASTA_FILE" "$PATH_OUTPUT_CHROM_SIZES_TXT_FILE"
echo ""
# Warning, no conda deactivate here since step 2 uses the same conda environment


### --- Second Step -----------------
### ---------------------------------
echo "2nd step: Total mapped reads"
echo "------------------------------"

# Define the name of the intermediate file that will contain the mapped (effective) reads
# This is the output of the second analysis (no need to customize it)
PATH_OUTPUT_MAPPED_READS="$INTERMEDIATE_DIR/total_mapped_reads.csv"   # output of step 2

# Verifying the presence of total_read_counts.xlsx in the Input folder
# HRIBO should have produced it. Adjust its name here if needed.
FILENAME_EXCEL_READ_COUNT="total_read_counts.xlsx"
PATH_INPUT_EXCEL_READ_COUNT_FILE="$INPUT_DIR/$FILENAME_EXCEL_READ_COUNT"
if [ ! -f "$PATH_INPUT_EXCEL_READ_COUNT_FILE" ]; then
    echo "❌ Error: Required file $FILENAME_EXCEL_READ_COUNT not found in $INPUT_DIR"
    exit 1
fi
echo -e "Excel file successfully found: $PATH_INPUT_EXCEL_READ_COUNT_FILE \n"

# Warning, no conda activate here since the previous step 1 was using the same conda environment
echo "- Running python step2_total_mapped_reads.py..."
python -u "$SCRIPTS_DIR/step2_total_mapped_reads.py" "$PATH_INPUT_EXCEL_READ_COUNT_FILE" "$PATH_OUTPUT_MAPPED_READS"
echo ""
# Warning, no conda deactivate here since step 3 uses the same conda environment


### --- Third (Mini) Step -----------------
### ---------------------------------------
echo "3rd step: Exporting the Annotated tab"
echo "--------------------------------------------"

# Verifying the presence of overview.xlsx in the Input folder
# HRIBO should have produced it. Adjust its name here if needed.
EXCEL_OVERVIEW_FILENAME="overview.xlsx"
EXCEL_OVERVIEW_FILE="$INPUT_DIR/$EXCEL_OVERVIEW_FILENAME"
if [ ! -f "$EXCEL_OVERVIEW_FILE" ]; then
    echo "❌ Error: Required file $EXCEL_OVERVIEW_FILENAME not found in $INPUT_DIR"
    exit 1
fi
echo -e "Excel file successfully found: $EXCEL_OVERVIEW_FILE \n"

PATH_OUTPUT_CSV_ANNOTATED="$INTERMEDIATE_DIR/overview.csv"

# Warning, no conda activate here since the previous step 2 was using the same conda environment
echo "- Running python step3_exportAnnotatedTab.py..."
python -u "$SCRIPTS_DIR/step3_exportAnnotatedTab.py" "$EXCEL_OVERVIEW_FILE" "$PATH_OUTPUT_CSV_ANNOTATED"
echo ""

# Now deactivating this conda env, since the step 4 uses a different one
conda deactivate


### --- Fourth Step -----------------
### ---------------------------------
echo "4th step: Reads per nucleotide and Filtering"
echo "--------------------------------------------"

# Define the name of the intermediate files that will be be produced
PATH_OUTPUT_CSV_WITH_UNIQUE_IDS="$INTERMEDIATE_DIR/overview_nd.csv"                          # output of step 4.1
PATH_OUTPUT_REPORT="$INTERMEDIATE_DIR/overview_modifications.txt"                            # output of step 4.1
PATH_OUTPUT_RPN="$INTERMEDIATE_DIR/short_overview_annotated_with_rpn.csv"                    # output of step 4.2
PATH_OUTPUT_OVERLAPPING_CDS="$INTERMEDIATE_DIR/overlapping_cds.csv"                          # output of step 4.3
PATH_OUTPUT_FILTERED_OVERLAP="$INTERMEDIATE_DIR/filtered_overlap_overview.csv"               # output of step 4.4

PATH_OUTPUT_FILTERED_LOW_RPN="$INTERMEDIATE_DIR/RPN-filteredGenesToAnalyze.csv"              # output of step 4.5
PATH_OUTPUT_FILTERED_OUT_LOW_RPN="$INTERMEDIATE_DIR/filtered_out_Lowrpn_genes.csv"           # (secondary) output of step 4.5

PATH_OUTPUT_FULLY_FILTERED="$INTERMEDIATE_DIR/GenesToAnalyze.csv"                       # output of step 4.6 ("_temp" removed in name)
PATH_OUTPUT_FILTERED_OUT_EMPTY_AMINOACID_SEQ="$INTERMEDIATE_DIR/empty_aminoacid_seq.csv"     # (secondary) output of step 4.6
PATH_OUTPUT_FILTERED_OUT_CODON_COUNT_NOT_INT="$INTERMEDIATE_DIR/non_integer_codon_count.csv" # (secondary) output of step 4.6

PATH_INPUT_MAPPED_READS="$PATH_OUTPUT_MAPPED_READS"     # This is the output from the step 2, which is now an input to this step

PATH_INPUT_CSV_ANNOTATED="$PATH_OUTPUT_CSV_ANNOTATED"   # This is the output from the mini step 3, which is now an input to this step
# IMPORTANT:
# During development, to validate the final outputs exactly with 'diff', comment the line just above
# and uncomment the next one. That will bypass the mini step 3 and instead use the overview_Run2_annotated.csv that
# we had exported manually from Excel in the old pipeline. Otherwise we get very slightly different values
# due to the fact that pandas exports with slightly different precision/rules than Excel does.
# PATH_INPUT_CSV_ANNOTATED="$INTERMEDIATE_DIR/../../REF-FOR-VALIDATION-DURING-DEV/REFCOPIED_overview_Run2_annotated.csv"

# Verifying the presence of the gff in the Input folder
# Adjust its name here if necessary
GFF_FILENAME="NCIB3610.gff"
GFF_FILE="$INPUT_DIR/$GFF_FILENAME"
if [ ! -f "$GFF_FILE" ]; then
    echo "❌ Error: Required file $GFF_FILENAME not found in $INPUT_DIR"
    exit 1
fi
echo -e "GFF file successfully found: $GFF_FILE \n"

echo "- Activating conda environment CodonPausingKit-steps4-to-5_condaEnv..."
conda activate CodonPausingKit-steps4-to-5_condaEnv
echo "- Running python script step4_rpn_and_filtering_genes_to_analyze.py..."
python -u "$SCRIPTS_DIR/step4_rpn_and_filtering_genes_to_analyze.py" "$PATH_INPUT_CSV_ANNOTATED" "$PATH_INPUT_MAPPED_READS" "$GFF_FILE" \
            "$PATH_OUTPUT_CSV_WITH_UNIQUE_IDS" "$PATH_OUTPUT_REPORT" \
            "$PATH_OUTPUT_RPN" \
            "$PATH_OUTPUT_OVERLAPPING_CDS" \
            "$PATH_OUTPUT_FILTERED_OVERLAP" \
            "$PATH_OUTPUT_FILTERED_LOW_RPN" "$PATH_OUTPUT_FILTERED_OUT_LOW_RPN" \
            "$PATH_OUTPUT_FULLY_FILTERED" "$PATH_OUTPUT_FILTERED_OUT_EMPTY_AMINOACID_SEQ" "$PATH_OUTPUT_FILTERED_OUT_CODON_COUNT_NOT_INT"
echo ""
# Warning, no conda deactivate here since step 5 uses the same conda environment


### --- Fifth Step -----------------
### ---------------------------------
echo "5th step: Retrieving codon info of genes to analyze"
echo "---------------------------------------------------"

# Define the name of the intermediate files that will be be produced
PATH_OUTPUT_CODON_INFO="$INTERMEDIATE_DIR/codon_info-GenestoAnalyze.csv"      # output of step 5.1
PATH_OUTPUT_CODONS_TO_ANALYZE="$INTERMEDIATE_DIR/codons_to_analyze.csv"       # output of step 5.2
PATH_OUTPUT_BED_CODONS_TO_ANALYZE="$INTERMEDIATE_DIR/codons_to_analyze.bed"   # output of step 5.3 

PATH_INPUT_FULLY_FILTERED="$PATH_OUTPUT_FULLY_FILTERED" # This is the output from the step 4, which is now an input to this step

echo "- Running python script step5_retrieve_codon_info.py..."
python -u "$SCRIPTS_DIR/step5_retrieve_codon_info.py" "$PATH_INPUT_FULLY_FILTERED" "$PATH_OUTPUT_CODON_INFO" "$PATH_OUTPUT_CODONS_TO_ANALYZE" "$PATH_OUTPUT_BED_CODONS_TO_ANALYZE"
echo ""

# Now deactivating this conda env, since the step 6 uses a different one
conda deactivate


### --- Sixth Step -----------------
### ---------------------------------
echo "6th step: Bigwig processing"
echo "---------------------------------------------------"

# HRIBO should have produced the bigwig files in the result output > genome-browser > coverage > fiveprime > min
# (they are "min normalized")
# The user is expected to move these bigwig files (obtained from HRIBO) in this folder:
BIGWIG_FOLDER_NAME="$INPUT_DIR/bigwig_files"

# Will be used to verify the presence of at least one specific bigwig file (.bw) in this folder
BIGWIG_FILENAME="RIBO-Dt1-1.min.reverse.fiveprime.bw"
JUST_A_BIGWIG_FILE="$BIGWIG_FOLDER_NAME/$BIGWIG_FILENAME"
if [ ! -f "$JUST_A_BIGWIG_FILE" ]; then
    echo "❌ Error: Required file $BIGWIG_FILENAME not found in $BIGWIG_FOLDER_NAME"
    exit 1
fi
echo -e "Bigwig file successfully found: $JUST_A_BIGWIG_FILE \n"

WIG_FOLDER_NAME="$INTERMEDIATE_DIR/wiggle_files"

CSV_SHIFT_FILENAME="shifts.csv"
CSV_SHIFT_FILE="$INPUT_DIR/$CSV_SHIFT_FILENAME"
if [ ! -f "$CSV_SHIFT_FILE" ]; then
    echo "❌ Error: Required file $CSV_SHIFT_FILENAME not found in $INPUT_DIR"
    exit 1
fi
echo -e "CSV of shifts successfully found: $CSV_SHIFT_FILE \n"

SHIFTED_WIG_FOLDER_NAME="$INTERMEDIATE_DIR/shifted_wiggles"
SHIFTED_BEDGRAPHS_FOLDER_NAME="$INTERMEDIATE_DIR/shifted_bedgraphs"
MERGED_BEDGRAPHS_FOLDER_NAME="$INTERMEDIATE_DIR/FR-Merged_files/RIBO/merged_bedgraphs"
SORTED_MERGED_BEDGRAPHS_FOLDER_NAME="$INTERMEDIATE_DIR/FR-Merged_files/RIBO/sorted_bedgraphs"
CLEANED_BEDGRAPHS_FOLDER_NAME="$INTERMEDIATE_DIR/FR-Merged_files/RIBO/cleaned_bedgraphs"

PATH_INPUT_CHROM_SIZES_TXT_FILE="$PATH_OUTPUT_CHROM_SIZES_TXT_FILE" # This is the output from the step 1, which is now an input to this step

BIGWIG_OUTPUT_FOLDER_NAME="$INTERMEDIATE_DIR/FR-Merged_files/RIBO/bigwig_output"


echo "- Activating conda environment CodonPausingKit-steps6-to-9_condaEnv..."
conda activate CodonPausingKit-steps6-to-9_condaEnv
echo "- Running python script step6_bigwig_processing.py..."
python -u "$SCRIPTS_DIR/step6_bigwig_processing.py" "$JUST_A_BIGWIG_FILE" "$BIGWIG_FOLDER_NAME" "$WIG_FOLDER_NAME" "$CSV_SHIFT_FILE" "$SHIFTED_WIG_FOLDER_NAME" "$SHIFTED_BEDGRAPHS_FOLDER_NAME" "$MERGED_BEDGRAPHS_FOLDER_NAME" "$SORTED_MERGED_BEDGRAPHS_FOLDER_NAME" "$CLEANED_BEDGRAPHS_FOLDER_NAME" "$PATH_INPUT_CHROM_SIZES_TXT_FILE" "$BIGWIG_OUTPUT_FOLDER_NAME"
echo ""
# Warning, no conda deactivate here since step 7 uses the same conda environment


### --- Seventh Step -----------------
### ---------------------------------
echo "7th step: Computations of codon coverage"
echo "----------------------------------------"

PATH_BIGWIG_CLEANED_FROM_STEP6_FOLDER_NAME="$BIGWIG_OUTPUT_FOLDER_NAME" # Obtained from step 6
PATH_INPUT_BED_CODONS_TO_ANALYZE="$PATH_OUTPUT_BED_CODONS_TO_ANALYZE"   # Obtained from step 5
PATH_COMPUTATIONS_FOLDER_NAME="$INTERMEDIATE_DIR/codon_computations"
PATH_RESULTS_CODON_COVERAGE_FOLDER_NAME="$INTERMEDIATE_DIR/Results-CodonCoverage"
PATH_CODON_COVERAGE_FOLDER_NAME="$INTERMEDIATE_DIR/Normalized-CodonCoverage"


echo "- Running python script step7_computation_codon_coverage.py..."
python -u "$SCRIPTS_DIR/step7_computation_codon_coverage.py" "$PATH_BIGWIG_CLEANED_FROM_STEP6_FOLDER_NAME" "$PATH_INPUT_BED_CODONS_TO_ANALYZE" "$PATH_COMPUTATIONS_FOLDER_NAME" "$PATH_RESULTS_CODON_COVERAGE_FOLDER_NAME" "$PATH_CODON_COVERAGE_FOLDER_NAME"
echo ""
# Warning, no conda deactivate here since step 8 still uses the same conda environment


### --- Eighth Step -----------------
### ---------------------------------
echo "8th step: Getting RNA coverage"
echo "-------------------------------"

# HRIBO should have produced the RNA bigwig files in (where?)
# The user is expected to move these bigwig files (obtained from HRIBO) in this folder:
RNA_BIGWIG_FOLDER_NAME="$INPUT_DIR/RNAbigwig-global"

RNA_WIG_FOLDER_NAME="$INTERMEDIATE_DIR/RNAwiggle_files"
BEDGRAPHS_FOLDER_NAME="$INTERMEDIATE_DIR/RNAbedgraph_files"
RNA_MERGED_BEDGRAPHS_FOLDER_NAME="$INTERMEDIATE_DIR/FR-Merged_files/RNA/merged_bedgraphs"
RNA_SORTED_MERGED_BEDGRAPHS_FOLDER_NAME="$INTERMEDIATE_DIR/FR-Merged_files/RNA/sorted_bedgraphs"
RNA_CLEANED_BEDGRAPHS_FOLDER_NAME="$INTERMEDIATE_DIR/FR-Merged_files/RNA/cleaned_bedgraphs"
RNA_BIGWIG_OUTPUT_FOLDER_NAME="$INTERMEDIATE_DIR/FR-Merged_files/RNA/RNAbigwig-output"


echo "- Running python script step8_getting_rna_coverage.py..."
python -u "$SCRIPTS_DIR/step8_getting_rna_coverage.py" "$RNA_BIGWIG_FOLDER_NAME" "$RNA_WIG_FOLDER_NAME" "$BEDGRAPHS_FOLDER_NAME" "$RNA_MERGED_BEDGRAPHS_FOLDER_NAME" "$RNA_SORTED_MERGED_BEDGRAPHS_FOLDER_NAME" "$RNA_CLEANED_BEDGRAPHS_FOLDER_NAME" "$PATH_INPUT_CHROM_SIZES_TXT_FILE" "$RNA_BIGWIG_OUTPUT_FOLDER_NAME"
echo ""
# Warning, no conda deactivate here since step 9 still uses the same conda environment


### --- Ninth Step -----------------
### ---------------------------------
echo "9th step: Computations of codon coverage for RNA"
echo "------------------------------------------------"

RNA_PATH_BIGWIG_CLEANED_FROM_STEP6_FOLDER_NAME="$RNA_BIGWIG_OUTPUT_FOLDER_NAME" # Obtained from step 8
PATH_INPUT_BED_CODONS_TO_ANALYZE="$PATH_OUTPUT_BED_CODONS_TO_ANALYZE"           # Obtained from step 5
RNA_PATH_COMPUTATIONS_FOLDER_NAME="$INTERMEDIATE_DIR/codon_RNAcomputations"
RNA_PATH_RESULTS_CODON_COVERAGE_FOLDER_NAME="$INTERMEDIATE_DIR/Results-CodonCoverage-RNA"
RNA_PATH_CODON_COVERAGE_FOLDER_NAME="$INTERMEDIATE_DIR/Normalized-CodonCoverage-RNA"


echo "- Running python script step9_computation_codon_coverageRNA.py..."
python -u "$SCRIPTS_DIR/step9_computation_codon_coverageRNA.py" "$RNA_PATH_BIGWIG_CLEANED_FROM_STEP6_FOLDER_NAME" "$PATH_INPUT_BED_CODONS_TO_ANALYZE" "$RNA_PATH_COMPUTATIONS_FOLDER_NAME" "$RNA_PATH_RESULTS_CODON_COVERAGE_FOLDER_NAME" "$RNA_PATH_CODON_COVERAGE_FOLDER_NAME"
echo ""
conda deactivate


### --- Tenth Step -----------------
### ---------------------------------
echo "10th step: Pre data plotting"
echo "----------------------------"


# We will copy the files that will be prepared for plotting into the Output folder
# - On the Ribo side, there's all the Asite files post normalization that have been produced
# in $PATH_CODON_COVERAGE_FOLDER_NAME (that was the final production of step 7)
# - On the RNA side, there's all the all the files produced in 
# $RNA_PATH_CODON_COVERAGE_FOLDER_NAME (that was the final production of step 9)
OUTPUTS_ASITE_PAUSING_CALCULATION_FOLDER_NAME="$OUTPUT_DIR/Asite_pausing_calculations"
mkdir -p "$OUTPUTS_ASITE_PAUSING_CALCULATION_FOLDER_NAME"
cp "$PATH_CODON_COVERAGE_FOLDER_NAME/"*_Asite_*.txt "$OUTPUTS_ASITE_PAUSING_CALCULATION_FOLDER_NAME"
cp "$RNA_PATH_CODON_COVERAGE_FOLDER_NAME/"*.txt "$OUTPUTS_ASITE_PAUSING_CALCULATION_FOLDER_NAME"


echo "- Activating conda environment CodonPausingKit-steps1-to-3-and-10_condaEnv..."
conda activate CodonPausingKit-steps1-to-3-and-10_condaEnv
echo "- Running python script step10_pre_data_plotting.py..."
python -u "$SCRIPTS_DIR/step10_stats_and_pre_data_plotting.py" "$OUTPUTS_ASITE_PAUSING_CALCULATION_FOLDER_NAME"
echo ""
conda deactivate


### --- Conclusion -----------------
### ---------------------------------
echo -e "\n ✅ All tasks completed successfully!"
echo "Intermediates, Outputs, and full log are located in: $RUN_DIR"
