# --------------------------------------------------------
# --- CodonPausingKit                                  ---
# --- Copyright (c) 2025-2026 Aude Trinquier           ---
# --- All rights reserved.                             ---
# --- PolyForm Noncommercial License 1.0.0.            ---
# --- step3_exportAnnotatedTab.py                      ---
# --------------------------------------------------------

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas as pd


# Definition of Mini Step 3 (New 12th Dec 2025, was done manually earlier)
# (Moved to its own file on 14th Mar 2026 because this one still needs the previous conda env)
def export_annotated_tab(excel_overview_input, csv_output):
    """
    Load an Excel file with multiple tabs and export only the 'annotated' tab
    to a CSV file. Will be used on the overview.xlxs provided by HRIBO

    Parameters
    ----------
    excel_file : str
        Path to the input Excel (.xlsx) file.
    csv_output : str
        Path to the output CSV file.
    """
    # Read only the 'annotated' sheet
    df = pd.read_excel(excel_overview_input, sheet_name="annotated")

    # Save to CSV
    df.to_csv(csv_output, index=False, float_format="%.8f") # Trying with float_format to match what Excel was doing when we were exporting this tab manually

    print(f"'annotated' tab exported from {excel_overview_input} to {csv_output}")


def main():
    if len(sys.argv) != 3:
        print("Usage: step3_exportAnnotatedTab.py <excel_input> <output_csv_intermediate>")
        sys.exit(1)

    excel_input = sys.argv[1]
    output_csv_intermediate = sys.argv[2]

    # Calling substep Mini step 3
    export_annotated_tab(excel_input, output_csv_intermediate)


if __name__ == "__main__":
    main()
