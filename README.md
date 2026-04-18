# --------------------------------------------------------
# --- CodonPausingKit                                  ---
# --- Copyright (c) 2025-2026 Aude Trinquier           ---
# --- All rights reserved.                             ---
# --- PolyForm Noncommercial License 1.0.0.            ---
# --- README.md                                        ---
# --------------------------------------------------------

**CodonPausingKit** automates the analysis of codon pausing.
More details and documentation will be added soon.

---

## 1. Installing the Conda environments

You only need to install the Conda environments once.

Please refer to the `README.md` file in the `CondaEnvs/` subfolder for full details.
In most cases, you can simply run:

```bash
cd PATH/TO/CodonPausingKit/CondaEnvs/
conda env create -n CodonPausingKit-steps1-to-3-and-10_condaEnv -f CodonPausingKit-steps1-to-3-and-10_condaEnv.yml
conda env create -n CodonPausingKit-steps4-to-5_condaEnv -f CodonPausingKit-steps4-to-5_condaEnv.yml
conda env create -n CodonPausingKit-steps6-to-9_condaEnv -f CodonPausingKit-steps6-to-9_condaEnv.yml
```

---

## 2. Setting up your input files

Place your input files in the `Inputs/` subfolder.

*(More detailed instructions on input formats and requirements will be added soon.)*

---

## 3. Running the pipeline

Once your input files are ready, run the full pipeline with:

```bash
./codon_pausing_kit.sh
```

This script is the main entry point and automates the entire workflow.

Each run of CodonPausingKit creates a new output folder with a timestamp, for example:

```
Run-20260417-180311
```

This folder contains:

* Intermediate results (in `Intermediates/`)
* Final outputs (in `Outputs/`)
* A full log of execution of this run

---

## 4. Logs and execution status

Information about the success or failure of each step is:

* Printed to the terminal during execution
* Saved to the log file: `output_run.log` inside the corresponding run folder

---

## Notes

* The pipeline automatically handles Conda environment activation/deactivation.
* Make sure all environments are successfully created before running the pipeline.
