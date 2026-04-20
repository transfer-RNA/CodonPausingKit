=================================
**Codon Pausing Kit**
=================================

**Copyright © 2025–2026 Aude Trinquier**  
All rights reserved.

**License:** PolyForm Noncommercial 1.0.0

---

## README.md

**CodonPausingKit** performs genome-wide codon pausing from *Bacillus subtilis* Ribo-seq data.

One prerequisite for our Codon Pausing Kit workflow is that the Ribo-seq data is processed using the HRIBO workflow (https://github.com/RickGelhausen/HRIBO).
Several files from your HRIBO "Output" folder will serve as inputs for your Codon Pausing Kit run.

Our **CodonPausingKit** is optimized for *Bacillus subtilis* Ribo-seq data as it is using custom offsets for 5´-end position to gain nucleotide resolution.
The offsets used are listed in `shifts.csv` provided in the `Inputs/` subfolder and have been optimized for *Bacillus subtilis* data.
More information about determination of offset values can be found in the related publication (doi link pending).
(Note that this method diverges from the 3´-end assignment used for *Escherichia coli* Ribo-seq data).

More information about input files required for a **CodonPausingKit** can be found in section 2 below.

We provide our inputs as examples on our Mendeley data repository (Add link when published) and for reproducibility of our analysis.


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

One prerequisite for our Codon Pausing Kit workflow is that the Ribo-seq data is processed using the HRIBO workflow (https://github.com/RickGelhausen/HRIBO).
Several files from your HRIBO "Output" folder will serve as inputs for your Codon Pausing Kit run.

Place your input files in the `Inputs/` subfolder:
The folder already contains `shifts.csv`: the .csv files with 5´-end offset values optimized for getting A, P & E-sites positions from *Bacillus subtilis* Ribo-seq data.
You will need to add to the `Inputs/` subfolder:
    - The *Bacillus subtilis* genome reference files:
          - `3610.fasta`: genome fasta sequence downloaded from NCBI.
          - `NCIB3610.gff`: genome annotation GFF downloaded from NCBI.
    - The `overview.xlsx` file obtained from the HRIBO workflow (in the HRIBO `output` folder).
    - The `total_read_counts.xlsx` file obtained from the HRIBO workflow (in the `quality-control` subfolder the `output` folder). 
    - Bigwig files from the HRIBO workflow:
          - RNA coverage files 
          - Ribosome coverage files

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
