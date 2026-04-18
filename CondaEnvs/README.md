# --------------------------------------------------------
# --- CodonPausingKit                                  ---
# --- Copyright (c) 2025-2026 Aude Trinquier           ---
# --- All rights reserved.                             ---
# --- PolyForm Noncommercial License 1.0.0.            ---
# --- CondaEnvs/README.md                              ---
# --------------------------------------------------------

## Lightweight Conda Environments

The Conda environment files provided in this folder (`CondaEnvs/`) are **lightweight exports** that have been created using:

```bash
conda env export --from-history > CodonPausingKit-steps1-to-3-and-10_condaEnv.yml
conda env export --from-history > CodonPausingKit-steps4-to-5_condaEnv.yml
conda env export --from-history > CodonPausingKit-steps6-to-9_condaEnv.yml
```

(i.e., **with** the `--from-history` option).

These files include only the packages that were explicitly installed during development and **exclude automatically resolved dependencies**.

---

### Why use these environments?

* Faster to solve and install
* More portable across systems
* Easier to maintain and update

These are the **recommended default environments** for running the Codon Pausing Kit.

---

### Creating the environments

Run the following commands to create the environments:

```bash
cd PATH/TO/CodonPausingKit/CondaEnvs/
conda env create -n CodonPausingKit-steps1-to-3-and-10_condaEnv -f CodonPausingKit-steps1-to-3-and-10_condaEnv.yml
conda env create -n CodonPausingKit-steps4-to-5_condaEnv -f CodonPausingKit-steps4-to-5_condaEnv.yml
conda env create -n CodonPausingKit-steps6-to-9_condaEnv -f CodonPausingKit-steps6-to-9_condaEnv.yml
```

---

### Environment activation

You do **not** need to manually activate or deactivate the Conda environments.
The main script `codon_pausing_kit.sh` automatically handles environment switching as needed.

The environments simply need to be installed beforehand.

---

### If you encounter issues

Because these are lightweight environments, Conda resolves dependencies at install time. In some cases, this may lead to:

* Version conflicts
* Slightly different dependency versions
* Installation failures on certain systems

If you run into problems:

1. Remove the problematic environment:

   ```bash
   conda deactivate
   conda env remove -n CodonPausingKit-steps1-to-3-and-10_condaEnv
   conda env remove -n CodonPausingKit-steps4-to-5_condaEnv
   conda env remove -n CodonPausingKit-steps6-to-9_condaEnv
   ```

2. Retry creation (this can resolve transient solver issues)

---

### Still not working?

If problems persist, use the **fully specified environments** provided in:

```
CondaEnvs/Full-Envs-If-Issues/
```

These environments include all dependencies with fixed versions and are more reproducible, though heavier and slower to install.

Please refer to the README in that folder for detailed instructions.
