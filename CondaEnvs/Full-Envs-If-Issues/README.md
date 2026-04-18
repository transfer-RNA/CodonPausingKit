# --------------------------------------------------------
# --- CodonPausingKit                                  ---
# --- Copyright (c) 2025-2026 Aude Trinquier           ---
# --- All rights reserved.                             ---
# --- PolyForm Noncommercial License 1.0.0.            ---
# --- CondaEnvs/Full-Envs-If-Issues/README.md          ---
# --------------------------------------------------------

## Full and Heavy Conda Environments (use only if the lightweight versions fail)

These Conda environments (located in `CondaEnvs/Full-Envs-If-Issues/`) include **all packages**, including dependencies that were not explicitly installed.

They were exported using:

```bash
conda env export > CodonPausingKit-steps1-to-3-and-10_condaEnv.yml
conda env export > CodonPausingKit-steps4-to-5_condaEnv.yml
conda env export > CodonPausingKit-steps6-to-9_condaEnv.yml
```

(i.e., **without** using the `--from-history` option).

These environments are heavier than the lightweight versions in `CondaEnvs/`, but they improve reproducibility and can help resolve dependency issues.

**Recommendation:**
Only use these environments if you encounter problems with the standard (lightweight) ones.

---

### If lightweight environment creation fails

First, remove any partially created or broken environments:

```bash
conda deactivate
conda env remove -n CodonPausingKit-steps1-to-3-and-10_condaEnv
conda env remove -n CodonPausingKit-steps4-to-5_condaEnv
conda env remove -n CodonPausingKit-steps6-to-9_condaEnv
```

You can list existing environments with:

```bash
conda env list
```

---

### Recreate environments using full specifications

Then recreate them using the full environment files:

```bash
cd PATH/TO/CodonPausingKit/CondaEnvs/Full-Envs-If-Issues
conda env create -n CodonPausingKit-steps1-to-3-and-10_condaEnv -f CodonPausingKit-steps1-to-3-and-10_condaEnv.yml
conda env create -n CodonPausingKit-steps4-to-5_condaEnv -f CodonPausingKit-steps4-to-5_condaEnv.yml
conda env create -n CodonPausingKit-steps6-to-9_condaEnv -f CodonPausingKit-steps6-to-9_condaEnv.yml
```

---

### Notes

* If an environment does not exist, `conda env remove` will return an error; this is safe to ignore.
