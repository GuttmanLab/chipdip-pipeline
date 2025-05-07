#!/bin/bash

# Prevent locally installed user site-packages from superseding conda environments
# - This environment variable will be propagated to jobs launched by snakemake
# - Reference: https://github.com/conda/conda/issues/448
export PYTHONNOUSERSITE=True

# add --workflow-profile workflow/profiles/slurm for SLURM cluster job submission
snakemake "$@"
