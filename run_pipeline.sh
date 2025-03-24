#!/bin/bash

# Prevent locally installed user site-packages from superseding conda environments
# - This environment variable will be propagated to jobs launched by snakemake
# - Reference: https://github.com/conda/conda/issues/448
export PYTHONNOUSERSITE=True

snakemake \
--snakefile Snakefile \
--configfile config.yaml \
--use-conda \
--conda-frontend conda \
--printshellcmds \
--rerun-incomplete \
-j 50 \
--cluster-config cluster.yaml \
--cluster "sbatch --parsable -c {cluster.cpus} \
-t {cluster.time} -N {cluster.nodes} \
--mem {cluster.mem} \
--output {cluster.output} \
--error {cluster.error}" \
--cluster-cancel scancel \
"$@"
