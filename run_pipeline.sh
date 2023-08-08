#!/bin/bash

snakemake \
--snakefile Snakefile \
--configfile config.yaml \
--use-conda \
--conda-frontend conda \
--printshellcmds \
--rerun-incomplete \
-j 50 \
--cluster-config cluster.yaml \
--cluster "sbatch -c {cluster.cpus} \
-t {cluster.time} -N {cluster.nodes} \
--mem {cluster.mem} \
--output {cluster.output} \
--error {cluster.error}" \
--cluster-cancel scancel \
"$@"
