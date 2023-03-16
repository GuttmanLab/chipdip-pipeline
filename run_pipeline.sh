#!/bin/bash

snakemake \
--snakefile Snakefile \
--use-conda \
--rerun-incomplete \
-j 50 \
--cluster-config cluster.yaml \
--configfile config.yaml \
--cluster "sbatch -c {cluster.cpus} \
-t {cluster.time} -N {cluster.nodes} \
--mem {cluster.mem} \
--output {cluster.output} \
--error {cluster.error}"
