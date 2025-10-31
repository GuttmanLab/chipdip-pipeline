# History

Major releases
- [v1](https://github.com/GuttmanLab/chipdip-pipeline/releases/tag/v1.0_publication), 2024-10: release for Nature Genetics paper
  - Snakemake v7.31
- [v2](https://github.com/GuttmanLab/chipdip-pipeline/releases/tag/v2.0.0), 2025-05: Update to Snakemake 8+; re-structure repository to a [Snakemake-standardized workflow](https://snakemake.github.io/snakemake-workflow-catalog/docs/catalog.html#standardized-usage-workflows)
- [v3](https://github.com/GuttmanLab/chipdip-pipeline/releases/tag/v3.0.0), 2025-10: Use BAM files with custom tags instead of text-based cluster files
  - Pipeline implements streaming approach to demultiplexing so that only 1 cluster needs to be fully loaded in memory at a time.

# Complexity analysis

Variables
- $k$: number of antibody targets
- $m$: number of oligo reads
- $n$: number of chromatin reads
- $s$: number of splits (parallelized chunks)
- $c$: number of clusters

Pipeline v1, v2
- Time complexity: overall $m+n$
	- Oligo processing
		- Trimming: $m/s$
		- FASTQ to BAM: $m/s$
		- Sorting by coordinate (UMI): $(m/s) \log (m/s)$
		- Merge across splits: $m$
	- Chromatin processing
		- Trimming: $n/s$
		- Alignment: $n/s$
		- Renaming, masking: $n/s$
		- Merge across splits: $n$
	- Clustering
		- Make clusters: $(m+n)/s$
		- Sort clusters: $c \log c$
		- Merge clusters: $m+n$
	- Demultiplexing: $m+n$
- Memory complexity: $m+n$ (from storing the entire cluster file in memory!)

Pipeline v3
- Time complexity: overall $m+n$
	- Oligo processing
		- Trimming: $m/s$
		- FASTQ to BAM: $m/s$
		- Sort by barcode: $(m/s) \log (m/s)$
		- Merge across splits and deduplicate: $m$
	- Chromatin processing
		- Trimming: $n/s$
		- Alignment: $n/s$
		- Renaming, masking: $n/s$
		- Extract barcode to tags: $n/s$
		- Sort by barcode: $(n/s) \log (n/s)$
		- Merge across splits and deduplicate: $n$
	- Demultiplexing: $m+n$
- Memory complexity: $c$ (only store 1 entire cluster in memory at a time)

# Benchmarks

## Environments

chipdip.yaml

| pipeline version | Operating System | conda platform | Disk space |
| ---------------- | ---------------- | -------------- | ---------- |
| v3 | Red Hat Enterprise Linux 9.5 | linux-64 | 2.1 GB |
| v3 | Ubuntu 24.04.2 LTS | linux-64 | 2.2 GB |

How disk space was checked:

```
conda env create -f workflow/envs/chipdip.yaml -p /tmp/chipdip
du -sh /tmp/chipdip
```

## Resource utilization

| Dataset | Pipeline version | CPU cores and architecture | Operating System | Pipeline wall time (HH:MM:SS) | Core-hours utilized (HH:MM:SS) | Maximum RAM utilization | Disk space of output files | Benchmark times include creating `chipdip` conda environment? | Notes |
| ------- | ---------------- | -------------------------- | ---------------- | ----------------------------- | ------------------------------ | ----------------------- | --------------------- | ------------------------------------------------------------- | ----- |
| example dataset provided in `data/`: 4,532 paired reads (≤ 130 bp read 1, ≤ 150 bp read 2) | v3 | 2 cores of AMD EPYC 7763 | Ubuntu 24.04.2 LTS | 00:04:06 | 00:04:45 | 2.91 GB | 12 M | No | run "locally" on a single node on [GitHub Codespaces](https://github.com/features/codespaces) with `--cores 2`; profiled using GNU time |
| example dataset provided in `data/` | v3 | 8 cores of Intel Xeon Gold 6130 | RHEL 9.5 | 00:02:21 | 00:04:27 | 2.76 GB | 13 MB | No | run "locally" on a single node on [Caltech's HPC cluster](https://www.hpc.caltech.edu/) with `--cores 8`; profiled using SLURM |
| example dataset provided in `data/` | v1 | 2 cores of AMD EPYC 7763 | Ubuntu 20.04.6 LTS | 00:13:01 | 00:20:28 | 2.91 GB | 20 MB | No | run "locally" on a single node on [GitHub Codespaces](https://github.com/features/codespaces) |
| example dataset provided in `data/` | v1 | Intel Core i5-4300U | Ubuntu 22.04.4 LTS via WSL | 00:12:23 | 00:41:31 | 2.76 GB | 14 MB | No | Hyperthreading enabled; `snakemake` run with `--jobs 4` |
| example dataset provided in `data/` | v1 | Apple M2 (4 performance + 4 efficiency cores) | macOS Sequoia 15.1.1 | 00:04:45 | 00:14:06 | 2.33 GB | 16 MB | No | conda environment targeting `osx-64` platform; `snakemake` run with `--jobs 4` |
| example dataset provided in `data/` | v1 | 2 cores of Intel Xeon Gold 6130 | [RHEL9.3](https://docs.redhat.com/en/documentation/red_hat_enterprise_linux/9/html/9.3_release_notes/index) | 00:21:20 | 00:26:35 | 827.49 MB | 20 MB | Yes | run "locally" on a single node on [Caltech's HPC cluster](https://www.hpc.caltech.edu/) |
| 49,222,185 paired reads (89 bp read 1, 209 bp read 2) | v1 | 8 cores of Intel Xeon Gold 6130 | [RHEL9.3](https://docs.redhat.com/en/documentation/red_hat_enterprise_linux/9/html/9.3_release_notes/index) | 02:43:38 | 09:49:34 | 10.65 GB | 45 GB | Yes | run "locally" on a single node on [Caltech's HPC cluster](https://www.hpc.caltech.edu/) |

Using [GNU time](https://www.gnu.org/software/time/):
```
time -v -o time.log ./run_pipeline.sh
```

Using SLURM:
```
jobID=$(sbatch --parsable -c <n_cpus> --mem=20G --time=<time> run_pipeline.sh)
# wait for job to finish
sacct -j $jobid$ --format=JobID,Elapsed,AllocCPUS,TotalCPU,ReqMem,MaxRSS
```
