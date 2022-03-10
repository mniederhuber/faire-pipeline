#!/bin/bash

# Spencer Nystrom (snystrom@email.unc.edu)
# 2016-01-12
# Submits snakemake pipeline to slurm cluster

if [[ ! -d slurmOut ]]; then
	mkdir slurmOut
fi

snakemake --rerun-incomplete --snakefile Snakefile --use-envmodules --cluster-config clusterConfig/slurmConfig.json --latency-wait 60 --cluster "sbatch -J {rule} -o slurmOut/slurm-%j.out -N1 -n {cluster.threads} --time {cluster.time} --mem={cluster.mem} -A {cluster.account}" --jobs 100

