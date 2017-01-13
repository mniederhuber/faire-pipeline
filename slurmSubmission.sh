#!/bin/bash

# Spencer Nystrom (snystrom@email.unc.edu)
# 2016-01-12
# Submits snakemake pipeline to slurm cluster

pipePath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
	# This grabs the full path the submission script lies in to detect where the faire pipeline actually is on the filesystem
	# makes it easier to tell snakemake where the Snakefile and config.json files actually are
	# (from http://stackoverflow.com/questions/59895/getting-the-source-directory-of-a-bash-script-from-within)


if [[ ! -d slurmOut ]]; then
	mkdir slurmOut
fi

snakemake --snakefile $pipePath/Snakefile --cluster-config $pipePath/clusterConfig/slurmConfig.json --cluster "sbatch -J {rule} -o slurmOut/slurm-%j.out -N1 -n {cluster.threads} --time {cluster.time} --mem={cluster.mem} -A {cluster.account}" --jobs 100
