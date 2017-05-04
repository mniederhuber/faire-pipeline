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

while getopts ":f:u:C" opt; do
	case $opt in
		f)
		rm $(snakemake --snakefile $pipePath/Snakefile --summary | tail -n+2 | cut -f1)
		snakemake --snakefile $pipePath/Snakefile --cluster-config $pipePath/clusterConfig/slurmConfig.json -R all --cluster "sbatch -J {rule} -o slurmOut/slurm-%j.out -N1 -n {cluster.threads} --time {cluster.time} --mem={cluster.mem} -A {cluster.account}" --jobs 100 
		exit
		;;
		u)
		snakemake --snakefile $pipePath/Snakefile --unlock
		exit
		;;
		C)
		# take input from file named 'config' in cwd
		snakemake --snakefile $pipePath/Snakefile --config $(cat config) --cluster-config $pipePath/clusterConfig/slurmConfig.json --cluster "sbatch -J {rule} -o slurmOut/slurm-%j.out -N1 -n {cluster.threads} --time {cluster.time} --mem={cluster.mem} -A {cluster.account}" --jobs 100
		exit
		;;
		\?)
		echo "Invalid flag: -$OPTARG" >&2
		exit
		;;
	esac
done

snakemake --snakefile $pipePath/Snakefile --cluster-config $pipePath/clusterConfig/slurmConfig.json --cluster "sbatch -J {rule} -o slurmOut/slurm-%j.out -N1 -n {cluster.threads} --time {cluster.time} --mem={cluster.mem} -A {cluster.account}" --jobs 100

