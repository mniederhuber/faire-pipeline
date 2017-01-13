# McKayLab FAIRE Pipeline version 3.0
## Authors: Spencer Nystrom, Chris Uyehara, Jayashree Kumar
## Date: 2016-01-12

### Description
**Pipeline Steps:**

![](docs/dag.png)


### Usage:
**Ideal Directory Structure**
```{bash}
Project_Dir
├── <genotype>-<time>-<tissue> 
│   ├── <sample>_rep1.fastq.gz
│   └── <sample>_rep2.fastq.gz
└── src
    └──	faire-pipeline/
		├── clusterConfig
		│   └── slurmConfig.json
		├── docs
		│   └── dag.png
		├── ReadMe.md
		├── slurmSubmission.sh
		├── Snakefile
		├── tester.py
		├── Tester_sub.sh
		├── zNorm.r
		└── z_norm_v2.py
  
```

1. Make project directory
1. Clone repository into src/
1. change `GenomeAssembly` if necessary (default = 'dm3')
1. Create directories for each sample
	* Copy or symlink fastq.gz files (pool technical replicate fastq.gz or do read trimming first)
1. Inside each sample directory run: ` sh ../src/faire-pipeline/slurmSubmission.sh ` 
	- To manage job submission in background run ` nohup sh ../src/faire-pipeline/slurmSubmission.sh & disown`. Not recommended.

**NOTE**: If not using submission script, you *must* provide the path to the snakefile as the first call to snakemake:
`snakemake --snakefile <path/to/Snakefile>`

### Requirements:
	- pysam (python2.7)
	- pyBigWig (python2.7)
	- snakemake (python3)
	- multiqc (python3)
	- rtracklayer (r, bioconductor)	
**Note**: 

On longleaf python packages can be installed locally with ` pip install --user <packageName> ` (use pip3 for python3 modules)

rtracklayer is preinstalled on longleaf


# ToDo:
	- Generate pooled BigWigs
	- Call range of peaks for QC analysis
		- Ideally, perhaps another pipeline for determining optimal number of peaks
			as determined by the sample with lowest ideal peak # in a set of peaks that will be compared
