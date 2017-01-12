# McKayLab FAIRE Pipeline version 3.0
## Authors: Spencer Nystrom, Chris Uyehara, Jayashree Kumar
## Date: 2016-01-12

### Description
#### Pipeline Steps:
![]('docs/dag.png')


### Usage:



### Requirements:
	- pysam (python2.7)
	- pyBigWig (python2.7)
	- snakemake (python3)
	- multiqc (python3)
	- rtracklayer (r, bioconductor)	
**Note**: 

On longleaf python pacakges can be installed locally with ` pip install --user <pacakgeName> ` (use pip3 for python3 modules)

rtracklayer is preinstalled on longleaf


# ToDo:
	- Generate pooled BigWigs
	- Call range of peaks for QC analysis
		- Ideally, perhaps another pipeline for determining optimal number of peaks
			as determined by the sample with lowest ideal peak # in a set of peaks that will be compared
