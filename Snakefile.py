#!/usr/bin/python3

import getpass
import os


USR=getpass.getuser()
DIR=os.getcwd()

							
REFGENEPATH='/proj/mckaylab/genomeFiles/dm3/RefGenome/dm3'	# Point directly to the refgeneome file you want to use

CTRLPATH='/proj/mckaylab/genomeFiles/dm3/ControlGenomicDNA/ControlGenomicDNA_q5_sorted_dupsRemoved_noYUHet.bed' # Point directly to negative control genomic DNA input

#PIPEPATH=str(DIR + '/faire-pipeline')
PIPEPATH=str('/pine/scr/s/n/snystrom/Bitbucket/faire-pipeline')

stdOUT=str(DIR + '/OutputFiles/')				# standard output directory, end path with '/'
stdERR=str(DIR + '/ErrorFiles/')				# standard error directory, end path with '/'

# SLURM Params:
#numThreads=8
#runTime=1:00:00
#maxMem=16G
#group=rc_dmckay1_pi

##############################
# Module Versions:

bowtie2Ver = str('bowtie2/2.2.8')
samtoolsVer = str('samtools/1.3.1')
bedtoolsVer = str('bedtools/2.25.0')

picardVer = str('2.2.4')
picardPath = str('/nas02/apps/picard-' + picardVer + '/picard-tools-' + picardVer + '/picard.jar')

deeptoolsVer = str('deeptools/2.4.1')
macsVer = str('macs/2016-02-15')
ucscVer = str('ucsctools/320')
rVer = str('r/3.3.1')

##############################


# Parse commandline flags
#STRAIN=${1%%.*}
#ALIGN=$2
#PEAK=$3


# This series of if statements checks that the control .bed file defined by $CTRLPATH exists 
# and creates standard out/error and Stats directories if they don't already exist
# Also checks whether stdOUT, stdERR, and NETSCR end with '/', exits with errorcode 1 if so

if os.path.exists(CTRLPATH) == False:
	print(CTRLPATH + ' does not exist. Is CTRLPATH set correctly?')
	quit()

if stdOUT[-1] != '/':
	print('ERROR: stdOUT (' + stdOUT + ') must end in /')
	quit()

if stdERR[-1] != '/':
	print('ERROR: stdERR (' + stdERR + ') must end in /')
	quit()

if os.path.isdir(stdOUT) == False:
	os.mkdir(stdOUT)

if os.path.isdir(stdERR) == False:
	os.mkdir(stdERR)

outdirs = ['Stats', 'Bam/', 'Sam/', 'Peakfiles/', 'BigWigs/', 'PCRdups/', 'BigWigs/ZNormalized/']

for d in outdirs:
	if os.path.isdir(d) == False:
		os.mkdir(d)

import glob
FASTQ = glob.glob('*.fastq.gz')
SAMPLE = FASTQ[0].split('.')[0]

rule all:
	input:
		expand("Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bam", sample = SAMPLE),
		expand("Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bam.bai", sample = SAMPLE),
		expand("Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bed", sample = SAMPLE),
		expand("BigWigs/ZNormalized/{sample}_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC_zNorm.bw", sample = SAMPLE)
		
rule align:
	input:
		expand("{sample}.fastq.gz", sample = SAMPLE)
	output:
		temp(expand("Sam/{sample}.sam", sample = SAMPLE))
	threads: 8
	params: module = bowtie2Ver
	shell:
		"""
		module purge && module load {params.module}
		bowtie2 --seed 123 -x {REFGENEPATH} -p {threads} -U {input} -S {output}
		"""

rule sam2bam:
	input:
		"Sam/{sample}.sam"
	output:
		"Bam/{sample}.bam"	
	params: moduleVer = samtoolsVer 
	shell:
		"""
		module purge && module load {params.moduleVer}
		samtools view -@ 4 -b {input} > {output}
		"""
rule qFilter:
	input:
		"Bam/{sample}.bam"
	output:
		temp("Bam/{sample}_q5.bam")
	params: moduleVer = samtoolsVer 
	shell:
		"""
		module purge && module load {params.moduleVer}
		samtools view -@ 4 -bq 5 {input} > {output}
		"""
rule bamSort:
	input:
		"Bam/{sample}_q5.bam"
	output:
		temp("Bam/{sample}_q5_sorted.bam")
	params: moduleVer = samtoolsVer 
	shell:
		"""
		module purge && module load {params.moduleVer}
		samtools sort -@ 4 -o {output} {input}
		"""
rule markDups:
	input:
		"Bam/{sample}_q5_sorted.bam"

	output:
		markedDups = temp("Bam/{sample}_q5_sorted_dupsMarked.bam"),
		PCRdups = "PCRdups/{sample}_PCR_duplicates"
	params: moduleVer = str('picard/' + picardVer)

	shell:
		"""
		module purge && module load {params.moduleVer}
		java -Xmx4g -jar {picardPath} MarkDuplicates INPUT= {input} OUTPUT= {output.markedDups} METRICS_FILE= {output.PCRdups} REMOVE_DUPLICATES= false ASSUME_SORTED= true
		"""
# continue making rules here

rule removeDups:
# Remove the reads that are marked as pcr duplicates from the bam file using the bit flag for pcr dups
	input:
		"Bam/{sample}_q5_sorted_dupsMarked.bam"
	output:
		temp("Bam/{sample}_q5_sorted_dupsRemoved.bam")
	params: moduleVer = samtoolsVer

	shell:
		"""
		module purge && module load {params.moduleVer}
		samtools view -@ 4 -bF 0x400 {input} > {output} 
		"""

rule idxNoDups:
# Create index for bam file. This is needed for the next step to remove Chrm Y, U and Het
	input:
		"Bam/{sample}_q5_sorted_dupsRemoved.bam"
	output:
		temp("Bam/{sample}_q5_sorted_dupsRemoved.bam.bai")
	params: moduleVer = samtoolsVer
	shell:
		"""
		module purge && module load {params.moduleVer}
		samtools index {input}
		"""
rule noYUHet:
# Extract only reads from Chrm 2R,2L,3R,3L,4 and X
	input:
		bam = "Bam/{sample}_q5_sorted_dupsRemoved.bam",
		idx = "Bam/{sample}_q5_sorted_dupsRemoved.bam.bai"
	output:
		"Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bam"
	params: moduleVer = samtoolsVer
	shell:
		"""
		module purge && module load {params.moduleVer}
		samtools view -@ 4 -b {input.bam} chr2L chr2R chr3L chr3R chr4 chrX > {output} 
		"""

rule noYUHet_idx:
# Create new index for filtered bam file. Needed to convert bam file to bed file
	input:
		"Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bam"
	output:
		"Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bam.bai"
	params: moduleVer = samtoolsVer
	shell:
		"""
		module purge && module load {params.moduleVer}
		samtools index {input} 
		"""
rule noYUHet_toBed:
# Convert the bam file into a bed file
	input:
		bam = "Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bam",
		idx = "Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bam.bai"
	output:
		"Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bed"
	params: moduleVer = bedtoolsVer
	shell:
		"""
		module purge && module load {params.moduleVer}
		bedtools bamtobed -i {input.bam} > {output}
		"""

rule bigWig:
# Bam Coverage to output bigwig file normalized to genomic coverage
	input:
		"Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bam"
	output:
		"BigWigs/{sample}_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC.bw"
	threads: 8
	params: moduleVer = deeptoolsVer
	shell:
		"""
		module purge && module load {params.moduleVer}
		bamCoverage -b {input} -p {threads} --normalizeTo1x 121400000 --outFileFormat bigwig --binSize 10 -e 125 -o {output}
		"""

rule zNormBigWig:
# Z-Normalize Bigwig Files
	input:
		"BigWigs/{sample}_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC.bw"
	output:
		"BigWigs/ZNormalized/{sample}_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC_zNorm.bw"
	params: pipePath = PIPEPATH, moduleVer = rVer
	shell:
		"""
		module purge && module load {params.moduleVer}
		Rscript --vanilla {params.pipePath}/zNorm.r {input} {output}
		"""

	

## Create collected flagstats files
#srun for BAM in \$(ls Bam/${sample}*.bam); do
#	NAME=\${BAM##*/}
#	echo "\${NAME%%.*}: " >> Stats/${sample}_flagstats.txt
#	samtools flagstat \${BAM} | grep -v '^0 + 0' >> Stats/${sample}_flagstats.txt;
#	echo >> Stats/${sample}_flagstats.txt
#done
#
## Parse FlagStats
#srun python2.7 ${PIPEPATH}/parseflagstat.py Stats/

# Purge all currently loaded modules, load required modules for pipeline:

#module purge

