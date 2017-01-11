#!/usr/bin/python3

import getpass
import os

USR=getpass.getuser()
DIR=os.getcwd()

							
REFGENEPATH='/proj/mckaylab/genomeFiles/dm3/RefGenome/dm3'	# Point directly to the refgeneome file you want to use

CTRLPATH='/proj/mckaylab/genomeFiles/dm3/ControlGenomicDNA/ControlGenomicDNA_q5_sorted_dupsRemoved_noYUHet.bed' # Point directly to negative control genomic DNA input
BLACKLIST='/proj/mckaylab/genomeFiles/dm3/dm3Blacklist.bed'

GENOMESIZE = 121400000 						# effective genome size for RPGC normalization, might need to change depending on assembly/what you blacklist
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
	print('ERROR: ' + CTRLPATH + ' does not exist. Is CTRLPATH set correctly?')
	quit()

if os.path.exists(BLACKLIST) == False:
	print('ERROR: ' + BLACKLIST + ' does not exist. Is BLACKLIST set correctly?')

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

#outdirs = ['Stats', 'Bam/', 'Sam/', 'Peakfiles/', 'BigWigs/', 'PCRdups/', 'BigWigs/ZNormalized/']
#
#for d in outdirs:
#	if os.path.isdir(d) == False:
#		os.mkdir(d)

import glob
FASTQ = glob.glob('*.fastq.gz')

def getSamples(fastqList):
	sampleList = []
	for f in fastqList:
		sample = f.split('.')[0]
		sampleList.append(sample)
	return(sampleList)

SAMPLE = getSamples(FASTQ)

# Parse info for peak calling
dirID = os.getcwd().split('/')[-1] 	# get name of directory (should be sample Pool name)
nFiles = len(SAMPLE) 			# number of files in directory should be number of biological replicates (pool technical replicats first!)
peakDir = 'Peakfiles/' 			# be sure to include '/' 

if peakDir[-1] != '/':
	# die if peakDir set incorrectly
	print('ERROR: peakDir (' + peakDir + ') does not end with /')
	quit()

peakOutName = expand("{dirID}_{nFiles}Reps_PooledPeaks", dirID = dirID, nFiles = nFiles)

rule all:
	input:
		expand("Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bam", sample = SAMPLE),
		expand("Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bam.bai", sample = SAMPLE),
		expand("Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bed", sample = SAMPLE),
		expand("BigWigs/ZNormalized/{sample}_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC_zNorm.bw", sample = SAMPLE),
		touch("Peakfiles/peakCall.done")
		#expand("{outdir}{name}{fType}", outdir = peakDir, name = peakOutName, fType = ['_peaks.narrowPeak', '_peaks.xls', '_summits.bed'])
		
rule align:
	input:
		#expand("{sample}.fastq.gz", sample = SAMPLE)
		"{sample}.fastq.gz"
	output:
		#temp(expand("Sam/{sample}.sam", sample = SAMPLE))
		temp("Sam/{sample}.sam")
	threads: 8
	params: module = bowtie2Ver, time = "1:00:00"
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
		java -Xmx8g -jar {picardPath} MarkDuplicates INPUT= {input} OUTPUT= {output.markedDups} METRICS_FILE= {output.PCRdups} REMOVE_DUPLICATES= false ASSUME_SORTED= true
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
		bam = "Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bam",
		idx = "Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bam.bai"
	output:
		"BigWigs/{sample}_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC.bw"
	threads: 8 
	params: moduleVer = deeptoolsVer, blacklist = BLACKLIST, genomeSize = GENOMESIZE
	shell:
		"""
		module purge && module load {params.moduleVer}
		bamCoverage -b {input.bam} -p {threads} --blackListFileName {params.blacklist} --normalizeTo1x {params.genomeSize} --outFileFormat bigwig --binSize 10 -e 125 -o {output}
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

if os.path.isdir == False:
	os.mkdir('Peakfiles')

rule CallPooledPeaks:
# Call Peaks using MACS. Pools all input files first.
	input:
		expand("Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bed", sample = SAMPLE)
	params:
		outdir = peakDir, 
		control = CTRLPATH,
		name = expand("{dirID}_{nFiles}Reps_PooledPeaks", dirID = dirID, nFiles = nFiles),
		moduleVer = macsVer 
	output:
		touch("Peakfiles/peakCall.done")
		#expand("{outdir}/{name}{fType}", outdir = peakDir, name = peakOutName, fType = ['_peaks.narrowPeak', '_peaks.xls', '_summits.bed'])
	shell:
		"""
		module purge && module load {params.moduleVer}
		macs2 callpeak -t {input}  -c {params.control} -n {params.name} -g dm --nomodel --extsize 125 --seed 123 --outdir {params.outdir} 
		"""

#First sort peak file by q-value in decreasing order, then cut out chromosome, start and end coordinates, peak name and q-value and take the top 5000 peaks based on q-value
#sort -n -r -k9 ${STRAIN}_q5_sorted_dupsRemoved_noYUHet_peaks.narrowPeak | cut -f1,2,3,4,9 | head -5000 > ${STRAIN}_q5_sorted_dupsRemoved_noYUHet_Top5000Peaks.bed
	

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

