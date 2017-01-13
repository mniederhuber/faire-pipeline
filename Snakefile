#!/usr/bin/python3

import os
import glob

GenomeAssembly = 'dm3'

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

python3Ver = str('python/3.5.1')

##############################

# Dictionary of effective genome sizes. 
# Used as lookup table later for RPGC adjustments 
# based on GenomeAssembly
effectiveGenomeSizes = {'dm3':121400000}

# Symlink to pipeline path for calling scripts with commandline arguments (zNorm.r for example)
# Requires that snakemake call be structured as `snakemake --snakefile <path/to/snakefile>`
PIPEPATH = os.path.dirname(os.path.abspath(sys.argv[2]))
symLink = str(os.getcwd() + '/.faire-pipeline')

if os.path.isdir(symLink) == False:
	os.symlink(PIPEPATH, symLink) 

if os.path.exists(symLink + '/Snakefile') == False:
	print('ERROR: no Snakefile in .faire-pipeline. Ensure --snakefile <path/to/snakefile> is 3rd argument of snakemake call.')
	quit()


							
REFGENEPATH=str('/proj/mckaylab/genomeFiles/' + GenomeAssembly + '/RefGenome/' + GenomeAssembly)	# Point directly to the refgeneome file you want to use

# Point directly to negative control genomic DNA input
CTRLPATH=str('/proj/mckaylab/genomeFiles/' + GenomeAssembly + '/ControlGenomicDNA/ControlGenomicDNA_q5_sorted_dupsRemoved_noYUHet.bed') 

BLACKLIST=str('/proj/mckaylab/genomeFiles/'+ GenomeAssembly + '/' + GenomeAssembly + 'Blacklist.bed')

# effective genome size for RPGC normalization, might need to change depending on assembly/what you blacklist
GENOMESIZE = effectiveGenomeSizes[GenomeAssembly]

if os.path.exists(glob.glob(REFGENEPATH + '*')[0]) == False:
	print('ERROR: REFGENEPATH (' + REFGENEPATH + ') does not contain bowtie genome files.')
	quit()

if os.path.exists(BLACKLIST) == False:
	print('ERROR: BLACKLIST (' + BLACKLIST + ') does not exist.')
	quit()

if os.path.exists(CTRLPATH) == False:
	print('ERROR: ' + CTRLPATH + ' does not exist. Is CTRLPATH set correctly?')
	quit()


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
		#expand("Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bam", sample = SAMPLE),
		#expand("Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bam.bai", sample = SAMPLE),
		#expand("Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bed", sample = SAMPLE),
		expand("BigWigs/ZNormalized/{sample}_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC_zNorm.bw", sample = SAMPLE),
		"Peakfiles/.peakCall.done",
		"multiqc_report.html"
	
rule align:
	input:
		#expand("{sample}.fastq.gz", sample = SAMPLE)
		"{sample}.fastq.gz"
	output:
		#temp(expand("Sam/{sample}.sam", sample = SAMPLE))
		sam = temp("Sam/{sample}.sam"),
		logInfo = "logs/{sample}_bowtie2.txt"
	threads: 8
	params: module = bowtie2Ver
	shell:
		"""
		module purge && module load {params.module}
		bowtie2 --seed 123 -x {REFGENEPATH} -p {threads} -U {input} -S {output.sam} 2> {output.logInfo}
		"""

rule sam2bam:
	input:
		"Sam/{sample}.sam"
	output:
		bam = "Bam/{sample}.bam",
		flagstat = "logs/{sample}.flagstat"
	params: moduleVer = samtoolsVer 
	shell:
		"""
		module purge && module load {params.moduleVer}
		samtools view -@ 4 -b {input} > {output.bam} &&
		samtools flagstat {output.bam} > {output.flagstat}
		"""
rule qFilter:
	input:
		"Bam/{sample}.bam"
	output:
		bam = temp("Bam/{sample}_q5.bam"),
		flagstat = "logs/{sample}_q5.flagstat"
	params: moduleVer = samtoolsVer 
	shell:
		"""
		module purge && module load {params.moduleVer}
		samtools view -@ 4 -bq 5 {input} > {output.bam} &&
		samtools flagstat {output.bam} > {output.flagstat}
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
# Create index for bam file. This is needed for the next step to remove Chrm Y, U and Het
	input:
		"Bam/{sample}_q5_sorted_dupsMarked.bam"
	output:
		bam = temp("Bam/{sample}_q5_sorted_dupsRemoved.bam"),
		flagstat = "logs/{sample}_q5_sorted_dupsRemoved.flagstat",
		idx = temp("Bam/{sample}_q5_sorted_dupsRemoved.bam.bai")
	params: moduleVer = samtoolsVer

	shell:
		"""
		module purge && module load {params.moduleVer}
		samtools view -@ 4 -bF 0x400 {input} > {output.bam} &&
		samtools flagstat {output.bam} > {output.flagstat} &&
		samtools index {output.bam}
		"""

rule noYUHet:
# Extract only reads from Chrm 2R,2L,3R,3L,4 and X
# Create new index for filtered bam file. Needed to convert bam file to bed file
	input:
		bam = "Bam/{sample}_q5_sorted_dupsRemoved.bam",
		idx = "Bam/{sample}_q5_sorted_dupsRemoved.bam.bai"
	output:
		bam = "Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bam",
		flagstat =  "logs/{sample}_q5_sorted_dupsRemoved_noYUHet.flagstat",
		idx = "Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bam.bai"
	params: moduleVer = samtoolsVer
	shell:
		"""
		module purge && module load {params.moduleVer}
		samtools view -@ 4 -b {input.bam} chr2L chr2R chr3L chr3R chr4 chrX > {output.bam} &&
		samtools flagstat {output.bam} > {output.flagstat} &&
		samtools index {output.bam}
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
		bw = "BigWigs/{sample}_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC.bw"
	output:
		zNorm = "BigWigs/ZNormalized/{sample}_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC_zNorm.bw",
		zStats = "logs/{sample}.zNorm"
	params: pipePath = '.faire-pipeline', moduleVer = rVer
#	script:
#		"zNorm_from_Snakemake.R"
	shell:
		"""
		module purge && module load {params.moduleVer}
		Rscript --vanilla {params.pipePath}/zNorm.r {input} {output.zNorm} > {output.zStats}
		"""

if os.path.isdir == False:
	os.mkdir('Peakfiles')

rule CallPooledPeaks:
# Call Peaks using MACS. Pools all input files first.
# Because macs2 outputs multiple files and parsing it is annoying, 
# snakemake will just create a hidden file ('Peakfiles/.peakCall.done')
# at the end, which is requested by rule all
	input:
		expand("Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bed", sample = SAMPLE)
	params:
		outdir = peakDir, 
		control = CTRLPATH,
		name = expand("{dirID}_{nFiles}Reps_PooledPeaks", dirID = dirID, nFiles = nFiles),
		moduleVer = macsVer 
	output:
		touch("Peakfiles/.peakCall.done")
		#expand("{outdir}/{name}{fType}", outdir = peakDir, name = peakOutName, fType = ['_peaks.narrowPeak', '_peaks.xls', '_summits.bed'])
	shell:
		"""
		module purge && module load {params.moduleVer}
		macs2 callpeak -t {input}  -c {params.control} -n {params.name} -g dm --nomodel --extsize 125 --seed 123 --outdir {params.outdir} 
		"""

rule qcReport:
# use multiqc to generate report of pipeline run.
# Ignores output in any files appended as *.out or *.err, which by convention should be 
# any stdout or stderr files generated by the job scheduler (not the pipeline)
# The idea here is that any logs we want information on will be explicitly asked for by the 
# pipeline itself and stored in the 'logs/' directory or elsewhere as appropriate
	input:
		expand("Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.{ext}", sample = SAMPLE, ext = ['bam', 'bam.bai', 'bed']),
		expand("PCRdups/{sample}_PCR_duplicates", sample = SAMPLE)
	output:
		"multiqc_report.html"
	params: moduleVer = python3Ver , reportName = dirID
	shell:
		"""
		module purge && module load {params.moduleVer}
		multiqc . -f -x *.out -x *.err --filename {params.dirID}_report
		"""


#First sort peak file by q-value in decreasing order, then cut out chromosome, start and end coordinates, peak name and q-value and take the top 5000 peaks based on q-value
#sort -n -r -k9 ${STRAIN}_q5_sorted_dupsRemoved_noYUHet_peaks.narrowPeak | cut -f1,2,3,4,9 | head -5000 > ${STRAIN}_q5_sorted_dupsRemoved_noYUHet_Top5000Peaks.bed
	

