#!/usr/bin/python3
import pandas as pd
import preProcessSampleConfig as pre

configfile: 'config.json'

file_info_path = config['sampleInfo']
basename_columns = config['baseNameColumns']
pool_basename_columns = config['poolBaseNameColumns']
is_paired_end = config['pairedEnd']

REFGENOME = config['refGenome']
SPIKEGENOME = config['spikeGenome']
REFGENOMEPATH = config['genome'][REFGENOME]['bowtie']
SPIKEGENOMEPATH = config['genome'][SPIKEGENOME]['bowtie']
controlDNAPath  = config['genome'][REFGENOME]['controlDNAPath']
chromSize_Path  = config['genome'][REFGENOME]['chrSize']

GENOMESIZE = config['genome'][REFGENOME]['genomeSize']
readLen = config['readLen']

modules = config['module']
#########
# Validation

if is_paired_end:
	sys.exit("paired end mode is not fully supported yet")

if os.path.exists(file_info_path) == False:
	print('Error: {name} does not exist. Be sure to set `sampleInfo` in config.json.'.format(name = file_info_path))

#########
# Generating sampleSheet outputs

speciesList  = [REFGENOME, SPIKEGENOME]
indexDict    = {REFGENOME: REFGENOMEPATH, SPIKEGENOME: SPIKEGENOMEPATH}

normTypeList = ['unNorm', 'rpgcNorm']

sampleInfo, sampleSheet = pre.makeSampleSheets(file_info_path, basename_columns, "-", fileDelimiter = config['sampleInfoDelimiter'])
# Create Pooled Sample Sheet
## leveraging the fact that sampleSheet is already cleaned
poolSampleSheet = sampleSheet[pool_basename_columns].copy()
poolSampleSheet = pre.addBaseName(poolSampleSheet, pool_basename_columns, config['sampleInfoDelimiter']).drop_duplicates()

sampleSheet['fastq_trim_r1'] = expand("Fastq/{sample}_R{num}_trim.fastq.gz", sample = sampleSheet.baseName, num = ['1'])
# TODO: add support for paired end
#sampleSheet['fastq_trim_r2'] = expand("Fastq/{sample}_R{num}_trim.fastq.gz", sample = sampleSheet.baseName, num = ['2'])
sampleSheet['bam']           = expand("Bam/{sample}_{species}_trim_q5_dupsRemoved.{ftype}", sample = sampleSheet.baseName, species = REFGENOME, ftype = {"bam"})
sampleSheet['peaks']         = expand("Peaks/{sample}_{species}_trim_q5_dupsRemoved_peaks.narrowPeak", sample = sampleSheet.baseName, species = REFGENOME)
sampleSheet['bed']           = expand('Bed/{sample}_{species}_trim_q5_dupsRemoved.bed', sample = sampleSheet.baseName, species = REFGENOME)

poolSampleSheet['bam']           = expand("Bam/{sample}_{species}_trim_q5_dupsRemoved_POOL.{ftype}", sample = poolSampleSheet.baseName, species = REFGENOME, ftype = {"bam"})
poolSampleSheet['peaks']         = expand("Peaks/{sample}_{species}_trim_q5_dupsRemoved_peaks_POOL.narrowPeak", sample = poolSampleSheet.baseName, species = REFGENOME)
poolSampleSheet['bed']           = expand('Bed/{sample}_{species}_trim_q5_dupsRemoved_POOL.bed', sample = poolSampleSheet.baseName, species = REFGENOME)

for norm in normTypeList:

	# Add column per bigwig
	bw_colName = 'bigwig_{norm}'.format(norm = norm)
	sampleSheet[bw_colName] = expand("BigWig/{sample}_{species}_trim_q5_dupsRemoved_{norm}.bw", sample = sampleSheet.baseName, species = REFGENOME, norm = norm)
	poolSampleSheet[bw_colName] = expand("BigWig/{sample}_{species}_trim_q5_dupsRemoved_POOL_{norm}.bw", sample = poolSampleSheet.baseName, species = REFGENOME, norm = norm)

	# Add column per zNorm bigwig
	bw_colName = 'bigwig_{norm}_zNorm'.format(norm = norm)
	sampleSheet[bw_colName] = expand("BigWig/{sample}_{species}_trim_q5_dupsRemoved_{norm}_zNorm.bw", sample = sampleSheet.baseName, species = REFGENOME, norm = norm)
	poolSampleSheet[bw_colName] = expand("BigWig/{sample}_{species}_trim_q5_dupsRemoved_POOL_{norm}_zNorm.bw", sample = poolSampleSheet.baseName, species = REFGENOME, norm = norm)

sampleSheet.to_csv('sampleSheet.tsv', sep = "\t", index = False)
poolSampleSheet.to_csv('sampleSheetPooled.tsv', sep = "\t", index = False)

#################################################
# OLD SECTION
##################################################
#import os
#import glob
#
#GenomeAssembly = 'dm3'
#
###############################
## Module Versions:
#
#bowtie2Ver = str('bowtie2/2.2.8')
#samtoolsVer = str('samtools/1.3.1')
#bedtoolsVer = str('bedtools/2.25.0')
#
#picardVer = str('2.2.4')
##picardPath = str('/nas02/apps/picard-' + picardVer + '/picard-tools-' + picardVer + '/picard.jar')
#picardPath = str('/nas/longleaf/apps/picard/' + picardVer + '/picard-tools-' + picardVer + '/picard.jar')
#
#deeptoolsVer = str('deeptools/2.4.1')
#macsVer = str('macs/2016-02-15')
#ucscVer = str('ucsctools/320')
#rVer = str('r/3.3.1')
#
#python3Ver = str('python/3.5.1')
#
###############################
#
## Dictionary of effective genome sizes.
## Used as lookup table later for RPGC adjustments
## based on GenomeAssembly
#effectiveGenomeSizes = {'dm3':121400000}
#
## Symlink to pipeline path for calling scripts with commandline arguments (zNorm.r for example)
## Requires that snakemake call be structured as `snakemake --snakefile <path/to/snakefile>`
#PIPEPATH = os.path.dirname(os.path.abspath(sys.argv[2]))
#symLink = str(os.getcwd() + '/.faire-pipeline')
#
#if os.path.isdir(symLink) == False:
#	os.symlink(PIPEPATH, symLink)
#
#if os.path.exists(symLink + '/Snakefile') == False:
#	print('ERROR: no Snakefile in .faire-pipeline. Ensure --snakefile <path/to/snakefile> is 3rd argument of snakemake call.')
#	quit()
#
#
#
#REFGENEPATH=str('/proj/mckaylab/genomeFiles/' + GenomeAssembly + '/RefGenome/' + GenomeAssembly)	# Point directly to the refgeneome file you want to use
#
## Point directly to negative control genomic DNA input unless ctrlPath variable is set in config
#if 'ctrlPath' in config:
#	# set with --config ctrlPath="path/to/file"
#	# Takes working directory, adds path to config as relative to entire filestructure.
#	# I should actually change this back to being a relative path thing,
#	# Because I fixed the error in the submission script
#	rootPath = os.getcwd().split('/')[0:-1]
#	rootPath.append(config["ctrlPath"])
#
#	CTRLPATH = '/'.join(rootPath)
#
#else:
#	CTRLPATH = str('/proj/mckaylab/genomeFiles/' + GenomeAssembly + '/ControlGenomicDNA/ControlGenomicDNA_q5_sorted_dupsRemoved_noYUHet.bed')
#
#
#BLACKLIST=str('/proj/mckaylab/genomeFiles/'+ GenomeAssembly + '/' + GenomeAssembly + 'Blacklist.bed')
#
## effective genome size for RPGC normalization, might need to change depending on assembly/what you blacklist
#GENOMESIZE = effectiveGenomeSizes[GenomeAssembly]
#
#if os.path.exists(glob.glob(REFGENEPATH + '*')[0]) == False:
#	print('ERROR: REFGENEPATH (' + REFGENEPATH + ') does not contain bowtie genome files.')
#	quit()
#
#if os.path.exists(BLACKLIST) == False:
#	print('ERROR: BLACKLIST (' + BLACKLIST + ') does not exist.')
#	quit()
#
#if os.path.exists(CTRLPATH) == False:
#	print('ERROR: ' + CTRLPATH + ' does not exist. Is CTRLPATH set correctly?')
#	quit()
#
#
#
#def getSamples(fastqList):
#	sampleList = []
#	for f in fastqList:
#		sample = f.split('.')[0]
#		sampleList.append(sample)
#	return(sampleList)
#
#FASTQ = glob.glob('*.fastq.gz')
#
#SAMPLE = getSamples(FASTQ)
#
## Parse info for peak calling
#dirID = os.getcwd().split('/')[-1] 	# get name of directory (should be sample Pool name)
#nFiles = len(SAMPLE) 			# number of files in directory should be number of biological replicates (pool technical replicats first!)
#peakDir = 'Peakfiles/' 			# be sure to include '/'
#
#if peakDir[-1] != '/':
#	# die if peakDir set incorrectly
#	print('ERROR: peakDir (' + peakDir + ') does not end with /')
#	quit()

#################################################
# END OLD SECTION
#################################################

output_files = []
output_files.append(sampleSheet['fastq_trim_r1'])
output_files.append(sampleSheet['bam'])
output_files.append(sampleSheet['peaks'])
output_files.append(sampleSheet['bed'])
output_files.append(sampleSheet['bigwig_unNorm_zNorm'])
output_files.append(sampleSheet['bigwig_rpgcNorm_zNorm'])
output_files.append(poolSampleSheet['bam'])
output_files.append(poolSampleSheet['peaks'])
output_files.append(poolSampleSheet['bed'])
output_files.append(poolSampleSheet['bigwig_unNorm_zNorm'])
output_files.append(poolSampleSheet['bigwig_rpgcNorm_zNorm'])
output_files.append("multiqc_report.html")


rule all:
    input: output_files



#rule all:
## Only creates pooled files if there are files to pool.
#	input:
#		expand("BigWigs/ZNormalized/{sample}_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC_zNorm.bw", sample = SAMPLE),
#		expand("BigWigs/ZNormalized/{dirID}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC_zNorm.bw", dirID = dirID, nFiles = nFiles),
#		expand("Peakfiles/.{nFiles}Reps_peakCall.done", nFiles = nFiles),
#		expand("Peakfiles/.{nFiles}Reps_peakSort.done", nFiles = nFiles),
#                expand("Peakfiles/.{sample}_peakCall.done", sample = SAMPLE),
#		expand("{dirID}_report.html", dirID = dirID)
#		if nFiles > 1 else
#		expand("BigWigs/ZNormalized/{sample}_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC_zNorm.bw", sample = SAMPLE),
#		expand("BigWigs/ZNormalized/{dirID}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC_zNorm.bw", dirID = dirID, nFiles = nFiles),
#		expand("Peakfiles/.{nFiles}Reps_peakCall.done", nFiles = nFiles),
#		expand("Peakfiles/.{nFiles}Reps_peakSort.done", nFiles = nFiles),
#                expand("Peakfiles/.{sample}_peakCall.done", sample = SAMPLE),
#		expand("{dirID}_report.html", dirID = dirID),
#		expand("Bam/{dirID}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet.bed", dirID = dirID, nFiles = nFiles)

# TODO: get this working for paired-end mode
rule combine_technical_reps:
	input:
		#r1 = lambda wildcards : sampleInfo[sampleInfo.baseName == wildcards.sample].fastq_r1,
		#r2 = lambda wildcards : sampleInfo[sampleInfo.baseName == wildcards.sample].fastq_r2
		#if is_paired_end else:
		r1 = lambda wildcards : sampleInfo[sampleInfo.baseName == wildcards.sample].fastq_r1,
	output:
		#r1 = 'Fastq/{sample}_R1.fastq.gz',
		#r2 = 'Fastq/{sample}_R2.fastq.gz'
		#if is_paired_end else:
		r1 = 'Fastq/{sample}_R1.fastq.gz'
	run:
		if is_paired_end:
			shell("cat {input.r1} > {output.r1} && cat {input.r2} > {output.r2}")
		else:
			shell("cat {input.r1} > {output.r1}")

rule trim_adapter:
	input:
		r1 = "Fastq/{sample}_R1.fastq.gz"
	output:
		r1 = "Fastq/{sample}_R1_trim.fastq.gz"
	log:
		adapterStats = 'Logs/{sample}_adapterStats',
		trimStats = 'Logs/{sample}_trimStats'
	envmodules:
		modules['bbmapVer']
	shell:
		"""
		bbduk.sh in={input.r1} out={output.r1} ktrim=r ref=adapters rcomp=t tpe=t tbo=t hdist=1 mink=11 stats={log.adapterStats} > {log.trimStats}
		"""

rule align:
	input:
		"{sample}_R1_trim.fastq.gz"
	output:
		sam = "Sam/{sample}_{species}_{species}.sam",
		logInfo = "logs/{sample}_{species}_bowtie2.txt"
	threads: 8
	params:
		refgenome = lambda wildcards: indexDict[wildcards.species]
	envmodules:
		modules['bowtie2Ver']
	shell:
		"""
		bowtie2 --seed 123 -x {params.refgenome} -p {threads} -U {input} -S {output.sam} 2> {output.logInfo}
		"""
# TODO: combine w/ alignment step:
# bowtie2 --flags 2> {log} | samtools view -@ 4 -b - > {output.bam} && {flagstat_command}
rule sam2bam:
	input:
		"Sam/{sample}_{species}.sam"
	output:
		bam = "Bam/{sample}_{species}.bam",
		flagstat = "logs/{sample}_{species}.flagstat"
	envmodules:
		modules['samtoolsVer']
	shell:
		"""
		samtools view -@ 4 -b {input} > {output.bam} &&
		samtools flagstat {output.bam} > {output.flagstat}
		"""
rule qFilter:
	input:
		"Bam/{sample}_{species}.bam"
	output:
		bam = "Bam/{sample}_{species}_trim_q5.bam",
		flagstat = "logs/{sample}_{species}_trim_q5.flagstat"
	envmodules:
		modules['samtoolsVer']
	shell:
		"""
		samtools view -@ 4 -bq 5 {input} > {output.bam} &&
		samtools flagstat {output.bam} > {output.flagstat}
		"""
rule bamSort:
	input:
		"Bam/{sample}_{species}_trim_q5.bam"
	output:
		"Bam/{sample}_{species}_trim_q5_sorted.bam"
	envmodules:
		modules['samtoolsVer']
	shell:
		"""
		samtools sort -@ 4 -o {output} {input}
		"""
rule removeDups:
# Remove the reads that are marked as pcr duplicates from the bam file using the bit flag for pcr dups
# Create index for bam file. This is needed for the next step to remove Chrm Y, U and Het
	input:
		"Bam/{sample}_{species}_trim_q5_sorted.bam"
	output:
		markedDups = "Bam/{sample}_{species}_trim_q5_sorted_dupsMarked.bam",
		PCRdups = "PCRdups/{sample}_{species}_PCR_duplicates",
		bam = "Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved.bam",
		flagstat = "logs/{sample}_{species}_trim_q5_sorted_dupsRemoved.flagstat",
		idx = "Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved.bam.bai"
	params: picardPath = modules['picardPath']
	envmodules:
		modules['picardVer']
	shell:
		"""
		java -Xmx8g -jar {params.picardPath} MarkDuplicates INPUT= {input} OUTPUT= {output.markedDups} METRICS_FILE= {output.PCRdups} REMOVE_DUPLICATES= false ASSUME_SORTED= true &&
		samtools view -@ 4 -bF 0x400 {input} > {output.bam} &&
		samtools flagstat {output.bam} > {output.flagstat} &&
		samtools index {output.bam}
		"""

rule noYUHet:
# Extract only reads from Chrm 2R,2L,3R,3L,4 and X
# Create new index for filtered bam file. Needed to convert bam file to bed file
	input:
		bam = "Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved.bam",
		idx = "Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved.bam.bai"
	output:
		bam = "Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved_noYUHet.bam",
		flagstat =  "logs/{sample}_{species}_trim_q5_sorted_dupsRemoved_noYUHet.flagstat",
		idx = "Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved_noYUHet.bam.bai"
	envmodules:
		modules['samtoolsVer']
	shell:
		"""
		samtools view -@ 4 -b {input.bam} chr2L chr2R chr3L chr3R chr4 chrX > {output.bam} &&
		samtools flagstat {output.bam} > {output.flagstat} &&
		samtools index {output.bam}
		"""

rule noYUHet_toBed:
# Convert the bam file into a bed file
	input:
		bam = "Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved_noYUHet.bam",
		idx = "Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved_noYUHet.bam.bai"
	output:
		"Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved_noYUHet.bed"
	envmodules:
		modules['bedtoolsVer']
	shell:
		"""
		bedtools bamtobed -i {input.bam} > {output}
		"""

rule bigWig:
# Bam Coverage to output bigwig file normalized to genomic coverage
	input:
		bam = "Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved_noYUHet.bam",
		idx = "Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved_noYUHet.bam.bai"
	output:
		"BigWigs/{sample}_{species}_trim_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC.bw"
	threads: 8 
	params: genomeSize = GENOMESIZE
	envmodules:
		modules['deeptoolsVer']
	shell:
		"""
		bamCoverage -b {input.bam} -p {threads} --normalizeTo1x {params.genomeSize} --outFileFormat bigwig --binSize 10 -e 125 -o {output}
		"""
# TODO: fix names
rule zNormBigWig:
# Z-Normalize Bigwig Files
	input:
		bw = "BigWigs/{sample}_{species}_trim_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC.bw"
	output:
		zNorm = "BigWigs/ZNormalized/{sample}_{species}_trim_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC_zNorm.bw",
		zStats = "logs/{sample}_{species}.zNorm"
	envmodules:
		modules['rVer']
	shell:
		"""
		Rscript --vanilla scripts/zNorm.r {input} {output.zNorm} > {output.zStats}
		"""

#TODO: fix names & nFiles crap here0
rule CallPooledPeaks:
# Call Peaks using MACS. Pools all input files first.
# Then sort peak file by q-value in decreasing order, then output chromosome, start and end coordinates, peak name and q-value to bedfile (for subsetting into top x# peaks) 
# Because macs2 outputs multiple files and parsing it is annoying, 
# snakemake will just create a hidden file ('Peakfiles/.{nFiles}Reps_peakCall.done')
# at the end, which is requested by rule all. nFiles is included to allow recalling when more samples are added to directory upon rerunning pipeline.
	input:
		expand("Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved_noYUHet.bed", sample = poolSampleSheet.baseName)
	output:
		"Peaks/{sample}_{species}_trim_q5_dupsRemoved_peaks_POOL.narrowPeak"
	params:
		control = CTRLPATH,
		name = "Peaks/{sample}_{species}_trim_q5_dupsRemoved_peaks_POOL"
	envmodules:
		modules['macsVer']
	shell:
		"""
		macs2 callpeak -t {input}  -c {params.control} -n {params.name} -g dm --nomodel --extsize 125 --seed 123
		"""

rule CallRepPeaks:
    input:
	"Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved_noYUHet.bed"
    output:
	"Peaks/{sample}_{species}_trim_q5_dupsRemoved_peaks.narrowPeak"
    params:
    	prefix = "Peaks/{sample}_{species}_trim_q5_dupsRemoved_peaks",
	control = controlDNAPath
    envmodules:
        modules['bedtoolsVer']
    shell:
    	"""
	macs2 callpeak -t {input} -c {params.control} -n {params.prefix} -g dm --nomodel --extsize 125 --seed 123
	"""

rule sortPooledPeaks:
	input:
		"Peaks/{sample}_{species}_trim_q5_dupsRemoved_peaks_POOL.narrowPeak"
	output:
		"Peaks/{sample}_{species}_trim_q5_dupsRemoved_peaks_POOL_qSorted.narrowPeak"
	shell:
		"sort -n -r -k9 {input} | cut -f1,2,3,4,9 > {output}"
		

rule mergeBams:
	input:
		expand("Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved_noYUHet.bam", sample = sampleSheet.baseName)
	output:
		bam = expand("Bam/{samplePool}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet.bam", samplePool = poolSampleSheet.baseName, nFiles = nFiles),
		idx = expand("Bam/{samplePool}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet.bam.bai", samplePool = poolSampleSheet.baseName, nFiles = nFiles)
	envmodules:
		modules['samtoolsVer']
	shell:
		"""
		samtools merge {output.bam} {input} &&
		samtools index {output.bam}
		"""
rule mergeBamtoBed:
# Convert the bam file into a bed file
	input:
		bam = expand("Bam/{samplePool}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet.bam", samplePool = poolSampleSheet.baseName, nFiles = nFiles),
		idx = expand("Bam/{samplePool}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet.bam.bai", samplePool = poolSampleSheet.baseName, nFiles = nFiles)
	output:
		"Bam/{samplePool}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet.bed"
	envmodules:
		modules['bedtoolsVer']
	shell:
		"""
		bedtools bamtobed -i {input.bam} > {output}
		"""

rule mergeBigWig:
# Bam Coverage to output bigwig file normalized to genomic coverage
	input:
		bam = expand("Bam/{dirID}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet.bam", dirID = dirID, nFiles = nFiles),
		idx = expand("Bam/{dirID}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet.bam.bai", dirID = dirID, nFiles = nFiles)
	output:
		expand("BigWigs/{dirID}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC.bw", dirID = dirID, nFiles = nFiles)
	threads: 8 
	params: genomeSize = GENOMESIZE
	envmodules:
		modules['deeptoolsVer']
	shell:
		"""
		bamCoverage -b {input.bam} -p {threads} --blackListFileName {params.blacklist} --normalizeTo1x {params.genomeSize} --outFileFormat bigwig --binSize 10 -e 125 -o {output}
		"""

rule mergeZNormBigWig:
# Z-Normalize Bigwig Files
	input:
		bw = expand("BigWigs/{dirID}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC.bw", dirID = dirID, nFiles = nFiles)
	output:
		zNorm = expand("BigWigs/ZNormalized/{dirID}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC_zNorm.bw", dirID = dirID, nFiles = nFiles),
		zStats = expand("logs/{dirID}_{nFiles}Reps_POOLED.zNorm", dirID = dirID, nFiles = nFiles)
	envmodules:
		modules['rVer']
	shell:
		"""
		Rscript --vanilla scripts/zNorm.r {input.bw} {output.zNorm} > {output.zStats}
		"""

rule qcReport:
# use multiqc to generate report of pipeline run.
# Ignores output in any files appended as *.out or *.err, which by convention should be 
# any stdout or stderr files generated by the job scheduler (not the pipeline)
# The idea here is that any logs we want information on will be explicitly asked for by the 
# pipeline itself and stored in the 'logs/' directory or elsewhere as appropriate
	input:
		expand("Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved_noYUHet.{ext}", sample = SAMPLE, ext = ['bam', 'bam.bai']),
		expand("PCRdups/{sample}_{species}_PCR_duplicates", sample = SAMPLE)
	output:
		expand("multiqc_report.html", dirID = dirID)
	envmodules:
		modules['pythonVer']
	shell:
		"""
		module purge && module load {params.moduleVer}
		multiqc . -f -x *.out -x *.err
		"""


	

