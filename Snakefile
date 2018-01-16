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

# Point directly to negative control genomic DNA input unless ctrlPath variable is set in config
if 'ctrlPath' in config:
	# set with --config ctrlPath="path/to/file"
	# Takes working directory, adds path to config as relative to entire filestructure.
	# I should actually change this back to being a relative path thing,
	# Because I fixed the error in the submission script
	rootPath = os.getcwd().split('/')[0:-1]
	rootPath.append(config["ctrlPath"])

	CTRLPATH = '/'.join(rootPath)

else:
	CTRLPATH = str('/proj/mckaylab/genomeFiles/' + GenomeAssembly + '/ControlGenomicDNA/ControlGenomicDNA_q5_sorted_dupsRemoved_noYUHet.bed') 
	

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



def getSamples(fastqList):
	sampleList = []
	for f in fastqList:
		sample = f.split('.')[0]
		sampleList.append(sample)
	return(sampleList)

FASTQ = glob.glob('*.fastq.gz')

SAMPLE = getSamples(FASTQ)

# Parse info for peak calling
dirID = os.getcwd().split('/')[-1] 	# get name of directory (should be sample Pool name)
nFiles = len(SAMPLE) 			# number of files in directory should be number of biological replicates (pool technical replicats first!)
peakDir = 'Peakfiles/' 			# be sure to include '/' 

if peakDir[-1] != '/':
	# die if peakDir set incorrectly
	print('ERROR: peakDir (' + peakDir + ') does not end with /')
	quit()

rule all:
# Only creates pooled files if there are files to pool.
	input: 
		expand("BigWigs/ZNormalized/{sample}_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC_zNorm.bw", sample = SAMPLE),
		expand("BigWigs/ZNormalized/{dirID}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC_zNorm.bw", dirID = dirID, nFiles = nFiles),
		expand("Peakfiles/.{nFiles}Reps_peakCall.done", nFiles = nFiles),
		expand("Peakfiles/.{nFiles}Reps_peakSort.done", nFiles = nFiles),
                expand("Peakfiles/.{sample}_peakCall.done", sample = SAMPLE),
		expand("{dirID}_report.html", dirID = dirID)
		if nFiles > 1 else
		expand("BigWigs/ZNormalized/{sample}_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC_zNorm.bw", sample = SAMPLE),
		expand("BigWigs/ZNormalized/{dirID}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC_zNorm.bw", dirID = dirID, nFiles = nFiles),
		expand("Peakfiles/.{nFiles}Reps_peakCall.done", nFiles = nFiles),
		expand("Peakfiles/.{nFiles}Reps_peakSort.done", nFiles = nFiles),
                expand("Peakfiles/.{sample}_peakCall.done", sample = SAMPLE),
		expand("{dirID}_report.html", dirID = dirID),
		expand("Bam/{dirID}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet.bed", dirID = dirID, nFiles = nFiles)

rule align:
	input:
		"{sample}.fastq.gz"
	output:
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
	shell:
		"""
		module purge && module load {params.moduleVer}
		Rscript --vanilla {params.pipePath}/zNorm.r {input} {output.zNorm} > {output.zStats}
		"""

if os.path.isdir == False:
	os.mkdir('Peakfiles')

rule CallPooledPeaks:
# Call Peaks using MACS. Pools all input files first.
# Then sort peak file by q-value in decreasing order, then output chromosome, start and end coordinates, peak name and q-value to bedfile (for subsetting into top x# peaks) 
# Because macs2 outputs multiple files and parsing it is annoying, 
# snakemake will just create a hidden file ('Peakfiles/.{nFiles}Reps_peakCall.done')
# at the end, which is requested by rule all. nFiles is included to allow recalling when more samples are added to directory upon rerunning pipeline.
	input:
		expand("Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bed", sample = SAMPLE)
	params:
		outdir = peakDir, 
		control = CTRLPATH,
		name = expand("{dirID}_{nFiles}Reps_FAIRE_PooledPeaks", dirID = dirID, nFiles = nFiles) if nFiles > 1 else expand("{sample}_FAIRE", sample = SAMPLE),
		moduleVer = macsVer 
	output:
		touch(expand("Peakfiles/.{nFiles}Reps_peakCall.done", nFiles = nFiles))
	shell:
		"""
		module purge && module load {params.moduleVer}
		macs2 callpeak -t {input}  -c {params.control} -n {params.name} -g dm --nomodel --extsize 125 --seed 123 --outdir {params.outdir}
		"""

rule CallRepPeaks:
    input:
            "Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bed"
    params:
            outdir = peakDir,
            control = CTRLPATH,
            name = "{sample}_ATAC_Peaks",
            moduleVer = macsVer
    output:
            touch("Peakfiles/.{sample}_peakCall.done")
    shell:
            """
            module purge && module load {params.moduleVer}
            macs2 callpeak -t {input} -c {params.control} -n {params.name} -g dm --nomodel --extsize 125 --seed 123 --outdir {params.outdir}
            """
rule sortPooledPeaks:
	input:
		expand("Peakfiles/.{nFiles}Reps_peakCall.done", nFiles = nFiles)
	output:
		touch(expand("Peakfiles/.{nFiles}Reps_peakSort.done", nFiles = nFiles)),
		name = expand("Peakfiles/{dirID}_{nFiles}Reps_FAIRE_PooledPeaks_peaks_qSorted.bed", dirID = dirID, nFiles = nFiles) if nFiles > 1 else expand("{sample}_FAIRE", sample = SAMPLE),
	params:
		name = expand("Peakfiles/{dirID}_{nFiles}Reps_FAIRE_PooledPeaks_peaks.narrowPeak", dirID = dirID, nFiles = nFiles)
	shell:
		"sort -n -r -k9 {params.name} | cut -f1,2,3,4,9 > {output.name}"
		

rule mergeBams:
	input:
		expand("Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.bam", sample = SAMPLE)
	params: moduleVer = samtoolsVer
	output:
		bam = expand("Bam/{dirID}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet.bam", dirID = dirID, nFiles = nFiles),
		idx = expand("Bam/{dirID}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet.bam.bai", dirID = dirID, nFiles = nFiles)
	shell:
		"""
		module purge && module load {params.moduleVer}
		samtools merge {output.bam} {input} &&
		samtools index {output.bam}
		"""
rule mergeBamtoBed:
# Convert the bam file into a bed file
	input:
		bam = expand("Bam/{dirID}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet.bam", dirID = dirID, nFiles = nFiles),
		idx = expand("Bam/{dirID}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet.bam.bai", dirID = dirID, nFiles = nFiles)
	output:
		expand("Bam/{dirID}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet.bed", dirID = dirID, nFiles = nFiles),
	params: moduleVer = bedtoolsVer
	shell:
		"""
		module purge && module load {params.moduleVer}
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
	params: moduleVer = deeptoolsVer, blacklist = BLACKLIST, genomeSize = GENOMESIZE
	shell:
		"""
		module purge && module load {params.moduleVer}
		bamCoverage -b {input.bam} -p {threads} --blackListFileName {params.blacklist} --normalizeTo1x {params.genomeSize} --outFileFormat bigwig --binSize 10 -e 125 -o {output}
		"""

rule mergeZNormBigWig:
# Z-Normalize Bigwig Files
	input:
		bw = expand("BigWigs/{dirID}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC.bw", dirID = dirID, nFiles = nFiles)
	output:
		zNorm = expand("BigWigs/ZNormalized/{dirID}_{nFiles}Reps_POOLED_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC_zNorm.bw", dirID = dirID, nFiles = nFiles),
		zStats = expand("logs/{dirID}_{nFiles}Reps_POOLED.zNorm", dirID = dirID, nFiles = nFiles)
	params: pipePath = '.faire-pipeline', moduleVer = rVer
	shell:
		"""
		module purge && module load {params.moduleVer}
		Rscript --vanilla {params.pipePath}/zNorm.r {input.bw} {output.zNorm} > {output.zStats}
		"""

rule qcReport:
# use multiqc to generate report of pipeline run.
# Ignores output in any files appended as *.out or *.err, which by convention should be 
# any stdout or stderr files generated by the job scheduler (not the pipeline)
# The idea here is that any logs we want information on will be explicitly asked for by the 
# pipeline itself and stored in the 'logs/' directory or elsewhere as appropriate
	input:
		expand("Bam/{sample}_q5_sorted_dupsRemoved_noYUHet.{ext}", sample = SAMPLE, ext = ['bam', 'bam.bai']),
		expand("PCRdups/{sample}_PCR_duplicates", sample = SAMPLE)
	output:
		expand("{dirID}_report.html", dirID = dirID)
	params: moduleVer = python3Ver , reportName = dirID
	shell:
		"""
		module purge && module load {params.moduleVer}
		multiqc . -f -x *.out -x *.err --filename {params.reportName}_report
		"""


	

