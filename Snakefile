#!/usr/bin/python3
import pandas as pd
import preProcessSampleConfig as pre

configfile: 'config.json'

file_info_path = config['sampleInfo']
basename_columns = config['baseNameColumns']
pool_basename_columns = config['poolBaseNameColumns']
is_paired_end = config['pairedEnd']
READS = ['1','2'] if is_paired_end else ['1']

REFGENOME = config['refGenome']
SPIKEGENOME = config['spikeGenome']
REFGENOMEPATH = config['genome'][REFGENOME]['bowtie']
BLACKLISTPATH = config['genome'][REFGENOME]['blacklistPath']
SPIKEGENOMEPATH = config['genome'][SPIKEGENOME]['bowtie']
controlDNAPath  = config['genome'][REFGENOME]['controlDNAPath']
chromSize_Path  = config['genome'][REFGENOME]['chrSize']

GENOMESIZE = config['genome'][REFGENOME]['genomeSize']
readLen = config['readLen']

modules = config['module']
#########
# Validation

#if is_paired_end:
#	sys.exit("paired end mode is not fully supported yet")

if os.path.exists(file_info_path) == False:
	print('Error: {name} does not exist. Be sure to set `sampleInfo` in config.json.'.format(name = file_info_path))

#########
# Generating sampleSheet outputs

speciesList  = [REFGENOME, SPIKEGENOME]
indexDict    = {REFGENOME: REFGENOMEPATH, SPIKEGENOME: SPIKEGENOMEPATH}

normTypeList = ['rpgcNorm']

sampleInfo, sampleSheet = pre.makeSampleSheets(file_info_path, basename_columns, "-", fileDelimiter = config['sampleInfoDelimiter'])
# add pool name column
sampleSheet = pre.addBaseName(sampleSheet, pool_basename_columns, "-", "poolBaseName")
# Create Pooled Sample Sheet
## leveraging the fact that sampleSheet is already cleaned
poolSampleSheet = sampleSheet[pool_basename_columns].copy()
poolSampleSheet = pre.addBaseName(poolSampleSheet, pool_basename_columns, "-").drop_duplicates()

sampleSheet['fastq_trim_r1'] = expand("Fastq/{sample}_R{num}_trim.fastq.gz", sample = sampleSheet.baseName, num = ['1'])
if is_paired_end: sampleSheet['fastq_trim_r2'] = expand("Fastq/{sample}_R{num}_trim.fastq.gz", sample = sampleSheet.baseName, num = ['2']) 
sampleSheet['bam']           = expand("Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved.{ftype}", sample = sampleSheet.baseName, species = REFGENOME, ftype = {"bam"})
sampleSheet['peaks']         = expand("Peaks/{sample}_{species}_trim_q5_sorted_dupsRemoved_peaks.narrowPeak", sample = sampleSheet.baseName, species = REFGENOME)
sampleSheet['bed']           = expand('Bed/{sample}_{species}_trim_q5_sorted_dupsRemoved.bed', sample = sampleSheet.baseName, species = REFGENOME)

poolSampleSheet['bam']           = expand("Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved_POOLED.{ftype}", sample = poolSampleSheet.baseName, species = REFGENOME, ftype = {"bam"})
poolSampleSheet['peaks']         = expand("Peaks/{sample}_{species}_trim_q5_sorted_dupsRemoved_POOLED_peaks.narrowPeak", sample = poolSampleSheet.baseName, species = REFGENOME)
poolSampleSheet['bed']           = expand('Bed/{sample}_{species}_trim_q5_sorted_dupsRemoved_POOLED.bed', sample = poolSampleSheet.baseName, species = REFGENOME)

for norm in normTypeList:

	# Add column per bigwig
	bw_colName = 'bigwig_{norm}'.format(norm = norm)
	sampleSheet[bw_colName] = expand("BigWig/{sample}_{species}_trim_q5_sorted_dupsRemoved_{norm}.bw", sample = sampleSheet.baseName, species = REFGENOME, norm = norm)
	poolSampleSheet[bw_colName] = expand("BigWig/{sample}_{species}_trim_q5_sorted_dupsRemoved_POOLED_{norm}.bw", sample = poolSampleSheet.baseName, species = REFGENOME, norm = norm)

	# Add column per zNorm bigwig
	bw_colName = 'bigwig_{norm}_zNorm'.format(norm = norm)
	sampleSheet[bw_colName] = expand("BigWig/{sample}_{species}_trim_q5_sorted_dupsRemoved_{norm}_zNorm.bw", sample = sampleSheet.baseName, species = REFGENOME, norm = norm)
	poolSampleSheet[bw_colName] = expand("BigWig/{sample}_{species}_trim_q5_sorted_dupsRemoved_POOLED_{norm}_zNorm.bw", sample = poolSampleSheet.baseName, species = REFGENOME, norm = norm)

sampleSheet.to_csv('sampleSheet.tsv', sep = "\t", index = False)
poolSampleSheet.to_csv('sampleSheetPooled.tsv', sep = "\t", index = False)

output_files = []
output_files.append(sampleSheet['fastq_trim_r1'])
if is_paired_end: output_files.append(sampleSheet['fastq_trim_r2'])
output_files.append(sampleSheet['bam'])
output_files.append(sampleSheet['peaks'])
output_files.append(sampleSheet['bed'])
output_files.append(sampleSheet['bigwig_rpgcNorm_zNorm'])
output_files.append(poolSampleSheet['bam'])
output_files.append(poolSampleSheet['peaks'])
output_files.append(poolSampleSheet['bed'])
output_files.append(poolSampleSheet['bigwig_rpgcNorm_zNorm'])

# Unpack nested lists
output_files = [item for sublist in output_files for item in sublist]
# Add additional files not from data.frame
output_files.append("multiqc_report.html")
#print(output_files)

rule all:
	input:
    		output_files,
#				expand("Fastq/{sample}_R{num}.fastq.gz", sample = sampleSheet.baseName, num = READS)

#constraing the num wildcard for identifying R1, R2 
# was running into ambigious wildcard issues between combine_tech_reps and trim ruls that I wasn't able to resolve -MN
wildcard_constraints:
	num = "[0-9]"

#function to get the /proj/HTSF paths for either r1 or r1+r2 
def getFastq(wildcards):
	if is_paired_end:
		r1 = sampleInfo[sampleInfo.baseName == wildcards.sample].fastq_r1.iloc[0]
		r2 = sampleInfo[sampleInfo.baseName == wildcards.sample].fastq_r2.iloc[0]
		return [r1, r2]
	else:
		r1 = sampleInfo[sampleInfo.baseName == wildcards.sample].fastq_r1.iloc[0]
		return [r1]

#passes getFastq function as input -- either 1 or 2 reads depending on is_paired_end
rule combine_technical_reps:
	input:
		getFastq
	output:
		expand("Fastq/{{sample}}_R{num}.fastq.gz", num = READS)
	run:
		if is_paired_end:
			shell("cat {input[0]} > {output[0]} && cat {input[1]} > {output[1]}")
		else:
			shell("cat {input[0]} > {output[0]}")

#using conditional run to trim either PE or SE reads
#paired-end reads contained in separate files (which they are in our normal case) must be trimmed together 
#pairs are either kept or tossed, and never separated, otherwise you end up with unequal read numbers between reads and alignment fails
rule trim_adapter:
	input:
		expand("Fastq/{{sample}}_R{num}.fastq.gz", num = READS)
	output:
		expand("Fastq/{{sample}}_R{num}_trim.fastq.gz", num = READS)
	log:
		adapterStats = 'Logs/{sample}_adapterStats',
		trimStats = 'Logs/{sample}_trimStats'
	params:
		bbmap = modules['bbmapVer']
	run:
		if is_paired_end:
			shell("module purge && module load {params.bbmap} && bbduk.sh in1={input[0]} in2={input[1]} out1={output[0]} out2={output[1]} ktrim=r ref=adapters rcomp=t tpe=t tbo=t hdist=1 mink=11 stats={log.adapterStats} > {log.trimStats}")
		else:
			shell("module purge && module load {params.bbmap} && bbduk.sh in={input} out={output} ktrim=r ref=adapters rcomp=t tpe=t tbo=t hdist=1 mink=11 stats={log.adapterStats} > {log.trimStats}")

#will only run if is_paired_end == false - as determined by input called for by sam2bam
rule align_se:
	input:
		"Fastq/{sample}_R1_trim.fastq.gz"
	output:
		sam = "Sam/{sample}_{species}_trim_se.sam",
		logInfo = "logs/{sample}_{species}_trim_se_bowtie2.txt"
	threads: 8
	params:
		refgenome = lambda wildcards: indexDict[wildcards.species]
	envmodules:
		modules['bowtie2Ver']
	shell:
		"""
		bowtie2 --seed 123 --very-sensitive -x {params.refgenome} -p {threads} -U {input} -S {output.sam} 2> {output.logInfo}
		"""

#will only run if is_paired_end == true - as determined by input called for by sam2bam
rule align_pe:
	input:
		r1 = "Fastq/{sample}_R1_trim.fastq.gz",
		r2 = "Fastq/{sample}_R2_trim.fastq.gz"
	output:
		sam = "Sam/{sample}_{species}_trim_pe.sam",
		logInfo = "logs/{sample}_{species}_trim_pe_bowtie2.txt"
	threads: 8
	params:
		refgenome = lambda wildcards: indexDict[wildcards.species]
	envmodules:
		modules['bowtie2Ver']
	shell:
		"""
		bowtie2 --seed 123 -x {params.refgenome} -p {threads} -1 {input.r1} -2 {input.r2} -S {output.sam} 2> {output.logInfo}
		"""


# TODO: combine w/ alignment step:
# bowtie2 --flags 2> {log} | samtools view -@ 4 -b - > {output.bam} && {flagstat_command}

rule sam2bam:
	input:
		"Sam/{sample}_{species}_trim_pe.sam" if is_paired_end else "Sam/{sample}_{species}_trim_se.sam"
	output:
		bam = "Bam/{sample}_{species}_trim.bam",
		flagstat = "logs/{sample}_{species}_trim.flagstat"
	params:
		blklist = BLACKLISTPATH
	envmodules:
		modules['samtoolsVer'],
		modules['bedtoolsVer']
	shell:
		"""
		samtools view -@ 4 -b {input} | bedtools intersect -a stdin -b {params.blklist} -v > {output.bam} &&
		samtools flagstat {output.bam} > {output.flagstat}
		"""
rule qFilter:
	input:
		"Bam/{sample}_{species}_trim.bam"
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
# TODO: remove indexing step if no longer needed?
# TODO: samtools view -@ 4 ? don't think multithreading will speed it up...
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
		modules['picardVer'],
		modules['samtoolsVer']
	shell:
		"""
		java -Xmx8g -jar {params.picardPath} MarkDuplicates INPUT= {input} OUTPUT= {output.markedDups} METRICS_FILE= {output.PCRdups} REMOVE_DUPLICATES= false ASSUME_SORTED= true &&
		samtools view -bF 0x400 {output.markedDups} > {output.bam} &&
		samtools flagstat {output.bam} > {output.flagstat} &&
		samtools index {output.bam}
		"""

rule bamToBed:
# Convert the bam file into a bed file
	input:
		bam = "Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved.bam",
		idx = "Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved.bam.bai"
	output:
		"Bed/{sample}_{species}_trim_q5_sorted_dupsRemoved.bed"
	envmodules:
		modules['bedtoolsVer']
	shell:
		"""
		bedtools bamtobed -i {input.bam} > {output}
		"""

rule bigWig:
# Bam Coverage to output bigwig file normalized to genomic coverage
	input:
		bam = "Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved.bam",
		idx = "Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved.bam.bai"
	output:
		"BigWig/{sample}_{species}_trim_q5_sorted_dupsRemoved_rpgcNorm.bw"
	threads: 8 
	params: genomeSize = GENOMESIZE
	envmodules:
		modules['deeptoolsVer']
	shell:
		"""
		bamCoverage -b {input.bam} -p {threads} --normalizeTo1x {params.genomeSize} --outFileFormat bigwig --binSize 10 -e 125 -o {output}
		"""

rule zNormBigWig:
# Z-Normalize Bigwig Files
	input:
		bw = "BigWig/{sample}_{species}_trim_q5_sorted_dupsRemoved_rpgcNorm.bw"
	output:
		zNorm = "BigWig/{sample}_{species}_trim_q5_sorted_dupsRemoved_rpgcNorm_zNorm.bw",
		zStats = "logs/{sample}_{species}.zNorm"
	envmodules:
		modules['rVer']
	shell:
		"""
		Rscript --vanilla scripts/zNorm.r {input} {output.zNorm} > {output.zStats}
		"""

rule mergeBams:
	input:
		lambda wildcards: ["Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved.bam".format(sample = s, species = REFGENOME) for s in sampleSheet[sampleSheet.poolBaseName == wildcards.samplePool].baseName]
	output:
		bam = "Bam/{samplePool}_{species}_trim_q5_sorted_dupsRemoved_POOLED.bam",
		idx = "Bam/{samplePool}_{species}_trim_q5_sorted_dupsRemoved_POOLED.bam.bai"
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
		bam = "Bam/{samplePool}_{species}_trim_q5_sorted_dupsRemoved_POOLED.bam",
		idx = "Bam/{samplePool}_{species}_trim_q5_sorted_dupsRemoved_POOLED.bam.bai"
	output:
		"Bed/{samplePool}_{species}_trim_q5_sorted_dupsRemoved_POOLED.bed"
	envmodules:
		modules['bedtoolsVer']
	shell:
		"""
		bedtools bamtobed -i {input.bam} > {output}
		"""

rule CallRepPeaks:
    input:
    	"Bed/{sample}_{species}_trim_q5_sorted_dupsRemoved.bed"
    output:
    	"Peaks/{sample}_{species}_trim_q5_sorted_dupsRemoved_peaks.narrowPeak"
    params:
    	prefix = "Peaks/{sample}_{species}_trim_q5_sorted_dupsRemoved",
	control = controlDNAPath
    envmodules:
        modules['macsVer']
    shell:
    	"""
			macs2 callpeak -t {input} -c {params.control} -n {params.prefix} -g dm --nomodel --extsize 125 --seed 123
			"""

rule CallPooledPeaks:
# Call Peaks using MACS. Pools all input files first.
# Then sort peak file by q-value in decreasing order, then output chromosome, start and end coordinates, peak name and q-value to bedfile (for subsetting into top x# peaks)
# Because macs2 outputs multiple files and parsing it is annoying,
# snakemake will just create a hidden file ('Peakfiles/.{nFiles}Reps_peakCall.done')
# at the end, which is requested by rule all. nFiles is included to allow recalling when more samples are added to directory upon rerunning pipeline.
	input:
		"Bed/{samplePool}_{species}_trim_q5_sorted_dupsRemoved_POOLED.bed"
	output:
		"Peaks/{samplePool}_{species}_trim_q5_sorted_dupsRemoved_POOLED_peaks.narrowPeak"
	params:
		control = controlDNAPath,
		name = "Peaks/{samplePool}_{species}_trim_q5_sorted_dupsRemoved_POOLED"
	envmodules:
		modules['macsVer']
	shell:
		"""
		macs2 callpeak -t {input}  -c {params.control} -n {params.name} -g dm --nomodel --extsize 125 --seed 123
		"""


rule sortPooledPeaks:
	input:
		"Peaks/{samplePool}_{species}_trim_q5_dupsRemoved_POOLED_peaks.narrowPeak"
	output:
		"Peaks/{samplePool}_{species}_trim_q5_dupsRemoved_POOLED_peaks_qSorted.narrowPeak"
	shell:
		"sort -n -r -k9 {input} | cut -f1,2,3,4,9 > {output}"
		

rule mergeBigWig:
# Bam Coverage to output bigwig file normalized to genomic coverage
	input:
		bam = "Bam/{samplePool}_{species}_trim_q5_sorted_dupsRemoved_POOLED.bam",
		idx = "Bam/{samplePool}_{species}_trim_q5_sorted_dupsRemoved_POOLED.bam.bai"
	output:
		"BigWig/{samplePool}_{species}_trim_q5_sorted_dupsRemoved_POOLED_rpgcNorm.bw"
	threads: 8 
	params: genomeSize = GENOMESIZE
	envmodules:
		modules['deeptoolsVer']
	shell:
		"""
		bamCoverage -b {input.bam} -p {threads} --normalizeTo1x {params.genomeSize} --outFileFormat bigwig --binSize 10 -e 125 -o {output}
		"""

rule mergeZNormBigWig:
# Z-Normalize Bigwig Files
	input:
		bw = "BigWig/{samplePool}_{species}_trim_q5_sorted_dupsRemoved_POOLED_rpgcNorm.bw"
	output:
		zNorm = "BigWig/{samplePool}_{species}_trim_q5_sorted_dupsRemoved_POOLED_rpgcNorm_zNorm.bw",
		zStats = "logs/{samplePool}_{species}_POOLED_rpgcNorm.zNorm"
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
		expand("Bam/{sample}_{species}_trim_q5_sorted_dupsRemoved.{ext}", sample = sampleSheet.baseName, ext = ['bam', 'bam.bai'], species = REFGENOME),
		expand("PCRdups/{sample}_{species}_PCR_duplicates", sample = sampleSheet.baseName, species = REFGENOME)
	output:
		"multiqc_report.html"
	envmodules:
		modules['multiqcVer']
	shell:
		"""
		multiqc . -f -x *.out -x *.err
		"""


	

