#!/bin/bash

##############################
# Set parameters:

USR=$USER						# Don't change this unless you have a good reason to

#NETSCR=/netscr/$USR/					# Location of input fastq. Be sure to end path with '/' otherwise pipeline will fail.

NETSCR=$(pwd)/						# Uncomment to use working directory as input & output dir
							
REFGENEPATH=~/RefGenome/dm3				# Point directly to the refgeneome file you want to use
CTRLPATH=$(pwd)/src/faire-pipeline/ControlGenomicDNA/ControlGenomicDNA_q5_sorted_dupsRemoved_noYUHet.bed # Point directly to negative control genomic DNA input
PIPEPATH=$(pwd)/src/faire-pipeline/

QUEUE=day						# BSUB Queue
stdOUT=$NETSCR/OutputFiles/				# standard output directory, end path with '/'
stdERR=$NETSCR/ErrorFiles/				# standard error directory, end path with '/'

##############################
# Module Versions:

bowtie2Ver=/2.2.8
samtoolsVer=/1.3.1
bedtoolsVer=/2.25.0
picardVer=/2.2.4
deeptoolsVer=/2.2.4
macsVer=/2015-06-03
ucscVer=/320

##############################

# Parses complete path to picard from version number
picardNum=$(echo $picardVer | cut -f2 -d/ )
picardPath=/nas02/apps/picard-$picardNum/picard-tools-$picardNum/picard.jar

# Parse commandline flags
STRAIN=${1%%.*}
ALIGN=$2
PEAK=$3

# Purge all currently loaded modules, load required modules for pipeline:
module bash purge
module bash load bowtie2$bowtie2Ver samtools$samtoolsVer bedtools$bedtoolsVer picard$picardVer deeptools$deeptoolsVer macs$macsVer ucsctools$ucscVer

# This series of if statements checks that the control .bed file defined by $CTRLPATH exists 
# and creates standard out/error and Stats directories if they don't already exist
# Also checks whether stdOUT, stdERR, and NETSCR end with '/', exits with errorcode 1 if so

if [[ ! -f $CTRLPATH ]]; then
	echo "Error: "$CTRLPATH" does not exist. Is \$CTRLPATH set correctly?"
	exit 1
fi

#stdOUT_tester=$(echo $stdOUT | rev)
#stdERR_tester=$(echo $stdERR | rev)
#NETSCR_tester=$(echo $NETSCR | rev)

#if [[ ${stdOUT_tester:0:1} != "/" ]]; then
if [[ ${stdOUT: -1} != "/" ]]; then
	echo "Error: "$stdOUT" must end in /"
	exit 1
fi

#if [[ ${stdERR_tester:0:1} != "/" ]]; then
if [[ ${stdERR: -1} != "/" ]]; then 
	echo "Error: "$stdERR" must end in /"
	exit 1
fi

#if [[ ${NETSCR_tester:0:1} != "/" ]]; then
if [[ ${NETSCR: -1} != "/" ]]
	echo "Error: "$NETSCR" must end in /"
	exit 1
fi

if [[ ! -d $stdOUT ]]; then
	mkdir $stdOUT
fi

if [[ ! -d $stdERR ]]; then
	mkdir $stdERR
fi

if [[ ! -d ./Stats ]]; then
	mkdir ./Stats
fi

if [[ ! -d ./Bam ]]; then
	mkdir ./Bam/
	mkdir ./Sam/
	mkdir ./Peakfiles/
	mkdir ./Bigwigs/
	mkdir ./PCRdups/
    mkdir ./BigWigs/ZNormalized/
fi

echo "

#!/bin/bash

#BSUB -q ${QUEUE}
#BSUB -o ${stdOUT}${STRAIN}_outfile.%J
#BSUB -e ${stdERR}${STRAIN}_errorfile.%J

">processFAIRESeqReadsAndCallMACSPeaks_${STRAIN}.bsub

if [ "${ALIGN}" = "Align" ]; then 

###echo "

#BSUB -n 8
#BSUB -R \"span[hosts=1]\"

# Change directory
# Create a directory in your scratch directory to store all data files 
# Then move fastq file to that directory

#mkdir ${NETSCR}${STRAIN}
#mv ${NETSCR}${STRAIN}.fastq.gz ${NETSCR}${STRAIN}
#cd ${NETSCR}${STRAIN}


# Execute commands

# Run bowtie2 to align fastq files to the reference genome
bowtie2 --seed 123 -x ${REFGENEPATH} -p 8 -U ${STRAIN}.fastq.gz -S ./Sam/${STRAIN}.sam 

# Convert sam file to a bam file
samtools view -@ 4 -b ./Sam/${STRAIN}.sam > ./Bam/${STRAIN}.bam && rm ./Sam/${STRAIN}.sam

# Only have alignments that have a mapq score greater than 5
samtools view -@ 4 -bq 5 ./Bam/${STRAIN}.bam > ./Bam/${STRAIN}_q5.bam

# Sort the bam file with mapq greater than 5
samtools sort -@ 4 -o ./Bam/${STRAIN}_q5_sorted.bam ./Bam/${STRAIN}_q5.bam
 
# Mark the reads that are PCR Duplicates
java -Xmx4g -jar ${picardPath} MarkDuplicates INPUT=./Bam/${STRAIN}_q5_sorted.bam OUTPUT=./Bam/${STRAIN}_q5_sorted_dupsMarked.bam METRICS_FILE= ./PCRdups/${STRAIN}_PCR_duplicates REMOVE_DUPLICATES= false ASSUME_SORTED= true

# Remove the reads that are marked as pcr duplicates from the bam file using the bit flag for pcr dups
samtools view -@ 4 -bF 0x400 ./Bam/${STRAIN}_q5_sorted_dupsMarked.bam > ./Bam/${STRAIN}_q5_sorted_dupsRemoved.bam

# Create index for bam file. This is needed for the next step to remove Chrm Y, U and Het
samtools index ./Bam/${STRAIN}_q5_sorted_dupsRemoved.bam

# Extract only reads from Chrm 2R,2L,3R,3L,4 and X
samtools view -@ 4 -b ./Bam/${STRAIN}_q5_sorted_dupsRemoved.bam chr2L chr2R chr3L chr3R chr4 chrX > ./Bam/${STRAIN}_q5_sorted_dupsRemoved_noYUHet.bam

# Create new index for filtered bam file. Needed to convert bam file to bed file
samtools index ./Bam/${STRAIN}_q5_sorted_dupsRemoved_noYUHet.bam

# Convert the bam file into a bed file
bedtools bamtobed -i ./Bam/${STRAIN}_q5_sorted_dupsRemoved_noYUHet.bam > ./Bam/${STRAIN}_q5_sorted_dupsRemoved_noYUHet.bed

# Bam Coverage to output bigwig file normalized to genomic coverage
bamCoverage -b ./Bam/${STRAIN}_q5_sorted_dupsRemoved_noYUHet.bam --numberOfProcessors max --normalizeTo1x 121400000 --outFileFormat bedgraph --binSize 10 -e 125 -o ./BigWigs/${STRAIN}_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC.wig

# Z-Normalize Bigwig Files and write stats to a stats file
python2.7 ./${PIPEPATH}/zNormV3.py ./BigWigs/${STRAIN}_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC.wig ./BigWigs/ZNormalized/${STRAIN}_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC_zNorm.wig > ./Stats/${STRAIN}_zNormStats.csv

# Convert z-normalized wig file to bigwig and remove wig file
wigToBigWig ./BigWigs/ZNormalized/${STRAIN}_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC_zNorm.wig ${PIPEPATH}dm3.chrom.sizes ./BigWigs/ZNormalized/${STRAIN}_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC_zNorm.bw
rm ./BigWigs/ZNormalized/${STRAIN}_q5_sorted_dupsRemoved_noYUHet_normalizedToRPGC_zNorm.wig

# Write collected ZNorm Statfile
zStatFiles=(ls ./Stats/*zNormStats.csv)
cat zStatFiles[1] > collected_zNorm_statfiles.csv
for ((i=2; i<${#zStatFiles[@]}; i++)); do
    sed -n '2,$p' < ${zStatFiles[$i]} >> collected_zNorm_statfiles.csv
done

# Create collected flagstats files
for BAM in \$(ls ./Bam/${STRAIN}*.bam | cut -d. -f1); do
	echo "\${BAM}: " >> ../Stats/${STRAIN}_flagstats.txt
	samtools flagstat \${BAM}.bam | grep -v '^0 + 0' >> ../Stats/${STRAIN}_flagstats.txt;
	echo >> ../Stats/${STRAIN}_flagstats.txt
done

# Parse FlagStats
python2.7 ./${PIPEPATH}parseflagstats.py ./Stats/

">>processFAIRESeqReadsAndCallMACSPeaks_${STRAIN}.bsub
fi

if [ "${PEAK}" = "CallPeaks" ]; then

echo "
# Change directory
#cd ${NETSCR}${STRAIN}

#Set variable for control
CONTROL=${CTRLPATH}

# Call Peaks using MACS
macs2 callpeak -t ./Bed/${STRAIN}_q5_sorted_dupsRemoved_noYUHet.bed -c \${CONTROL} -n ${STRAIN}_q5_sorted_dupsRemoved_noYUHet -g dm --nomodel --extsize 125 --seed 123 --outdir ./Peakfiles/

#First sort peak file by q-value in decreasing order, then cut out chromosome, start and end coordinates, peak name and q-value and take the top 5000 peaks based on q-value
#sort -n -r -k9 ${STRAIN}_q5_sorted_dupsRemoved_noYUHet_peaks.narrowPeak | cut -f1,2,3,4,9 | head -5000 > ${STRAIN}_q5_sorted_dupsRemoved_noYUHet_Top5000Peaks.bed

">>processFAIRESeqReadsAndCallMACSPeaks_${STRAIN}.bsub
fi

bsub < processFAIRESeqReadsAndCallMACSPeaks_${STRAIN}.bsub
