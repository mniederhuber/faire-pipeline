# Merge Peak Set by CUyehara - v2.0
    # Merges peaksets for some set of files based on some fraction of overlap
    # Merges ALL .bed files in folder (cannot input specific files (yet))
    # Assumes either/or overlap - i.e. if either peak A or peak B overlap by some fraction the peak is merged

# To Run: Copy script into folder with bed files to be merged. 
# Usage: <fractionOverlap> <newFilename>
# Example: bash mergePeakSet.sh A.bed B.bed 0.7 a_b_union_pt7

module bash purge
module bash load bedtools

# frac = fraction overlap
frac=$1

# mergedPeaks = filename for merged peak set
mergedPeaks=$2

echo "

#BSUB -q day

#BSUB -o ${mergedPeaks}_outFile.%J

#BSUB -e ${mergedPeaks}_errFile.%J

# get a list of all bed files in current directory
fileList=(ls *bed)

# copy sorted peaks from 1st dataset to merged peak file
cat \${fileList[1]} | sort -k 1,1 -k 2,2n | cut -f 1,2,3 >> ${mergedPeaks}

# count number of peaks in the initial file
peakCount=\$(cat ${mergedPeaks} | wc -l)

# Create a stat file w/ starting peaks
echo 'Name,Num_Peaks,non_intersecting_peaks,intersecting_peaks,Total_merged_peaks' > merged_peak_stats.csv
echo "\${fileList[1]},\${peakCount},NA,NA,\${peakCount}" >> merged_peak_stats.csv

# Create a peak file to keep the non-intersecting peaks
touch non_intersecting_peaks.bed

# iterate through file list
for ((i=2; i<\${#fileList[@]}; i++)); do
    
    ## Count number of non-intersecting peaks before adding additional peaks
    nonIntersectBefore=\$(cat non_intersecting_peaks.bed | wc -l)
    
    ## get basename of file
    baseName=\$(echo "\${fileList[\$i]}" | cut -d . -f 1)
    
    ## sort first peak_file
    cat \${fileList[\$i]} | sort -k 1,1 -k 2,2n | cut -f 1,2,3 > \${baseName}_sorted.bed
    
    ## get all peaks that DO overlap by some frac w/ merged peak file
    bedtools intersect -wa -e -f ${frac} -F ${frac} -a \${baseName}_sorted.bed -b ${mergedPeaks} > temp
    
    ## get peaks that do NOT intersect with any peaks in merged file
    bedtools intersect -v -wa -e -f ${frac} -F ${frac} -a \${baseName}_sorted.bed -b ${mergedPeaks} >> non_intersecting_peaks.bed
    
    ## Combine the newly intersected peaks w/ old merged peaks, sort, and store in a temp file
    cat temp ${mergedPeaks} | sort -k 1,1 -k 2,2n | cut -f 1,2,3 > temp2
    
    ## overwrite old merged peak file w/ new, sorted peak list
    cat temp2 > ${mergedPeaks}
    
    ## count number of non-intersecting peaks in file after adding new peaks
    nonIntersectAfter=\$(cat non_intersecting_peaks.bed | wc -l)

    ## delta = number of non-intersecting peaks in sample
    #echo "\${nonIntersectAfter} \${noneIntersectBefore} HERE"
    delta=\$((\${nonIntersectAfter}-\${nonIntersectBefore}))
    
    ## intersect = total peaks in sample - number that do NOT intersect = number of intersecting peaks in sample
    sampPeaks=\$(cat \${fileList[\$i]} | wc -l)
    intersect=\$((\${sampPeaks}-\${delta}))
    
    #echo \$sampPeaks
    #echo \$intersect    

    ## peakCount = running counter of total peaks = previous peakCount + delta (non-intersecting peaks)
    peakCount=\$((\${peakCount}+\${delta}))
    
    ## echo stats to stat file
    echo -e "\${fileList[\$i]},\${sampPeaks},\${delta},\${intersect},\${peakCount}" >> merged_peak_stats.csv

done

## merge all peaks in mergedPeak file
bedtools merge -d -1 -i ${mergedPeaks} > temp2

## combine non-intersecting peak file and merged peak file
cat non_intersecting_peaks.bed temp2 | sort -k 1,1 -k 2,2n > ${mergedPeaks}

## remove temp files
rm temp
rm temp2

" > merge_peaks.bsub

bsub < merge_peaks.bsub


# PSEUDOCODE
# for each bed files in folder:
    # get basename
    # sort bed file by position, cut fields 1, 2, 3 (chr, start, end)
 # file_list = (ls *bed_pattern)
 # concatenate \${file_list[1]} >> merged_peak set
 # for ((i=1; i<\${file_list[@]}; i++))
    # report all entries of file_list[1] that intersect w/ merged_peak > temp
    # concatenate all entries of merged_peak w/ file_list[i] intersects from above >> merged_peak
    # concatenate all non-intersecting regions of file_list[i] >> non-intersect_peak
    # count lines of non-intersect peak >> stat_file
        # Name    non-intersect_peaks   intersecting peaks     total_peaks (merged_peaks + non_intersect peaks)
        # samp1   1000                  6000                   6500
        # samp2   500                   
        # etc     etc                   etc
 # Merge merge_peak
 # cat non-intersect_peak >> merge_peak
 # sort merge_peak
 # 
#

#######################
# mergePeakSet.sh v1.0
########################

#module bash purge
#module bash load bedtools
#
#ABED=\$1 # First Bed file
#BBED=\$2 # Second Bed file (order of bed files does NOT matter)
#FRAC=\$3 # Fractional threshold for peak merging. Eg 0.7 will merge peaks if overlap is 0.7
#NEWFILE=\$4 # New filename
#
## Remove file extension from new filename (if there is a file extension)
#NEWFILE=\$(echo \${NEWFILE} | cut -d . -f 1)
#
#ABEDBASE=\$(echo "\${ABED}" | cut -d . -f 1) # basename of ABED (strips extension)
#BBEDBASE=\$(echo "\${BBED}" | cut -d . -f 1) # basename of BBED (strips extension)
#
##remove this comment # echo "
#
##BSUB -q day
#
##BSUB -o \${NEWFILE}_outfile.%J
#
##BSUB -e \${NEWFILE}_errorfile.%J
#
## Reformat and sort peaks lists
#cat \${ABED} | sort -k 1,1 -k 2,2n | cut -f 1,2,3 > \${ABEDBASE}_sorted.bed
#cat \${BBED} | sort -k 1,1 -k 2,2n | cut -f 1,2,3 > \${BBEDBASE}_sorted.bed
#
## Report all entries of A that intersect with B by some fraction
#bedtools intersect -wa -e -f \${FRAC} -F \${FRAC} -a \${ABEDBASE}_sorted.bed -b \${BBEDBASE}_sorted.bed > temp
#
## Concatenate all entries of B with regions of A that intersect and then merge all entries
#cat temp \${BBEDBASE}_sorted.bed | sort -k 1,1 -k 2,2n | bedtools merge -d -1 -i - > temp2
#
## Concatenate the non-intersecting entries of A to merged file
#bedtools intersect -v -wa -e -f \${FRAC} -F \${FRAC} -a \${ABEDBASE}_sorted.bed -b \${BBEDBASE}_sorted.bed >> temp2
#
## Sort file and write to a new file
#cat temp2 | sort -k 1,1 -k 2,2n > \${NEWFILE}.bed
#
## Remove temporary files
#rm temp temp2
#" > \${NEWFILE}.bsub
#
#bsub < \${NEWFILE}.bsub
