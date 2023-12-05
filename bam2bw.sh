#! /bin/bash

# Convert generate normalized bigwig genome coverage files based on alignments in BAM format
BAM_DIR=$1 # directory with bam files

for file in $BAM_DIR/*.bam
do 
    filename=`basename $file`
    samplename=`echo ${filename%%.*}`
    echo bamCoverage --bam $file -p 6 -o $samplename.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2150570000 --ignoreForNormalization chrX --extendReads
    time bamCoverage --bam $file -p 6 -o $samplename.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2150570000 --ignoreForNormalization chrX --extendReads

done


