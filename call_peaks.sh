#! /bin/bash

# Call peaks with macs2

BAM_DIR=$1 # directory containing bam files

# Iterate over bam files and call peaks with macs2
for file in $BAM_DIR/*.bam
do
    filename=`basename $file`
    samplename=`echo ${filename%%.bam}`
    echo macs2 callpeak -t $file -f BAMPE -n $samplename -g mm --broad --broad-cutoff 0.05 --keep-dup 1 --outdir $samplename
    time macs2 callpeak -t $file -f BAMPE -n $samplename --broad --broad-cutoff 0.05 -g mm --keep-dup 1 --outdir $samplename
done

