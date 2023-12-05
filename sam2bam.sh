#! /bin/bash

BAM_DIR=$1 # directory containing bam files
OUT_DIR=$2 # Directory to save output in

# This script assumes that sam files were already converted to bam

for file in $BAM_DIR/*.bam
do
    filename=`basename $file`
    samplename=`echo ${filename%%.bam}`
    echo time samtools sort --threads 6 $file -o $OUT_DIR/$samplename.sorted.bam
    echo time samtools index --threads 6 $OUT_DIR/$samplename.sorted.bam 
    time samtools sort --threads 6 $file -o $OUT_DIR/$samplename.sorted.bam
    time samtools index --threads 6 $OUT_DIR/$samplename.sorted.bam 
    
done

