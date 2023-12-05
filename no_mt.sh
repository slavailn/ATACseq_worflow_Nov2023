#! /bin/bash

# Remove mitochondrial aligments and re-index bam files

BAM_DIR=$1 # directory containing bam files
OUT_DIR=$2 # Directory to save output in

# This script removes mitochondrial reads

for file in $BAM_DIR/*.bam
do
    filename=`basename $file`
    samplename=`echo ${filename%%.bam}`
    echo time samtools idxstats --threads 6 $file | cut -f 1 | grep -v MT | xargs samtools view --threads 6 -b $file > $samplename.no_mt.bam
    echo time samtools index --threads 6 $OUT_DIR/$samplename.no_mt.bam  
    time samtools idxstats --threads 6 $file | cut -f 1 | grep -v MT | xargs samtools view --threads 6 -b $file > $samplename.no_mt.bam
    time samtools index --threads 6 $OUT_DIR/$samplename.no_mt.bam
    
done

