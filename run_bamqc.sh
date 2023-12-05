#! /bin/bash

# Run bamqc 

BAM_DIR=$1 # directory with bam files

for file in $BAM_DIR/*.bam
do 
    filename=`basename $file`
    samplename=`echo ${filename%%.*}`
    echo qualimap bamqc --java-mem-size=4G -nt 6   -bam  $file  -outdir  $samplename   -outfile   $samplename
    qualimap bamqc --java-mem-size=4G -nt 6   -bam  $file  -outdir  $samplename   -outfile   $samplename

done