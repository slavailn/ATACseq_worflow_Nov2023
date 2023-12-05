#! /bin/bash

# Parse and summarize the output for samtools flagstat command
# for paired-end reads

BAM_DIR=$1 # directory with bam files

echo -e "SampleID\ttotal_reads\tmapped_reads\tproperly_paired"

for file in $BAM_DIR/*.bam
do
    filename=`basename $file`
    samplename=${filename%%.*}

    stats[0]=$samplename
      
    samtools flagstat $file > stats.tmp 
    while IFS= read -r line
    do
        if [[ $line =~ "in total"  ]]
        then
             total_reads=`echo $line | cut -f 1 -d ' '`
             stats[1]="$total_reads"
        fi

        if [[ $line =~ 'mapped (' ]]
        then
             mapped=`echo $line | cut -f 1 -d ' '`
             stats[2]="$mapped"
        fi

        if [[ $line =~ "properly paired" ]]
        then
             properly_paired=`echo $line | cut -f 1 -d ' '`
             stats[3]="$properly_paired"
        fi
    done < "stats.tmp"

    echo -e "${stats[0]}\t${stats[1]}\t${stats[2]}\t${stats[3]}"

done

rm stats.tmp