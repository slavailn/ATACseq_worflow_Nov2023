#! /bin/bash

# Map reads to the reference genome with bowtie2

FASTQ_DIR=$1 # directory containing Read1 and Read2 *.fastq.gz files
REF_INDEX=$2 # bowtie2 index directory

# Create empty array for read1 and read2 files
list_files_r1=()
list_files_r2=()

# populate read1 and read2 array
for file in $FASTQ_DIR/*.fq.gz
do
    if [[ "$file" =~ ".R1_" ]]
    then
        list_files_r1+=($file)
    else
        list_files_r2+=($file)
    fi
done

echo 'Bowtie2 commands '
echo '##########################################################'

# Iterate over the arrays and extract read1 and read2 files as arguments to trim galore command
for i in $( seq 0 $(( ${#list_files_r1[*]}-1 )) )
do
    file_name=`basename ${list_files_r1[$i]}`
    sample_name=${file_name%%.*}
    bowtie2 --threads 6 --very-sensitive -x $REF_INDEX -1 "${list_files_r1[$i]}" -2 "${list_files_r2[$i]}" -S $sample_name.sam
    echo  bowtie2 --threads 6 --very-sensitive -x $REF_INDEX -1 "${list_files_r1[$i]}" -2 "${list_files_r2[$i]}" -S $sample_name.sam
done

echo ''
echo '##########################################################'