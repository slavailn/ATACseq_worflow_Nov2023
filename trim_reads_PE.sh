#! /bin/bash

FASTQ_DIR=$1 # directory containing Read1 and Read2 *.fastq.gz files

# Create empty array for read1 and read2 files
list_files_r1=()
list_files_r2=()

# populate read1 and read2 array
for file in $FASTQ_DIR/*.fastq.gz
do
    if [[ "$file" =~ ".R1." ]]
    then
        list_files_r1+=($file)
    else
        list_files_r2+=($file)
    fi
done

echo 'Trime Galore! commands '
echo '##########################################################'

# Iterate over the arrays and extract read1 and read2 files as arguments to trim galore command
for i in $( seq 0 $(( ${#list_files_r1[*]}-1 )) )
do 
    trim_galore --paired -j 8 --nextera --fastqc -q 30 "${list_files_r1[$i]}" "${list_files_r2[$i]}"
    echo trim_galore --paired -j 8 --nextera --fastqc -q 30 "${list_files_r1[$i]}" "${list_files_r2[$i]}"
done

echo ''
echo '##########################################################'