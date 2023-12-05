# ATACseq_worflow_Nov2023
Updated ATAC-seq workflow and accompanying scripts as of November 2023

### Workflow

1. Merge fastq files from different lanes

2. Run FastQC on merged fastq files

```
fastqc --threads 8 *.fastq.gz
``` 

3. Trim adapters with trim galore and produce FastQC report on trimmed reads

```
trim_galore --paired -j 8 --nextera --fastqc -q 30 <R1.fastq.gz> <R2.fastq.gz>

```

4. Map the reads to the reference genome with bowtie2

```
bowtie2 --threads 6 --very-sensitive -x <bowtie2_index> -1 <R1.fq.gz> -2 <R2.fq.gz> -S <sam>
```

5. Convert sam to bam, sort and index

```
samtools view -bS <SAM> > <BAM>
samtools sort <BAM> <sorted>
samtools index <sorted.bam>

```

6. Remove mitochondrial reads

```
time samtools idxstats --threads 6 <BAM> | cut -f 1 | grep -v MT | xargs samtools view --threads 6 -b <BAM> > <NO_MT.BAM>
time samtools index --threads 6 <NO_MT.BAM>

```

7. Call peaks with MACS2

```
macs2 callpeak -t <BAM> -f BAMPE -n <sample_name> -g mm -q 0.05 --keep-dup 1 --outdir <out_dir>

```

8. Exclude blacklisted (repeat) rich regions

Link to blacklisted regions: https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz

```
bedtools intersect  -a  <broadPeak>  -b <BLACKLIST.BED>
```

9. Run quality control analysis on initial and mitochondria filtered bam files using  

```
qualimap bamqc --java-mem-size=4G  -nt 6  -bam  <BAM> -outdir  <OUT>  -outfile  <OUT>
```
10. Use MACS2 output and bam files filtered of mitochondrial reads to compare peaks in DiffBind

11. Convert filtered bam files to normalized bigwigs

```
bamCoverage --bam <bam> -p 6 -o <out.bw> --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2150570000 --ignoreForNormalization chrX --extendReads

```

