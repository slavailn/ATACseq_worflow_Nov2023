library(ATACseqQC)

setwd("~/Projects/Carol_ATACseq_Nov2023/")

# Plot fragment size distributions
pdf("MBT1KO1_fragment_size_dist.pdf", width = 6, height = 6)
fragSizeDist("no_mt_bam/MBT1KO1.sorted.no_mt.bam", 
             "MBT1KO1")
dev.off()

pdf("MBT1KO2_fragment_size_dist.pdf", width = 6, height = 6)
fragSizeDist("no_mt_bam/MBT1KO2.sorted.no_mt.bam", 
             "MBT1KO2")
dev.off()

pdf("MBT1KO3_fragment_size_dist.pdf", width = 6, height = 6)
fragSizeDist("no_mt_bam/MBT1KO3.sorted.no_mt.bam", 
             "MBT1KO3")
dev.off()

pdf("WT1_fragment_size_dist.pdf", width = 6, height = 6)
fragSizeDist("no_mt_bam/WT1.sorted.no_mt.bam", 
             "WT1")
dev.off()

pdf("WT2_fragment_size_dist.pdf", width = 6, height = 6)
fragSizeDist("no_mt_bam/WT2.sorted.no_mt.bam", 
             "WT2")
dev.off()

pdf("WT3_fragment_size_dist.pdf", width = 6, height = 6)
fragSizeDist("no_mt_bam/WT3.sorted.no_mt.bam", 
             "WT3")
dev.off()





