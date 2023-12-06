library(ATACseqQC)

setwd(<WORKING_DIR>)

# Plot fragment size distributions
# Use this code on every bam files
pdf("fragment_size_dist.pdf", width = 6, height = 6)
fragSizeDist("path/to/bam", "SampleName")
dev.off()



