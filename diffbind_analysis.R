library(DiffBind)

setwd(<WORKING_DIR>)

# Create target file
SampleID <- c("KO1", "KO2", "KO3", "WT1", "WT2", "WT3")
Tissue <- rep("NA", 6) 
Factor <- rep("NA", 6)
Condition <- rep("NA", 6)
Treatment <- c("ko", "ko", "ko", "wt", "wt", "wt")
Replicate <- c(1,2,3,1,2,3)

bamReads <- c("no_mt_bam/KO1.sorted.no_mt.bam", 
              "no_mt_bam/KO2.sorted.no_mt.bam",
              "no_mt_bam/KO3.sorted.no_mt.bam", 
              "no_mt_bam/WT1.sorted.no_mt.bam",
              "no_mt_bam/WT2.sorted.no_mt.bam", 
              "no_mt_bam/WT3.sorted.no_mt.bam")

Peaks <- c("macs2_peaks/KO1.sorted.no_mt/MBT1KO1.sorted.no_mt_peaks.xls", 
           "macs2_peaks/KO2.sorted.no_mt/MBT1KO2.sorted.no_mt_peaks.xls",
           "macs2_peaks/KO3.sorted.no_mt/MBT1KO3.sorted.no_mt_peaks.xls", 
           "macs2_peaks/WT1.sorted.no_mt/WT1.sorted.no_mt_peaks.xls",
           "macs2_peaks/WT2.sorted.no_mt/WT2.sorted.no_mt_peaks.xls", 
           "macs2_peaks/WT3.sorted.no_mt/WT3.sorted.no_mt_peaks.xls")

PeakCaller <- rep("macs", 6)

samples <- data.frame(SampleID=SampleID, Tissue=Tissue, Factor=Factor,
                      Condition=Condition, Treatment=Treatment, Replicate=Replicate,
                      bamReads=bamReads, Peaks=Peaks, PeakCaller=PeakCaller)
samples

write.csv(samples, file="samples.csv")

config = data.frame(AnalysisMethod=DBA_DESEQ2, th=0.05, 
                    DataType=DBA_DATA_GRANGES, RunParallel=T,
                    bCorPlot=F, bUsePval=F,
                    doBlacklist=F, doGreyList=F)

# Read in hotspot data and create dba object
atac <- dba(sampleSheet=samples, minOverlap = 2, config = config)
save(atac, file="atac_DBA.RData")
tiff("Correlation_heatmap_using_occupancy_peakCaller_data.tiff", 
     width = 800, height = 800)
plot(atac)
dev.off()

# Obtain read counts for hotspots
atac <- dba.count(atac, minOverlap=2)
save(atac, file="ATAC_diffind_counts.DBA.RData")

# Correlation plot with count data
tiff("Correlation_heatmap_using_counts_data.tiff", width = 800, height = 800)
plot(atac)
dev.off()

# Show the number of reads overlapping peaks
info <- dba.show(atac)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP, 
                  PeakReads=round(info$Reads * info$FRiP))
write.csv(libsizes, file="FRiP.csv")

# Plot PCA
tiff("ATACseq_PCA_based_on_read_counts.tiff", width = 800, height = 800)
dba.plotPCA(atac, DBA_TREATMENT, label=DBA_ID)
dev.off()

# Get masks for specific experimental groups
Mutmask <- dba.mask(atac, DBA_TREATMENT, "ko")
WTmask <- dba.mask(atac, DBA_TREATMENT, "wt")

dba.object <- atac
group1 <- Mutmask
group2 <- WTmask
name1 <- "ko"
name2 <- "wt"

# Detect differentially accesible regions
dba.object <- dba.contrast(atac, group1=group1, group2=group2, name1=name1, name2=name2, minMembers=2)
dba.object <- dba.analyze(dba.object, bParallel = F)
dba.object.DB <- dba.report(dba.object,th=1)

# Save the results of differential analysis
save(dba.object.DB, file="KO_vs_WT.GRanges.RData")

# Plot MA
pdf("KO_vs_WT_MAplot.pdf", width=6, height=6)
dba.plotMA(dba.object)
dev.off()

# Plot PCA
pdf("KO_vs_WT_PCAplot.pdf", width=6, height=6)
dba.plotPCA(dba.object, contrast=1, th=0.2, label=DBA_TREATMENT)
dev.off()


































