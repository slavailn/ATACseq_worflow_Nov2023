library(DiffBind)

# Try different methods of normalization and library size estimation

setwd(<WORKING_DIR>)

# Load DiffBind object with read counts
load("ATAC_diffind_counts.DBA.RData")
atac

# Quick comparison of normalization methods
atac <- dba.normalize(atac, method=DBA_ALL_METHODS,
                             normalize=DBA_NORM_NATIVE,
                             background=TRUE)
atac <- dba.analyze(atac, method=DBA_ALL_METHODS)
dba.show(atac,bContrasts=TRUE)
par(mfrow=c(2,1))
dba.plotMA(atac, method=DBA_EDGER, sub="edgeR:TMM:background")
dba.plotMA(atac, method=DBA_DESEQ2, sub="DESeq2:RLE:background")
par(mfrow=c(1,1))

# Different methods to run DiffBind comparisons, the authors of DiffBind
# recommend using background to calculate library sizes
# 1. method = DBA_EGDER, normalzation = TMM, library = DBA_LIBSIZE_FULL
# Result 2 regions FDR < 0.05
load("ATAC_diffind_counts.DBA.RData")
atac <- dba.normalize(atac, method = DBA_EDGER,
                      normalize = DBA_NORM_TMM,
                      library = DBA_LIBSIZE_FULL,
                      background = F)
atac <- dba.analyze(atac, method = DBA_EDGER)
dba.show(atac,bContrasts=TRUE)
pdf("edgeR_TMM_DBA_LIBSIZE_FULL_MAplot.pdf")
dba.plotMA(atac, method=DBA_EDGER, sub="edgeR:TMM:DBA_LIBSIZE_FULL")
dev.off()

# 2. method = DBA_EGDER, normalization = TMM, library = DBA_LIBSIZE_FULL
# Result 2 regions FDR < 0.05
load("ATAC_diffind_counts.DBA.RData")
atac <- dba.normalize(atac, method = DBA_EDGER,
                      normalize = DBA_NORM_TMM,
                      library = DBA_LIBSIZE_PEAKREADS,
                      background = F)
atac <- dba.analyze(atac, method = DBA_EDGER)
dba.show(atac,bContrasts=TRUE)
pdf("edgeR_TMM_DBA_LIBSIZE_PEAKREADS_MAplot.pdf")
dba.plotMA(atac, method=DBA_EDGER, sub="edgeR:TMM:DBA_LIBSIZE_PEAKREADS")
dev.off()

# 3. method = DBA_EGDER, normalization = TMM, library = DBA_LIBSIZE_BACKGROUND
# Result 3 regions FDR < 0.05
load("ATAC_diffind_counts.DBA.RData")
atac <- dba.normalize(atac, method = DBA_EDGER,
                      normalize = DBA_NORM_TMM,
                      library = DBA_LIBSIZE_BACKGROUND,
                      background = T)
atac <- dba.analyze(atac, method = DBA_EDGER)
dba.show(atac,bContrasts=TRUE)
pdf("edgeR_TMM_DBA_LIBSIZE_BACKGROUND_MAplot.pdf")
dba.plotMA(atac, method=DBA_EDGER, sub="edgeR:TMM:DBA_LIBSIZE_PEAKREADS")
dev.off()

# 4. method = DBA_DESEQ2, normalization = RLE, library = DBA_LIBSIZE_FULL
# Result 2 regions FDR < 0.05
load("ATAC_diffind_counts.DBA.RData")
atac <- dba.normalize(atac, method = DBA_DESEQ2,
                      normalize = DBA_NORM_RLE,
                      library = DBA_LIBSIZE_FULL,
                      background = F)
atac <- dba.analyze(atac, method = DBA_DESEQ2)
dba.show(atac,bContrasts=TRUE)
pdf("DESEQ2_DBA_LIBSIZE_FULL_MAplot.pdf")
dba.plotMA(atac, method=DBA_DESEQ2, sub="edgeR:RLE:DBA_LIBSIZE_FULL")
dev.off()

# 5. method = DBA_DESEQ2, normalization = RLE, library = DBA_LIBSIZE_PEAKREADS
# Result 10 regions FDR < 0.05
load("ATAC_diffind_counts.DBA.RData")
atac <- dba.normalize(atac, method = DBA_DESEQ2,
                      normalize = DBA_NORM_RLE,
                      library = DBA_LIBSIZE_PEAKREADS,
                      background = F)
atac <- dba.analyze(atac, method = DBA_DESEQ2)
dba.show(atac,bContrasts=TRUE)
pdf("DESEQ2_DBA_LIBSIZEL_PEAKREADS_MAplot.pdf")
dba.plotMA(atac, method=DBA_DESEQ2, sub="edgeR:RLE:DBA_PEAKREADS")
dev.off()

# 6. method = DBA_DESEQ2, normalization = RLE, library = DBA_LIBSIZE_BACKGROUND
# Result 10 regions FDR < 0.05
load("ATAC_diffind_counts.DBA.RData")
atac <- dba.normalize(atac, method = DBA_DESEQ2,
                      normalize = DBA_NORM_RLE,
                      library = DBA_LIBSIZE_BACKGROUND,
                      background = T)
atac <- dba.analyze(atac, method = DBA_DESEQ2)
dba.show(atac,bContrasts=TRUE)
pdf("DESEQ2_DBA_LIBSIZE_BACKGROUND_MAplot.pdf")
dba.plotMA(atac, method=DBA_DESEQ2, sub="DESEQ2:RLE:DBA_BACKGROUND")
dev.off()






















