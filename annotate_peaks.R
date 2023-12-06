library(ChIPpeakAnno)
library(org.Mm.eg.db)
data("TSS.mouse.GRCm38")

# Annotate DiffBind results using ChIPPeakAnno

setwd(<WORKING_DIR>)

# Load GRanges object that stores the results
# of DiffBind analysis
load("KO_vs_WT.GRanges.RData")
ls()
dba.object.DB

# Convert DiffBind genomic ranges object to a data frame
diff <- as.data.frame(dba.object.DB)
head(diff)
write.csv(diff, "diffbind_results.csv")

# Annotate the results of DiffBind analysis
# annotate by promoters
promoterData <- promoters(TSS.mouse.GRCm38, upstream=5000, downstream=500)
dba.object.DB_annot <- annotatePeakInBatch(dba.object.DB, 
                                           AnnotationData=promoterData)
# Add gene additional gene info
dba.object.DB_annot <- addGeneIDs(dba.object.DB_annot, orgAnn="org.Mm.eg.db", 
                                  IDs2Add=c("symbol", "entrez_id"))
head(dba.object.DB_annot)
write.csv(as.data.frame(dba.object.DB_annot), 
          file="KO_vs_WT_results_annotated.csv" )




