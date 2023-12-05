library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)

# Run GO enrichment analysis using GSEA procedure
setwd("~/Projects/Carol_ATACseq_Nov2023/diffbind_results/")

res <- read.csv("KO_vs_WT_results_annotated.csv", header = T)
head(res)

# Perform GO and pathway enrichment analysis 
# This will be done on differential hotspots located 
# no further then 5000 bp from 
# the middle of the promoter
res <- res[abs(res$distancetoFeature) <= 5000,]
dim(res)
head(res)
res <- res[!is.na(res$entrez_id),]
head(res$Fold)
res <- res[!duplicated(res$entrez_id),]
gsea_in <- res$Fold 
names(gsea_in) <- res$entrez_id
head(gsea_in)
gsea_in <- sort(gsea_in, decreasing = T)
head(gsea_in)
head(names(gsea_in))

# Run GSEA against GO terms
gse <- gseGO(geneList = gsea_in, ont="BP", keyType = "ENTREZID",
             minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1,
             pAdjustMethod = "BH", by = "fgsea",
             OrgDb = org.Mm.eg.db)
head(gse@result)
write.csv(gse@result, file="GSEA_BP_FC.csv", row.names = F)

# Run GSEA against GO terms (CC)
gse <- gseGO(geneList = gsea_in, ont="CC", keyType = "ENTREZID",
             minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, 
             pAdjustMethod = "BH", by = "fgsea",
             OrgDb = org.Mm.eg.db)
head(gse@result)
write.csv(gse@result, file="GSEA_CC_FC.csv", row.names = F)

# Run GSEA against GO terms (MF)
gse <- gseGO(geneList = gsea_in, ont="MF", keyType = "ENTREZID",
             minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, 
             pAdjustMethod = "BH", by = "fgsea",
             OrgDb = org.Mm.eg.db)
head(gse@result)
write.csv(gse@result, file="GSEA_MF_FC.csv", row.names = F)

# Run GSEA against KEGG pathways
kk <- gseKEGG(geneList = gsea_in, organism = "mmu", nPerm = 1000,
              keyType = "kegg", minGSSize = 10, maxGSSize = 500,
              pAdjustMethod = "BH", pvalueCutoff = 1)
write.csv(kk@result, file="GSEA_KEGG_FC.csv", row.names = F)

# Run GSEA against Reactome pathways
react <- gsePathway(geneList = gsea_in, 
                    pvalueCutoff = 1,
                    pAdjustMethod = "BH",
                    organism = "mouse",
                    minGSSize = 10,
                    maxGSSize = 500,
                    by = "fgsea",
                    verbose = FALSE)
write.csv(react@result, file="GSEA_REACTOME_FC.csv", row.names = F)









