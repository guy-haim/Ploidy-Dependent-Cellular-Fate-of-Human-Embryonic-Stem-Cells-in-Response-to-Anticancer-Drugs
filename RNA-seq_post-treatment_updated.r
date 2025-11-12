library(stringr)
library(dplyr)
library(edgeR)
library(ggplot2)
library(ggfortify)
library(pheatmap)
library(data.table)
library(zoo)
library(gplots)
library(PCAtools)
library(reshape2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(FSA)
library(ggsignif)
library(ggpubr)

setwd(choose.dir())

#loading RNA-seq data from each sample into a data.frame and getting rid of unnecessary data (first 4 rows and second and third columns)
Diploid_rep1_5FU_2.5uM <- read.delim("10R_4c_p48_5_FU_2_5uM_48h_S8_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]
Diploid_rep1_Bleo_20nM <- read.delim("10R_4c_p48_Bleo_20nM_48h_S7_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]
Diploid_rep1_Cis_100nM <- read.delim("10R_4c_p48_Cis_100nM_48h_S3_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]
Diploid_rep1_Pacli_3.3nM <- read.delim("10R_4c_p48_Pacli_3_3nM_48h_S4_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]

Diploid_rep2_5FU_2.5uM <- read.delim("10R_4c_p48_5_FU_2_5uM_48h_S9_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]
Diploid_rep2_Bleo_20nM <- read.delim("10R_4c_p48_Bleo_20nM_48h_S10_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]
Diploid_rep2_Cis_100nM <- read.delim("10R_4c_p48_Cis_100nM_48h_S11_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]
Diploid_rep2_Pacli_3.3nM <- read.delim("10R_4c_p48_Pacli_3_3nM_48h_S12_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]

Diploid_rep3_5FU_2.5uM <- read.delim("10G_4c_p47_5_FU_2_5uM_48h_S17_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]
Diploid_rep3_Bleo_20nM <- read.delim("10G_4c_p47_Bleo_20nM_48h_S18_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]
Diploid_rep3_Cis_100nM <- read.delim("10G_4c_p47_Cis_100nM_48h_S19_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]
Diploid_rep3_Pacli_3.3nM <- read.delim("10G_4c_p47_Pacli_3_3nM_48h_S20_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]

Triploid_B_5FU_2.5uM <- read.delim("3n_clone_B_f10_5_FU_2_5uM_48h_S6_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]
Triploid_B_Bleo_20nM <- read.delim("3n_clone_B_f10_Bleo_20nM_48h_S5_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]
Triploid_B_Cis_100nM <- read.delim("3n_clone_B_f9_Cis_100nM_48h_S1_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]
Triploid_B_Pacli_3.3nM <- read.delim("3n_clone_B_f9_Pacli_3_3nM_48h_S2_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]

Triploid_E_5FU_2.5uM <- read.delim("3n_clone_E_f10_5_FU_2_5uM_48h_S13_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]
Triploid_E_Bleo_20nM <- read.delim("3n_clone_E_f10_Bleo_20nM_48h_S14_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]
Triploid_E_Cis_100nM <- read.delim("3n_clone_E_f10_Cis_100nM_48h_S15_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]
Triploid_E_Pacli_3.3nM <- read.delim("3n_clone_E_f10_Pacli_3_3nM_48h_S16_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]

Triploid_J_5FU_2.5uM <- read.delim("3n_clone_J_f12_5_FU_2_5uM_48h_S21_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]
Triploid_J_Bleo_20nM <- read.delim("3n_clone_J_f12_Bleo_20nM_48h_S22_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]
Triploid_J_Cis_100nM <- read.delim("3n_clone_J_f12_Cis_100nM_48h_S23_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]
Triploid_J_Pacli_3.3nM <- read.delim("3n_clone_J_f12_Pacli_3_3nM_48h_S24_L001_R1_001_counts.txt", header = T, skip = 1)[,c(1,7)]


###merging all RNA-seq data into one data.frame and adding gene names according to annotation file:
gene_exp_all <- Reduce(function(x, y) merge(x, y, all=T, by = "Geneid"), 
                       list(Diploid_rep1_5FU_2.5uM,Diploid_rep1_Bleo_20nM,Diploid_rep1_Cis_100nM,Diploid_rep1_Pacli_3.3nM,
                            Diploid_rep2_5FU_2.5uM,Diploid_rep2_Bleo_20nM,Diploid_rep2_Cis_100nM,Diploid_rep2_Pacli_3.3nM,
                            Diploid_rep3_5FU_2.5uM,Diploid_rep3_Bleo_20nM,Diploid_rep3_Cis_100nM,Diploid_rep3_Pacli_3.3nM,
                            Triploid_B_5FU_2.5uM,Triploid_B_Bleo_20nM,Triploid_B_Cis_100nM,Triploid_B_Pacli_3.3nM,
                            Triploid_E_5FU_2.5uM,Triploid_E_Bleo_20nM,Triploid_E_Cis_100nM,Triploid_E_Pacli_3.3nM,
                            Triploid_J_5FU_2.5uM,Triploid_J_Bleo_20nM,Triploid_J_Cis_100nM,Triploid_J_Pacli_3.3nM))

colnames(gene_exp_all) <- c("ENSEMBL","Diploid_rep1_5FU_2.5uM","Diploid_rep1_Bleo_20nM","Diploid_rep1_Cis_100nM","Diploid_rep1_Pacli_3.3nM",
                            "Diploid_rep2_5FU_2.5uM","Diploid_rep2_Bleo_20nM","Diploid_rep2_Cis_100nM","Diploid_rep2_Pacli_3.3nM",
                            "Diploid_rep3_5FU_2.5uM","Diploid_rep3_Bleo_20nM","Diploid_rep3_Cis_100nM","Diploid_rep3_Pacli_3.3nM",
                            "Triploid_B_5FU_2.5uM","Triploid_B_Bleo_20nM","Triploid_B_Cis_100nM","Triploid_B_Pacli_3.3nM",
                            "Triploid_E_5FU_2.5uM","Triploid_E_Bleo_20nM","Triploid_E_Cis_100nM","Triploid_E_Pacli_3.3nM",
                            "Triploid_J_5FU_2.5uM","Triploid_J_Bleo_20nM","Triploid_J_Cis_100nM","Triploid_J_Pacli_3.3nM")
gene_lengths <- read.csv("gencode.v34.transcriptLengths.csv")[,-8]

#removing the end of the ENSEMBL ID so that the same genes from different versions of the annotation file will be joined together  
gene_exp_all$ENSEMBL <- str_remove(gene_exp_all$ENSEMBL, "[.].*")
gene_lengths$ENSEMBL <- str_remove(gene_lengths$ENSEMBL, "[.].*")
gene_exp_all <- aggregate(gene_exp_all[-1], list(gene_exp_all$ENSEMBL), FUN = sum, na.rm = T)
colnames(gene_exp_all)[1] <- "ENSEMBL"

#merging all data frame to one count table with all the data needed:
gene_exp_all <- merge(gene_lengths, gene_exp_all, by = "ENSEMBL", all.y = T)
gene_exp_all <- gene_exp_all[rowSums(is.na(gene_exp_all[ ,2:31])) == 0, ] #removing empty rows
indRemoved <- which(apply(gene_exp_all[,8:31], 1, function(x) all(x == 0)) ) #saving all unexpressed genes in the tpm data.frame
gene_exp_all <- gene_exp_all[-indRemoved,]      #getting rid of rows (genes) that all of their columns are zeros
.rowNamesDF(gene_exp_all, make.names=T) <- gene_exp_all$SYMBOL


cl=c("blue1","blue2","blue3","blue4","deepskyblue1","deepskyblue3","royalblue1","royalblue3","dodgerblue1","dodgerblue2","dodgerblue3","dodgerblue4",
     "red1","red2","red3","red4","tomato1","tomato2","tomato3","tomato4","coral1","coral2","coral3","coral4")
names(cl) <- str_replace_all(colnames(gene_exp_all[,8:31]),"_", " ")

#### Differential expression (DE) analysis of all samples:

ploidy <- factor(c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2), labels = c("2n","3n"))
DE_object <- DGEList(gene_exp_all[,c(8:31)], group = ploidy, genes = gene_exp_all$gene)
keep <- filterByExpr(DE_object, min.count=5, group = ploidy)
DE_object <- DE_object[keep, , keep.lib.sizes=FALSE]
DE_TMM <- calcNormFactors(DE_object)

cpm_normalized <- cpm(DE_TMM, log = T) #log2(cpm) after TMM normalization
cpm_normalized <- data.frame("SYMBOL" = rownames(cpm_normalized), cpm_normalized)
cpm_normalized <- merge(gene_lengths, cpm_normalized, by = "SYMBOL",all.y = T)
.rowNamesDF(cpm_normalized, make.names=T) <- cpm_normalized$SYMBOL
indRemoved <- which(apply(cpm_normalized[,8:31], 1, function(x) all(x < 1)) ) #saving all unexpressed genes in the tpm data.frame
cpm_normalized <- cpm_normalized[-indRemoved,]      #getting rid of rows (genes) that all of their columns are zeros

write.csv(cpm_normalized, file = "cpm_normalized_all.csv")


#### Differential expression (DE) analysis of 5-FU:

ploidy <- factor(c(1,1,1,2,2,2), labels = c("2n","3n")) #3 diploid, and 3 triploid samples
DE_object <- DGEList(gene_exp_all[,c(8,12,16,20,24,28)], group = ploidy, genes = gene_exp_all$gene)
keep <- filterByExpr(DE_object, min.count=5, group = ploidy)
DE_object <- DE_object[keep, , keep.lib.sizes=FALSE]
DE_TMM <- calcNormFactors(DE_object)

###PCA analysis of log2(CPM) values of each sample
cpm_normalized_5FU <- cpm(DE_TMM, log = T) #log2(cpm) after TMM normalization
cpm_normalized_5FU <- data.frame("SYMBOL" = rownames(cpm_normalized_5FU), cpm_normalized_5FU)
cpm_normalized_5FU <- merge(gene_lengths, cpm_normalized_5FU, by = "SYMBOL",all.y = T)
.rowNamesDF(cpm_normalized_5FU, make.names=T) <- cpm_normalized_5FU$SYMBOL
indRemoved <- which(apply(cpm_normalized_5FU[,8:13], 1, function(x) all(x < 1)) ) #saving all unexpressed genes in the tpm data.frame
cpm_normalized_5FU <- cpm_normalized_5FU[-indRemoved,]      #getting rid of rows (genes) that all of their columns are zeros
write.csv(cpm_normalized_5FU, file = "cpm_normalized_5FU.csv")


#continue DE analysis
design <- model.matrix(~ploidy)
DE_TMM <- estimateDisp(DE_TMM, design, robust=TRUE)
fit_TMM <- glmFit(DE_TMM, design)   #fit the model
LRT_3n_Vs_2n_5FU <- glmTreat(fit_TMM)
summary(decideTests(LRT_3n_Vs_2n_5FU))

TMM_3n_5FU <- topTags(LRT_3n_Vs_2n_5FU, n=nrow(DE_TMM),adjust.method = "BH")$table #assigning fold change gene expression between diploid and triploid ESCs
TMM_3n_5FU <- data.frame("SYMBOL"=row.names(TMM_3n_5FU), TMM_3n_5FU)

colnames(TMM_3n_5FU)[2] <- "logFC_triploids"

TMM_3n_5FU <- merge(gene_lengths, TMM_3n_5FU, by = "SYMBOL", all.y = T)

sig_DE_3n_5FU <-TMM_3n_5FU[TMM_3n_5FU$FDR<=0.05,] #assigning only genes with significant fold change between triploid and diploid ESCs (by FDR) 
sig_DE_3n_5FU <- sig_DE_3n_5FU[-rowSums(is.na(sig_DE_3n_5FU[ ,2:7])) == 0, ] #removing empty rows
sig_DE_3n_5FU<- sig_DE_3n_5FU[order(sig_DE_3n_5FU$logFC_triploids, decreasing = F),]
write.csv(sig_DE_3n_5FU, file = "significant DE genes between 5-FU-treated diploid and triploid hESCs.csv")


############# Gene Set Enrichment Analysis 5FU ############# 
library(fgsea)
library(dplyr)

ranks <- TMM_3n_5FU
ranks[, 'score'] <-sign(ranks$logFC_triploids)* (-log(x = ranks$PValue,base = 10))
.rowNamesDF(ranks, make.names=T) <- ranks$SYMBOL
ranks <- setNames(ranks$score,ranks$SYMBOL)
ranks <- ranks[which(!duplicated(ranks))]


## GSEA of 5FU-treated triploid Vs diploid hESCs:
comb_Terms <- gmtPathways("combined_collection.gmt") 
msigdb_pathways <- comb_Terms
set.seed(42)
fgseaRes <- fgseaMultilevel(msigdb_pathways, ranks, minSize=15, maxSize=500, eps = 0)

fgseaResTidy_5FU <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy_5FU %>%
  dplyr::select(-ES) %>%
  arrange(padj) %>%
  DT::datatable()
write.csv(x = data.frame(fgseaResTidy_5FU[,-8],
                         "leadingEdge"=sapply(fgseaResTidy_5FU$leadingEdge, function(x) paste(x, collapse = ";"))),
          file = "All GSEA results of 5FU-treated triploid hESCs.csv")

write.csv(x = data.frame(fgseaResTidy_5FU[fgseaResTidy_5FU$padj<0.05,-8],
                         "leadingEdge"=sapply(fgseaResTidy_5FU[fgseaResTidy_5FU$padj<0.05,8], function(x) paste(x, collapse = ";"))),
          file = "Significant GSEA results of 5FU-treated triploid hESCs.csv")



#### Differential expression (DE) analysis of Bleo:

ploidy <- factor(c(1,1,1,2,2,2), labels = c("2n","3n")) #4 diploid, and 4 triploid samples
DE_object <- DGEList(gene_exp_all[,c(9,13,17,21,25,29)], group = ploidy, genes = gene_exp_all$gene)
keep <- filterByExpr(DE_object, min.count=5, group = ploidy)
DE_object <- DE_object[keep, , keep.lib.sizes=FALSE]
DE_TMM <- calcNormFactors(DE_object)

###PCA analysis of log2(CPM) values of each sample
cpm_normalized_Bleo <- cpm(DE_TMM, log = T) #log2(cpm) after TMM normalization
cpm_normalized_Bleo <- data.frame("SYMBOL" = rownames(cpm_normalized_Bleo), cpm_normalized_Bleo)
cpm_normalized_Bleo <- merge(gene_lengths, cpm_normalized_Bleo, by = "SYMBOL",all.y = T)
.rowNamesDF(cpm_normalized_Bleo, make.names=T) <- cpm_normalized_Bleo$SYMBOL
indRemoved <- which(apply(cpm_normalized_Bleo[,8:13], 1, function(x) all(x < 1)) ) #saving all unexpressed genes in the tpm data.frame
cpm_normalized_Bleo <- cpm_normalized_Bleo[-indRemoved,]      #getting rid of rows (genes) that all of their columns are zeros
write.csv(cpm_normalized_Bleo, file = "cpm_normalized_Bleo.csv")


#continue DE analysis
design <- model.matrix(~ploidy)
DE_TMM <- estimateDisp(DE_TMM, design, robust=TRUE)
fit_TMM <- glmFit(DE_TMM, design)   #fit the model
LRT_3n_Vs_2n_Bleo <- glmTreat(fit_TMM)
summary(decideTests(LRT_3n_Vs_2n_Bleo))

TMM_3n_Bleo <- topTags(LRT_3n_Vs_2n_Bleo, n=nrow(DE_TMM),adjust.method = "BH")$table #assigning fold change gene expression between diploid and triploid ESCs
TMM_3n_Bleo <- data.frame("SYMBOL"=row.names(TMM_3n_Bleo), TMM_3n_Bleo)

colnames(TMM_3n_Bleo)[2] <- "logFC_triploids"

TMM_3n_Bleo <- merge(gene_lengths, TMM_3n_Bleo, by = "SYMBOL", all.y = T)

sig_DE_3n_Bleo <-TMM_3n_Bleo[TMM_3n_Bleo$FDR<=0.05,] #assigning only genes with significant fold change between triploid and diploid ESCs (by FDR) 
sig_DE_3n_Bleo <- sig_DE_3n_Bleo[-rowSums(is.na(sig_DE_3n_Bleo[ ,2:7])) == 0, ] #removing empty rows
sig_DE_3n_Bleo<- sig_DE_3n_Bleo[order(sig_DE_3n_Bleo$FDR, decreasing = F),]
write.csv(sig_DE_3n_Bleo, file = "significant DE genes between Bleo-treated diploid and triploid hESCs.csv")


############# Gene Set Enrichment Analysis Bleo ############# 
library(fgsea)
library(dplyr)

ranks <- TMM_3n_Bleo
ranks[, 'score'] <-sign(ranks$logFC_triploids)* (-log(x = ranks$PValue,base = 10))
.rowNamesDF(ranks, make.names=T) <- ranks$SYMBOL
ranks <- setNames(ranks$score,ranks$SYMBOL)
ranks <- ranks[which(!duplicated(ranks))]


## GSEA of Bleo-treated triploid Vs diploid hESCs:
comb_Terms <- gmtPathways("combined_collection.gmt") 
msigdb_pathways <- comb_Terms
set.seed(42)
fgseaRes <- fgseaMultilevel(msigdb_pathways, ranks, minSize=15, maxSize=500, eps = 0)

fgseaResTidy_Bleo <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy_Bleo %>%
  dplyr::select(-ES) %>%
  arrange(padj) %>%
  DT::datatable()
write.csv(x = data.frame(fgseaResTidy_Bleo[,-8],
                         "leadingEdge"=sapply(fgseaResTidy_Bleo$leadingEdge, function(x) paste(x, collapse = ";"))),
          file = "All GSEA results of Bleo-treated triploid hESCs.csv")

write.csv(x = data.frame(fgseaResTidy_Bleo[fgseaResTidy_Bleo$padj<0.05,-8],
                         "leadingEdge"=sapply(fgseaResTidy_Bleo[fgseaResTidy_Bleo$padj<0.05,8], function(x) paste(x, collapse = ";"))),
          file = "Significant GSEA results of Bleo-treated triploid hESCs.csv")



#### Differential expression (DE) analysis of Cis:

ploidy <- factor(c(1,1,1,2,2,2), labels = c("2n","3n")) #4 diploid, and 4 triploid samples
DE_object <- DGEList(gene_exp_all[,c(10,14,18,22,26,30)], group = ploidy, genes = gene_exp_all$gene)
keep <- filterByExpr(DE_object, min.count=5, group = ploidy)
DE_object <- DE_object[keep, , keep.lib.sizes=FALSE]
DE_TMM <- calcNormFactors(DE_object)

###PCA analysis of log2(CPM) values of each sample
cpm_normalized_Cis <- cpm(DE_TMM, log = T) #log2(cpm) after TMM normalization
cpm_normalized_Cis <- data.frame("SYMBOL" = rownames(cpm_normalized_Cis), cpm_normalized_Cis)
cpm_normalized_Cis <- merge(gene_lengths, cpm_normalized_Cis, by = "SYMBOL",all.y = T)
.rowNamesDF(cpm_normalized_Cis, make.names=T) <- cpm_normalized_Cis$SYMBOL
indRemoved <- which(apply(cpm_normalized_Cis[,8:13], 1, function(x) all(x < 1)) ) #saving all unexpressed genes in the tpm data.frame
cpm_normalized_Cis <- cpm_normalized_Cis[-indRemoved,]      #getting rid of rows (genes) that all of their columns are zeros
write.csv(cpm_normalized_Cis, file = "cpm_normalized_Cis.csv")


#continue DE analysis
design <- model.matrix(~ploidy)
DE_TMM <- estimateDisp(DE_TMM, design, robust=TRUE)
fit_TMM <- glmFit(DE_TMM, design)   #fit the model
LRT_3n_Vs_2n_Cis <- glmTreat(fit_TMM)
summary(decideTests(LRT_3n_Vs_2n_Cis))

TMM_3n_Cis <- topTags(LRT_3n_Vs_2n_Cis, n=nrow(DE_TMM),adjust.method = "BH")$table #assigning fold change gene expression between diploid and triploid ESCs
TMM_3n_Cis <- data.frame("SYMBOL"=row.names(TMM_3n_Cis), TMM_3n_Cis)

colnames(TMM_3n_Cis)[2] <- "logFC_triploids"

TMM_3n_Cis <- merge(gene_lengths, TMM_3n_Cis, by = "SYMBOL", all.y = T)

sig_DE_3n_Cis <-TMM_3n_Cis[TMM_3n_Cis$FDR<=0.05,] #assigning only genes with significant fold change between triploid and diploid ESCs (by FDR) 
sig_DE_3n_Cis <- sig_DE_3n_Cis[-rowSums(is.na(sig_DE_3n_Cis[ ,2:7])) == 0, ] #removing empty rows
sig_DE_3n_Cis<- sig_DE_3n_Cis[order(sig_DE_3n_Cis$logFC_triploids, decreasing = T),]
write.csv(sig_DE_3n_Cis, file = "significant DE genes between Cis-treated diploid and triploid hESCs.csv")


############# Gene Set Enrichment Analysis Cis ############# 
library(fgsea)
library(dplyr)

ranks <- TMM_3n_Cis
ranks[, 'score'] <-sign(ranks$logFC_triploids)* (-log(x = ranks$PValue,base = 10))
.rowNamesDF(ranks, make.names=T) <- ranks$SYMBOL
ranks <- setNames(ranks$score,ranks$SYMBOL)
ranks <- ranks[which(!duplicated(ranks))]


## GSEA of Cis-treated triploid Vs diploid hESCs:
comb_Terms <- gmtPathways("combined_collection.gmt") 
msigdb_pathways <- comb_Terms
set.seed(42)
fgseaRes <- fgseaMultilevel(msigdb_pathways, ranks, minSize=15, maxSize=500, eps = 0)

fgseaResTidy_Cis <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy_Cis %>%
  dplyr::select(-ES) %>%
  arrange(padj) %>%
  DT::datatable()
write.csv(x = data.frame(fgseaResTidy_Cis[,-8],
                         "leadingEdge"=sapply(fgseaResTidy_Cis$leadingEdge, function(x) paste(x, collapse = ";"))),
          file = "All GSEA results of Cis-treated triploid hESCs.csv")

write.csv(x = data.frame(fgseaResTidy_Cis[fgseaResTidy_Cis$padj<0.05,-8],
                         "leadingEdge"=sapply(fgseaResTidy_Cis[fgseaResTidy_Cis$padj<0.05,8], function(x) paste(x, collapse = ";"))),
          file = "Significant GSEA results of Cis-treated triploid hESCs.csv")



#### Differential expression (DE) analysis of Pacli:

ploidy <- factor(c(1,1,1,2,2,2), labels = c("2n","3n")) #4 diploid, and 4 triploid samples
DE_object <- DGEList(gene_exp_all[,c(11,15,19,23,27,31)], group = ploidy, genes = gene_exp_all$gene)
keep <- filterByExpr(DE_object, min.count=5, group = ploidy)
DE_object <- DE_object[keep, , keep.lib.sizes=FALSE]
DE_TMM <- calcNormFactors(DE_object)

###PCA analysis of log2(CPM) values of each sample
cpm_normalized_Pacli <- cpm(DE_TMM, log = T) #log2(cpm) after TMM normalization
cpm_normalized_Pacli <- data.frame("SYMBOL" = rownames(cpm_normalized_Pacli), cpm_normalized_Pacli)
cpm_normalized_Pacli <- merge(gene_lengths, cpm_normalized_Pacli, by = "SYMBOL",all.y = T)
.rowNamesDF(cpm_normalized_Pacli, make.names=T) <- cpm_normalized_Pacli$SYMBOL
indRemoved <- which(apply(cpm_normalized_Pacli[,8:13], 1, function(x) all(x < 1)) ) #saving all unexpressed genes in the tpm data.frame
cpm_normalized_Pacli <- cpm_normalized_Pacli[-indRemoved,]      #getting rid of rows (genes) that all of their columns are zeros
write.csv(cpm_normalized_Pacli, file = "cpm_normalized_Pacli.csv")

#continue DE analysis
design <- model.matrix(~ploidy)
DE_TMM <- estimateDisp(DE_TMM, design, robust=TRUE)
fit_TMM <- glmFit(DE_TMM, design)   #fit the model
LRT_3n_Vs_2n_Pacli <- glmTreat(fit_TMM)
summary(decideTests(LRT_3n_Vs_2n_Pacli))

TMM_3n_Pacli <- topTags(LRT_3n_Vs_2n_Pacli, n=nrow(DE_TMM),adjust.method = "BH")$table #assigning fold change gene expression between diploid and triploid ESCs
TMM_3n_Pacli <- data.frame("SYMBOL"=row.names(TMM_3n_Pacli), TMM_3n_Pacli)

colnames(TMM_3n_Pacli)[2] <- "logFC_triploids"

TMM_3n_Pacli <- merge(gene_lengths, TMM_3n_Pacli, by = "SYMBOL", all.y = T)

sig_DE_3n_Pacli <-TMM_3n_Pacli[TMM_3n_Pacli$FDR<=0.05,] #assigning only genes with significant fold change between triploid and diploid ESCs (by FDR) 
sig_DE_3n_Pacli <- sig_DE_3n_Pacli[-rowSums(is.na(sig_DE_3n_Pacli[ ,2:7])) == 0, ] #removing empty rows
sig_DE_3n_Pacli<- sig_DE_3n_Pacli[order(sig_DE_3n_Pacli$FDR, decreasing = F),]
write.csv(sig_DE_3n_Pacli, file = "significant DE genes between Pacli-treated diploid and triploid hESCs.csv")


############# Gene Set Enrichment Analysis Pacli ############# 
library(fgsea)
library(dplyr)

ranks <- TMM_3n_Pacli
ranks[, 'score'] <-sign(ranks$logFC_triploids)* (-log(x = ranks$PValue,base = 10))
.rowNamesDF(ranks, make.names=T) <- ranks$SYMBOL
ranks <- setNames(ranks$score,ranks$SYMBOL)
ranks <- ranks[which(!duplicated(ranks))]


## GSEA of Pacli-treated triploid Vs diploid hESCs:
comb_Terms <- gmtPathways("combined_collection.gmt") 
msigdb_pathways <- comb_Terms
set.seed(42)
fgseaRes <- fgseaMultilevel(msigdb_pathways, ranks, minSize=15, maxSize=500, eps = 0)

fgseaResTidy_Pacli <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy_Pacli %>%
  dplyr::select(-ES) %>%
  arrange(padj) %>%
  DT::datatable()
write.csv(x = data.frame(fgseaResTidy_Pacli[,-8],
                         "leadingEdge"=sapply(fgseaResTidy_Pacli$leadingEdge, function(x) paste(x, collapse = ";"))),
          file = "All GSEA results of Pacli-treated triploid hESCs.csv")

write.csv(x = data.frame(fgseaResTidy_Pacli[fgseaResTidy_Pacli$padj<0.05,-8],
                         "leadingEdge"=sapply(fgseaResTidy_Pacli[fgseaResTidy_Pacli$padj<0.05,8], function(x) paste(x, collapse = ";"))),
          file = "Significant GSEA results of Pacli-treated triploid hESCs.csv")




##############################################################################
relevant_terms <- c("P53","TP53","G2","G1","_M_","_S_","PHASE","CHECKPOINT","GTSE","DNA","MYC","E2F","MTOR","OXIDATIVE_PHOSPHORYLATION","SENESCENCE",
                    "PATTERN_SPECIFICATION","SKELETAL_SYSTEM_DEVELOPMENT","UV","APOP","DAMAGE","MITOTIC","MITOSIS","SPINDLE","PI3","AKT","P16","CDKN1A",
                    "CDKN2A","CHROMOSOME_SEGREGATION","REPAIR","PCNA","P21","P27","TNFA","NECRO","DOUBLE_STRAND_BREAK","MICROTUBULE","PHASE_TRANSITION")

GSEA_5FU <- read.csv("All GSEA results of 5FU-treated triploid hESCs.csv")
GSEA_Bleo <- read.csv("All GSEA results of Bleo-treated triploid hESCs.csv")
GSEA_Cis <- read.csv("All GSEA results of Cis-treated triploid hESCs.csv")
GSEA_Pacli <- read.csv("All GSEA results of Pacli-treated triploid hESCs.csv")


combined_GSEA <- Reduce(function(x, y) merge(x, y, all=T, by = "pathway"),
                        list(GSEA_5FU[,c(2,4,7)],
                             GSEA_Bleo[,c(2,4,7)],
                             GSEA_Cis[,c(2,4,7)],
                             GSEA_Pacli[,c(2,4,7)]))
colnames(combined_GSEA) <- c("pathway","5FU_p.adj","5FU_NES","Bleo_p.adj","Bleo_NES","Cis_p.adj","Cis_NES","Pacli_p.adj","Pacli_NES")

pathway_matches <- sapply(relevant_terms, function(term) {
  grep(term, combined_GSEA$pathway, ignore.case = TRUE)
})
pathway_matches <- unlist(pathway_matches)
# Remove duplicates (ensures each pathway appears only once)
pathway_matches <- pathway_matches[!duplicated(pathway_matches)]



##### generate a unified GSEA terms list and display the relevant pathways in 3n vs 2n TP53-KO cells:
# Function to create summary of GSEA pathways by treatment
create_gsea_pathway_summary <- function(combined_gsea_df, pathway_groups) {
  
  # Initialize results dataframe
  summary_results <- data.frame()
  
  # Get unique treatments from column names (assuming format like "Treatment_5FU", "Treatment_Bleo", etc.)
  treatment_cols <- colnames(combined_gsea_df)
  treatments <- c("5FU", "Bleo", "Cis", "Pacli")  # Define your treatments
  
  # Process each pathway group
  for (group_name in names(pathway_groups)) {
    
    cat("Processing pathway group:", group_name, "\n")
    
    # Get pathways for this group
    group_pathways <- pathway_groups[[group_name]]
    
    # Find rows that match any of the pathways in this group
    matching_rows <- combined_gsea_df[grepl(paste(group_pathways, collapse = "|"), 
                                            combined_gsea_df$pathway, ignore.case = TRUE), ]
    
    if (nrow(matching_rows) == 0) {
      cat("  No pathways found for group:", group_name, "\n")
      next
    }
    
    cat("  Found", nrow(matching_rows), "pathways for group:", group_name, "\n")
    
    # Process each treatment
    for (treatment in treatments) {
      
      # Find columns for this treatment (padj and NES)
      padj_col <- paste0(treatment, "_p.adj")  # Adjust column naming as needed
      nes_col <- paste0(treatment, "_NES")    # Adjust column naming as needed
      
      # Alternative column naming patterns (adjust based on your actual column names)
      if (!padj_col %in% colnames(matching_rows)) {
        # Try alternative naming patterns
        potential_padj_cols <- grep(paste0(treatment, ".*padj|padj.*", treatment), 
                                    colnames(matching_rows), ignore.case = TRUE, value = TRUE)
        potential_nes_cols <- grep(paste0(treatment, ".*NES|NES.*", treatment), 
                                   colnames(matching_rows), ignore.case = TRUE, value = TRUE)
        
        if (length(potential_padj_cols) > 0) padj_col <- potential_padj_cols[1]
        if (length(potential_nes_cols) > 0) nes_col <- potential_nes_cols[1]
      }
      
      # Check if columns exist
      if (!padj_col %in% colnames(matching_rows) || !nes_col %in% colnames(matching_rows)) {
        cat("  Warning: Columns not found for treatment", treatment, "\n")
        cat("  Looking for:", padj_col, "and", nes_col, "\n")
        cat("  Available columns:", paste(colnames(matching_rows), collapse = ", "), "\n")
        next
      }
      
      # Get valid (non-NA) padj values for this treatment
      valid_rows <- !is.na(matching_rows[[padj_col]])
      
      if (sum(valid_rows) == 0) {
        cat("  No valid padj values for treatment", treatment, "in group", group_name, "\n")
        next
      }
      
      # Find row with minimum padj value
      min_padj_idx <- which.min(matching_rows[[padj_col]][valid_rows])
      # Get the actual row index in the original dataframe
      actual_idx <- which(valid_rows)[min_padj_idx]
      best_row <- matching_rows[actual_idx, ]
      
      # Extract the best values
      best_padj <- best_row[[padj_col]]
      best_nes <- best_row[[nes_col]]
      best_pathway <- best_row$pathway
      
      # Create summary row
      summary_row <- data.frame(
        General_Pathway = group_name,
        Treatment = treatment,
        Best_Pathway = best_pathway,
        Min_padj = best_padj,
        Corresponding_NES = best_nes,
        stringsAsFactors = FALSE
      )
      
      # Add to results
      summary_results <- rbind(summary_results, summary_row)
      
      cat("  ", treatment, "- Best pathway:", substr(best_pathway, 1, 50), "..., padj:", 
          round(best_padj, 6), ", NES:", round(best_nes, 3), "\n")
    }
  }
  
  return(summary_results)
}

# Define pathway groups with common themes
pathway_groups <- list(
  "Cell_Cycle_G1/S_G2/M_Phase_Transition" = c("REACTOME_MITOTIC_G1_PHASE_AND_G1_S_TRANSITION", "REACTOME_CYCLIN_A_B1_B2_ASSOCIATED_EVENTS_DURING_G2_M_TRANSITION",
                                              "REACTOME_M_PHASE","GOBP_CELL_CYCLE_G1_S_PHASE_TRANSITION","REACTOME_CYCLIN_A_CDK2_ASSOCIATED_EVENTS_AT_S_PHASE_ENTRY",
                                              "REACTOME_G1_S_SPECIFIC_TRANSCRIPTION", "HALLMARK_G2M_CHECKPOINT", "GOBP_REGULATION_OF_MITOTIC_CELL_CYCLE","GOBP_CELL_CYCLE_PHASE_TRANSITION",
                                              "GOBP_REGULATION_OF_MITOTIC_CELL_CYCLE_PHASE_TRANSITION","GOBP_CELL_CYCLE_G2_M_PHASE_TRANSITION","REACTOME_G0_AND_EARLY_G1",
                                              "GOBP_REGULATION_OF_CELL_CYCLE_PHASE_TRANSITION","GOBP_REGULATION_OF_CELL_CYCLE_G2_M_PHASE_TRANSITION","REACTOME_S_PHASE",
                                              "GOBP_REGULATION_OF_CELL_CYCLE_G1_S_PHASE_TRANSITION","GOBP_POSITIVE_REGULATION_OF_MITOTIC_CELL_CYCLE_PHASE_TRANSITION",
                                              "GOBP_POSITIVE_REGULATION_OF_MITOTIC_CELL_CYCLE","GOBP_POSITIVE_REGULATION_OF_CELL_CYCLE_PHASE_TRANSITION",
                                              "GOBP_NEGATIVE_REGULATION_OF_MITOTIC_CELL_CYCLE","GOBP_MITOTIC_CELL_CYCLE_PHASE_TRANSITION","REACTOME_G2_M_CHECKPOINTS",
                                              "REACTOME_MITOTIC_PROPHASE","REACTOME_MITOTIC_PROMETAPHASE","REACTOME_CONDENSATION_OF_PROPHASE_CHROMOSOMES",
                                              "REACTOME_APC_C_CDH1_MEDIATED_DEGRADATION_OF_CDC20_AND_OTHER_APC_C_CDH1_TARGETED_PROTEINS_IN_LATE_MITOSIS_EARLY_G1",
                                              "REACTOME_ACTIVATION_OF_APC_C_AND_APC_C_CDC20_MEDIATED_DEGRADATION_OF_MITOTIC_PROTEINS","REACTOME_G2_M_DNA_DAMAGE_CHECKPOINT",
                                              "REACTOME_THE_ROLE_OF_GTSE1_IN_G2_M_PROGRESSION_AFTER_G2_CHECKPOINT","REACTOME_MITOTIC_G2_G2_M_PHASES"),
  "Apoptosis/TNFA_signaling" = c("GOBP_APOPTOTIC_PROCESS_INVOLVED_IN_DEVELOPMENT","REACTOME_REGULATION_OF_APOPTOSIS","HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
  "TP53_Signaling" = c("HALLMARK_P53_PATHWAY","REACTOME_REGULATION_OF_TP53_ACTIVITY_THROUGH_PHOSPHORYLATION"),
  "DNA_Repair" = c("WP_DNA_REPAIR_PATHWAYS_FULL_NETWORK", "WP_DNA_IRDAMAGE_AND_CELLULAR_RESPONSE_VIA_ATR", "REACTOME_HOMOLOGY_DIRECTED_REPAIR","GOBP_RECOMBINATIONAL_REPAIR",
                   "REACTOME_DNA_DOUBLE_STRAND_BREAK_RESPONSE", "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR","KEGG_MISMATCH_REPAIR","GOBP_DNA_INTEGRITY_CHECKPOINT_SIGNALING",
                   "KEGG_MEDICUS_REFERENCE_MISMATCH_REPAIR","HALLMARK_DNA_REPAIR","GOBP_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING",
                   "GOBP_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR","GOBP_REGULATION_OF_DNA_REPAIR","GOBP_POSITIVE_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR","GOBP_POSTREPLICATION_REPAIR",
                   "GOBP_POSITIVE_REGULATION_OF_DNA_REPAIR","GOBP_NUCLEOTIDE_EXCISION_REPAIR","GOBP_INTERSTRAND_CROSS_LINK_REPAIR","GOBP_DOUBLE_STRAND_BREAK_REPAIR","GOBP_MISMATCH_REPAIR",
                   "WP_NUCLEOTIDE_EXCISION_REPAIR","WP_DNA_MISMATCH_REPAIR","WP_BASE_EXCISION_REPAIR","REACTOME_TRANSLESION_SYNTHESIS_BY_Y_FAMILY_DNA_POLYMERASES_BYPASSES_LESIONS_ON_DNA_TEMPLATE",
                   "REACTOME_TRANSCRIPTION_COUPLED_NUCLEOTIDE_EXCISION_REPAIR_TC_NER","REACTOME_TERMINATION_OF_TRANSLESION_DNA_SYNTHESIS","REACTOME_NUCLEOTIDE_EXCISION_REPAIR",
                   "REACTOME_RECOGNITION_OF_DNA_DAMAGE_BY_PCNA_CONTAINING_REPLICATION_COMPLEX","REACTOME_PROCESSING_OF_DNA_DOUBLE_STRAND_BREAK_ENDS","KEGG_NUCLEOTIDE_EXCISION_REPAIR",
                   "REACTOME_PCNA_DEPENDENT_LONG_PATCH_BASE_EXCISION_REPAIR","REACTOME_MISMATCH_REPAIR","REACTOME_GLOBAL_GENOME_NUCLEOTIDE_EXCISION_REPAIR_GG_NER",
                   "REACTOME_GAP_FILLING_DNA_REPAIR_SYNTHESIS_AND_LIGATION_IN_GG_NER","REACTOME_DNA_REPAIR","REACTOME_BASE_EXCISION_REPAIR","KEGG_BASE_EXCISION_REPAIR",
                   "GOBP_DNA_DOUBLE_STRAND_BREAK_PROCESSING","GOBP_BASE_EXCISION_REPAIR","GOBP_REGULATION_OF_DNA_DAMAGE_CHECKPOINT","WP_NUCLEOTIDE_EXCISION_REPAIR_IN_XERODERMA_PIGMENTOSUM",
                   "GOMF_DAMAGED_DNA_BINDING","GOCC_SITE_OF_DOUBLE_STRAND_BREAK","GOCC_SITE_OF_DNA_DAMAGE","GOBP_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR_VIA_HOMOLOGOUS_RECOMBINATION",
                   "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING","GOBP_DNA_TEMPLATED_DNA_REPLICATION_MAINTENANCE_OF_FIDELITY","GOBP_DNA_SYNTHESIS_INVOLVED_IN_DNA_REPAIR",
                   "GOBP_NEGATIVE_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR_VIA_HOMOLOGOUS_RECOMBINATION","REACTOME_DNA_DAMAGE_BYPASS"),
  "DNA_Replication" = c("REACTOME_DNA_REPLICATION", "GOBP_DNA_TEMPLATED_DNA_REPLICATION", "GOBP_DNA_REPLICATION","WP_DNA_REPLICATION","REACTOME_SYNTHESIS_OF_DNA",
                        "REACTOME_SUMOYLATION_OF_DNA_REPLICATION_PROTEINS","REACTOME_DNA_STRAND_ELONGATION","REACTOME_DNA_REPLICATION_PRE_INITIATION","KEGG_DNA_REPLICATION",
                        "GOMF_DNA_REPLICATION_ORIGIN_BINDING","GOMF_DNA_POLYMERASE_BINDING","GOMF_DNA_POLYMERASE_ACTIVITY","GOBP_REGULATION_OF_DNA_TEMPLATED_DNA_REPLICATION",
                        "GOBP_REGULATION_OF_DNA_REPLICATION","GOBP_REGULATION_OF_DNA_METABOLIC_PROCESS","GOBP_POSITIVE_REGULATION_OF_DNA_METABOLIC_PROCESS","GOBP_MITOTIC_DNA_REPLICATION",
                        "GOBP_DNA_STRAND_ELONGATION","GOBP_DNA_REPLICATION_INITIATION","GOBP_DNA_REPLICATION_CHECKPOINT_SIGNALING","GOBP_DNA_STRAND_ELONGATION_INVOLVED_IN_DNA_REPLICATION",
                        "GOBP_DNA_BIOSYNTHETIC_PROCESS","GOBP_CELL_CYCLE_DNA_REPLICATION"),
  "Senescence" = c("WP_SENESCENCEASSOCIATED_SECRETORY_PHENOTYPE_SASP", "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP", "REACTOME_OXIDATIVE_STRESS_INDUCED_SENESCENCE", 
                   "REACTOME_DNA_DAMAGE_TELOMERE_STRESS_INDUCED_SENESCENCE","REACTOME_CELLULAR_SENESCENCE"),
  "PI3K/AKT_Pathway" = c("WP_PI3KAKT_SIGNALING","WP_FOCAL_ADHESION_PI3KAKTMTORSIGNALING"),
  "Pattern_Specification_Process" = c("GOBP_PATTERN_SPECIFICATION_PROCESS"),
  "E2F_Pathway" = c("REACTOME_TRANSCRIPTION_OF_E2F_TARGETS_UNDER_NEGATIVE_CONTROL_BY_DREAM_COMPLEX","PID_E2F_PATHWAY","HALLMARK_E2F_TARGETS"),
  "MYC_Pathway" = c("PID_MYC_ACTIV_PATHWAY","HALLMARK_MYC_TARGETS_V2","HALLMARK_MYC_TARGETS_V1"),
  "Microtubule_Organization/Spindle_Function" = c("REACTOME_MITOTIC_SPINDLE_CHECKPOINT","GOBP_REGULATION_OF_MICROTUBULE_CYTOSKELETON_ORGANIZATION",
                                                  "REACTOME_INHIBITION_OF_THE_PROTEOLYTIC_ACTIVITY_OF_APC_C_REQUIRED_FOR_THE_ONSET_OF_ANAPHASE_BY_MITOTIC_SPINDLE_CHECKPOINT_COMPONENTS",
                                                  "GOBP_MICROTUBULE_POLYMERIZATION_OR_DEPOLYMERIZATION", "GOCC_MICROTUBULE","KEGG_MEDICUS_REFERENCE_SPINDLE_ASSEMBLY_CHECKPOINT_SIGNALING",
                                                  "KEGG_MEDICUS_REFERENCE_PROMOTION_OF_MICROTUBULE_GROWTH","KEGG_MEDICUS_REFERENCE_MICROTUBULE_NUCLEATION",
                                                  "KEGG_MEDICUS_REFERENCE_MICROTUBULE_DEPOLYMERIZATION_AT_THE_MINUS_ENDS","GOBP_MICROTUBULE_DEPOLYMERIZATION","GOCC_SPINDLE",
                                                  "WP_MICROTUBULE_CYTOSKELETON_REGULATION","GOBP_REGULATION_OF_MITOTIC_SISTER_CHROMATID_SEGREGATION","GOBP_MITOTIC_NUCLEAR_DIVISION",
                                                  "GOBP_REGULATION_OF_CHROMOSOME_SEGREGATION","GOBP_POSITIVE_REGULATION_OF_MITOTIC_SISTER_CHROMATID_SEPARATION","GOBP_CHROMOSOME_SEGREGATION",
                                                  "GOBP_MITOTIC_METAPHASE_CHROMOSOME_ALIGNMENT","GOBP_METAPHASE_CHROMOSOME_ALIGNMENT","GOBP_MEIOTIC_CHROMOSOME_SEGREGATION",
                                                  "REACTOME_MITOTIC_PROMETAPHASE","REACTOME_MITOTIC_METAPHASE_AND_ANAPHASE",
                                                  "GOBP_ATTACHMENT_OF_SPINDLE_MICROTUBULES_TO_KINETOCHORE","GOBP_NUCLEAR_CHROMOSOME_SEGREGATION","REACTOME_HOMOLOGOUS_DNA_PAIRING_AND_STRAND_EXCHANGE",
                                                  "GOBP_MITOTIC_SISTER_CHROMATID_SEPARATION","GOBP_MITOTIC_SISTER_CHROMATID_SEGREGATION","GOBP_MITOTIC_SISTER_CHROMATID_COHESION"),
  "Oxidative Phosphorylation" = c("WP_OXIDATIVE_PHOSPHORYLATION","KEGG_OXIDATIVE_PHOSPHORYLATION","HALLMARK_OXIDATIVE_PHOSPHORYLATION","GOBP_OXIDATIVE_PHOSPHORYLATION")
)

# Example usage:
# Assuming your combined_gsea_df has columns like:
# pathway, 5FU_padj, 5FU_NES, Bleo_padj, Bleo_NES, Cis_padj, Cis_NES, Pacli_padj, Pacli_NES

# Create the summary
gsea_summary <- create_gsea_pathway_summary(combined_GSEA, pathway_groups)

# View results
ggplot(data = gsea_summary, aes(General_Pathway, -log2(Min_padj)*sign(Corresponding_NES),group = Treatment)) +
  geom_col(aes(fill=Treatment),width = 0.5,position = position_dodge(0.8)) +
  coord_flip() +
    scale_fill_manual(values = c("goldenrod1", "orange", "tan3", "chocolate3")) +
labs(x="Pathway", y=expression(-log[2](FDR))) + 
  scale_x_discrete(limits = c("DNA_Replication",
                              "Cell_Cycle_G1/S_G2/M_Phase_Transition",
                              "Microtubule_Organization/Spindle_Function",
                              "DNA_Repair",
                              "Apoptosis/TNFA_signaling",
                              "Pattern_Specification_Process")) +
  theme_minimal()+ ggtitle("WT 3n/2n Unified GSEA")+
  geom_hline(yintercept = -log2(0.05), linetype=2,color = "grey", size=0.25)+
  geom_hline(yintercept = log2(0.05), linetype=2,color = "grey", size=0.25)+
  geom_hline(yintercept = -log2(1), linetype=1,color = "black", size=0.25)+
  theme(axis.text.y = element_text(size = 10))+scale_y_continuous(limits = c(-50,120))
# ggsave(filename="GSEA_relevant_general_terms_in_3n_WT_bars.pdf")