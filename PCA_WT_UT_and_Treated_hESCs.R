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
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(viridis)

setwd(choose.dir()) #chemotherpy experiment/RNA-seq samples post-treatment

#loading RNA-seq data from each sample into a data.frame and getting rid of unnecessary data (first 4 rows and second and third columns)
dsRed_4c_1 <- read.delim("dsRed-hygro4cp40_S3_counts.txt", skip = 1)[,c(1,7)]
dsRed_4c_2 <- read.delim("G1_10R_4c_p38_counts.txt", skip = 1)[,c(1,7)]
EGFP_4c_1 <- read.delim("EGFP-neo-4cp39_S10_counts.txt", skip = 1)[,c(1,7)]
EGFP_4c_2 <- read.delim("G2_10G_4c_p38_counts.txt", skip = 1)[,c(1,7)]
fusionI_3n_clone_H <- read.delim("fusionI-3n-f-4-cloneH_S8_counts.txt", skip = 1)[,c(1,7)]
fusionI_3n_clone_L_1 <- read.delim("fusion-3n-f5-cloneL_counts.txt", skip = 1)[,c(1,7)]
fusion_3n_Clone_L_2 <- read.delim("Clone-L-3n_S9_counts.txt", skip = 1)[,c(1,7)]
fusionI_3n_clone_I <- read.delim("fusionI-3n-f-4-cloneI_S6_counts.txt", skip = 1)[,c(1,7)]
fusion_3n_Clone_K <- read.delim("Clone-K-3n_S8_counts.txt", skip = 1)[,c(1,7)]
fusion_3n_Clone_J <- read.delim("Clone-J-3n_S12_counts.txt", skip = 1)[,c(1,7)]
fusion_3n_Clone_B <- read.delim("Clone-B-3n_S7_counts.txt", skip = 1)[,c(1,7)]

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
                            Triploid_J_5FU_2.5uM,Triploid_J_Bleo_20nM,Triploid_J_Cis_100nM,Triploid_J_Pacli_3.3nM,
                            dsRed_4c_1,dsRed_4c_2,EGFP_4c_1,EGFP_4c_2,fusionI_3n_clone_H,fusionI_3n_clone_I,fusionI_3n_clone_L_1,
                            fusion_3n_Clone_L_2,fusion_3n_Clone_K,fusion_3n_Clone_J,fusion_3n_Clone_B))

colnames(gene_exp_all) <- c("ENSEMBL","Diploid rep1 5FU 2.5uM","Diploid rep1 Bleo 20nM","Diploid rep1 Cis 100nM","Diploid rep1 Pacli 3.3nM",
                            "Diploid rep2 5FU 2.5uM","Diploid rep2 Bleo 20nM","Diploid rep2 Cis 100nM","Diploid rep2 Pacli 3.3nM",
                            "Diploid rep3 5FU 2.5uM","Diploid rep3 Bleo 20nM","Diploid rep3 Cis 100nM","Diploid rep3 Pacli 3.3nM",
                            "Triploid B 5FU 2.5uM","Triploid B Bleo 20nM","Triploid B Cis 100nM","Triploid B Pacli 3.3nM",
                            "Triploid E 5FU 2.5uM","Triploid E Bleo 20nM","Triploid E Cis 100nM","Triploid E Pacli 3.3nM",
                            "Triploid J 5FU 2.5uM","Triploid J Bleo 20nM","Triploid J Cis 100nM","Triploid J Pacli 3.3nM",
                            "Untreated Diploid dsRed rep1","Untreated Diploid dsRed rep2","Untreated Diploid EGFP rep1","Untreated Diploid EGFP rep2",
                            "Untreated Triploid clone H","Untreated Triploid clone I","Untreated Triploid clone L rep1","Untreated Triploid clone L rep2",
                            "Untreated Triploid clone K","Untreated Triploid clone J","Untreated Triploid clone B")


gene_lengths <- read.csv("gencode.v34.transcriptLengths.csv")[,-8]

#removing the end of the ENSEMBL ID so that the same genes from different versions of the annotation file will be joined together  
gene_exp_all$ENSEMBL <- str_remove(gene_exp_all$ENSEMBL, "[.].*")
gene_lengths$ENSEMBL <- str_remove(gene_lengths$ENSEMBL, "[.].*")
gene_exp_all <- aggregate(gene_exp_all[-1], list(gene_exp_all$ENSEMBL), FUN = sum, na.rm = T)
colnames(gene_exp_all)[1] <- "ENSEMBL"

#merging all data frame to one count table with all the data needed:
gene_exp_all <- merge(gene_lengths, gene_exp_all, by = "ENSEMBL", all.y = T)
gene_exp_all <- gene_exp_all[rowSums(is.na(gene_exp_all[ ,2:42])) == 0, ] #removing empty rows
indRemoved <- which(apply(gene_exp_all[,8:42], 1, function(x) all(x == 0)) ) #saving all unexpressed genes in the tpm data.frame
gene_exp_all <- gene_exp_all[-indRemoved,]      #getting rid of rows (genes) that all of their columns are zeros
.rowNamesDF(gene_exp_all, make.names=T) <- gene_exp_all$SYMBOL

cl=c("blue1","blue2","blue3","blue4","deepskyblue1","deepskyblue3","royalblue1","royalblue3","dodgerblue1","dodgerblue2","dodgerblue3","dodgerblue4",
     "red1","red2","red3","red4","orangered1","orangered2","orangered3","salmon3","firebrick1","firebrick2","firebrick3","brown",
     "cadetblue1","cadetblue2","cadetblue3","cadetblue4",
     "coral","coral1","coral2","coral3","indianred1","indianred2","indianred3")
names(cl) <- str_replace_all(colnames(gene_exp_all[,8:42]),"_", " ")



#### Differential expression (DE) analysis:

ploidy <- factor(c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,2,2,2,2,2,2,2), labels = c("2n","3n")) #4 diploid, and 4 triploid samples
treatment <-   factor(c(rep(2:5,6),rep(1,11)), labels = c("Untreated","5FU","Bleo","Cis","Pacli"))
DE_object <- DGEList(gene_exp_all[,8:42], group = ploidy, genes = gene_exp_all$gene)
keep <- filterByExpr(DE_object, min.count=5, group = ploidy)
DE_object <- DE_object[keep, , keep.lib.sizes=FALSE]
DE_TMM <- calcNormFactors(DE_object)

###PCA analysis of log2(CPM) values of each sample
cpm_normalized <- cpm(DE_TMM, log = T) #log2(cpm) after TMM normalization
cpm_normalized <- data.frame("SYMBOL" = rownames(cpm_normalized), cpm_normalized)
cpm_normalized <- merge(gene_lengths, cpm_normalized, by = "SYMBOL",all.y = T)
.rowNamesDF(cpm_normalized, make.names=T) <- cpm_normalized$SYMBOL

prin_comp <- prcomp(t(cpm_normalized[,8:42]), center = T, retx = T)
gg <- cbind.data.frame(prin_comp$x, "Sample" =colnames(gene_exp_all)[8:42])
gg <- data.frame(gg, "Ploidy"=ploidy,"Treatment"=treatment)
autoplot(prin_comp, data = gg,shape="Treatment", colour = cl,scale = T, label =T,size=2)+ scale_colour_manual(values=cl)+theme_bw()+
  guides(shape = guide_legend(order = 1),colour = guide_legend(override.aes = list(shape=c(15))))+theme(text = element_text(size = 16))
# ggsave("PCA_all_untreated_and_treated_samples.pdf",device = "pdf")

summary(prin_comp)

#not labeled PCA plot:
gg <- cbind.data.frame(prin_comp$x, "Sample" =names(cl))
gg <- data.frame(gg, "Ploidy"=ploidy,"Treatment"=treatment)
autoplot(prin_comp, data = gg,shape="Ploidy", colour = "Sample", label =F,size=3)+ scale_colour_manual(values=cl)+theme_bw()+
  guides(shape = guide_legend(order = 1),colour = guide_legend(override.aes = list(shape=c(15))))+theme(text = element_text(size = 16))
# ggsave("PCA_all_untreated_and_treated_samples_(not_labeled).pdf",device = "pdf")


autoplot(prin_comp, data = gg, shape = "Treatment", colour = "Ploidy",
         label = FALSE, size = 5) + scale_colour_manual(values = c("royalblue","red3")) +
  scale_shape_manual(values = c(1,15,16,17,18)) +  # Use 5 different shapes
  theme_classic() +
  guides(shape = guide_legend(nrow = 5, order = 1, title = "Treatment"),
         colour = guide_legend(override.aes = list(shape = 15), ncol = 2, title = "Ploidy")) +
  theme(text = element_text(size = 16))
# ggsave("PCA_all_untreated_and_treated_samples_(not_labeled)_2.pdf",device = "pdf")
