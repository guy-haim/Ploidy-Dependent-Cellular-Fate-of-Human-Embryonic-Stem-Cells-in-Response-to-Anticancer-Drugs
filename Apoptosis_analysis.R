library(ggplot2)
library(ggpubr)
library(tidyr)
library(rstatix)

setwd(choose.dir())
DsRed_2n_rep1 <- read.csv("2n DsRed rep1.csv")
DsRed_2n_rep1$Sample_name <- rep("DsRed_2n_rep1",nrow(DsRed_2n_rep1))
DsRed_2n_rep1$Ploidy <- rep("2n",nrow(DsRed_2n_rep1))
DsRed_2n_rep1 <- DsRed_2n_rep1[DsRed_2n_rep1$Depth=="> > ",c(5,6,2,3)]

DsRed_2n_rep2 <- read.csv("2n DsRed rep2.csv")
DsRed_2n_rep2$Sample_name <- rep("DsRed_2n_rep2",nrow(DsRed_2n_rep2))
DsRed_2n_rep2$Ploidy <- rep("2n",nrow(DsRed_2n_rep2))
DsRed_2n_rep2 <- DsRed_2n_rep2[DsRed_2n_rep2$Depth=="> > ",c(5,6,2,3)]

EGFP_2n <- read.csv("2n EGFP.csv")
EGFP_2n$Sample_name <- rep("EGFP_2n",nrow(EGFP_2n))
EGFP_2n$Ploidy <- rep("2n",nrow(EGFP_2n))
EGFP_2n <- EGFP_2n[EGFP_2n$Depth=="> > ",c(5,6,2,3)]

Clone_H_3n_rep1 <- read.csv("3n Clone H rep1.csv")
Clone_H_3n_rep1$Sample_name <- rep("Clone_H_3n_rep1",nrow(Clone_H_3n_rep1))
Clone_H_3n_rep1$Ploidy <- rep("3n",nrow(Clone_H_3n_rep1))
Clone_H_3n_rep1 <- Clone_H_3n_rep1[Clone_H_3n_rep1$Depth=="> > ",c(5,6,2,3)]

Clone_H_3n_rep2 <- read.csv("3n Clone H rep2.csv")
Clone_H_3n_rep2$Sample_name <- rep("Clone_H_3n_rep2",nrow(Clone_H_3n_rep2))
Clone_H_3n_rep2$Ploidy <- rep("3n",nrow(Clone_H_3n_rep2))
Clone_H_3n_rep2 <- Clone_H_3n_rep2[Clone_H_3n_rep2$Depth=="> > ",c(5,6,2,3)]

Clone_K_3n <- read.csv("3n Clone K.csv")
Clone_K_3n$Sample_name <- rep("Clone_K_3n",nrow(Clone_K_3n))
Clone_K_3n$Ploidy <- rep("3n",nrow(Clone_K_3n))
Clone_K_3n <- Clone_K_3n[Clone_K_3n$Depth=="> > ",c(5,6,2,3)]

Clone_L_3n <- read.csv("3n Clone L.csv")
Clone_L_3n$Sample_name <- rep("Clone_L_3n",nrow(Clone_L_3n))
Clone_L_3n$Ploidy <- rep("3n",nrow(Clone_L_3n))
Clone_L_3n <- Clone_L_3n[Clone_L_3n$Depth=="> > ",c(5,6,2,3)]

Apoptosis_results <- rbind(DsRed_2n_rep1,DsRed_2n_rep2,EGFP_2n,Clone_H_3n_rep1,Clone_H_3n_rep2,Clone_K_3n,Clone_L_3n)
Apoptosis_results$Cell_state <- rep(0,nrow(Apoptosis_results))
Apoptosis_results$Cell_state[grep(pattern = "Q1",Apoptosis_results$Name)] <- "PI.positive"
Apoptosis_results$Cell_state[grep(pattern = "Q2",Apoptosis_results$Name)] <- "PI+Annexin.positive"
Apoptosis_results$Cell_state[grep(pattern = "Q3",Apoptosis_results$Name)] <- "Annexin.positive"
Apoptosis_results$Cell_state[grep(pattern = "Q4",Apoptosis_results$Name)] <- "Negative"
Apoptosis_results$Name[grep(pattern = "Unstained",Apoptosis_results$Name)] <- "Unstained"
Apoptosis_results$Name[grep(pattern = "UT",Apoptosis_results$Name)] <- "Untreated"
Apoptosis_results$Name[grep(pattern = "Bleo",Apoptosis_results$Name)] <- "Bleo"
Apoptosis_results$Name[grep(pattern = "Cis",Apoptosis_results$Name)] <- "Cis"
Apoptosis_results$Name[grep(pattern = "Pacli",Apoptosis_results$Name)] <- "Pacli"
Apoptosis_results$Name[grep(pattern = "FU",Apoptosis_results$Name)] <- "5-FU"
colnames(Apoptosis_results)[3:4] <- c("Treatment","Cell_percentage")
Apoptosis_results$Cell_percentage <- as.numeric(format(Apoptosis_results$Cell_percentage))
Apoptosis_results <- Apoptosis_results[Apoptosis_results$Treatment %in% c("Untreated","Bleo","Pacli","5-FU","Cis"),]


#### Annexin-positive only analysis:
Annexin_only_positive_results <- Apoptosis_results[Apoptosis_results$Cell_state=="PI+Annexin.positive",]
Double_positive_results <- Apoptosis_results[Apoptosis_results$Cell_state=="Annexin.positive",]

Annexin_positive_results <- cbind(Annexin_only_positive_results,Double_positive_results)
Annexin_positive_results <- data.frame(Annexin_positive_results[,1:3],
                                       "Cell_percentage"=Annexin_positive_results[,4]+Annexin_positive_results[,9])


## Relative apoptotic cell percentage:
for (i in unique(Annexin_positive_results$Treatment)){
  Annexin_positive_results[Annexin_positive_results$Treatment==i,4] <- Annexin_positive_results[Annexin_positive_results$Treatment==i,4]/
    mean(Annexin_positive_results[Annexin_positive_results$Treatment==i & Annexin_positive_results$Ploidy=="2n",4])
}

stat=compare_means(Cell_percentage ~ Ploidy,data = Annexin_positive_results, method = "t.test", alternative= "greater", group.by="Treatment")
ggbarplot(data = Annexin_positive_results, x = "Treatment",y = "Cell_percentage",order = c("Untreated","5-FU","Pacli","Bleo","Cis"),
                position=position_dodge(0.8),fill="Ploidy",add="mean_se") +
        stat_pvalue_manual(data = stat, x = "Treatment", label = "p.format",y.position = c(8),
                            size = 4, tip.length = 0.01, hide.ns = F)+ylab("Relative apoptotic cell percentage")+
        theme_minimal()+theme(legend.position = "top")+scale_fill_manual(values=c("royalblue3","red3"))
# ggsave(filename="Relative apoptotic cells percentage for diploids and triploids.pdf")
