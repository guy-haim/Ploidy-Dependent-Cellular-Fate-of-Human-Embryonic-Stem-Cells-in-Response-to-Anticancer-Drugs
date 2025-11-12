library(stringr)
library(dplyr)
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
library(hrbrthemes)
library(rstatix)
library(viridis)

setwd(choose.dir())
## In this analysis, we assumed that near-3n CCLs will be more resistant to chemotherapis
## than Near-2n CCLs, so we test whether the near-3n CCLs' IC50 values are higher than 
## the near-2n CCLs' IC50 values. 
## To do that, when we use the function t_test() the alternative = "less", 
## and when use the function t.test() the alternative = "greater".
## We also assume p53-deficient CCLs are more resistant than p53-proficient CCLs.


GDSC_data <- read.csv("TableS4A - Whole set of natural log(IC50s) across all the screened compounds and cell lines, related to Figure 4.csv",skip = 4)
colnames(GDSC_data)[1] <- "COSMIC_ID"

CCLs_data <- read.delim("CellLinesProject_Sample_v99_GRCh38.tsv")
CCLs_data$COSMIC_SAMPLE_ID <- as.numeric(str_remove(CCLs_data$COSMIC_SAMPLE_ID,"COSS"))
colnames(CCLs_data)[1] <- "COSMIC_ID"
CCLs_data <- na.omit(CCLs_data[,c(1,2,12,19,20)])
CCLs_data <- CCLs_data[CCLs_data$AVERAGE_PLOIDY<=3.5,]

CCLs_data$Ploidy_group <- CCLs_data$AVERAGE_PLOIDY
CCLs_data$Ploidy_group[CCLs_data$AVERAGE_PLOIDY<=2.5] <- "Near-2n"
CCLs_data$Ploidy_group[CCLs_data$AVERAGE_PLOIDY>2.5 & CCLs_data$AVERAGE_PLOIDY<=3.5] <- "Near-3n"

CCLs_data <- merge(CCLs_data,GDSC_data, by= "COSMIC_ID")
cl <- c("Near-2n"="darkblue","Near-3n"="darkred")

models <- read.csv("model_list_latest.csv")
models <- models[models$model_type=="Cell Line",]
models <- models[,c(2,12,23:26)]
names(models)[1]<- "SAMPLE_NAME"
CCLs_data <- merge(CCLs_data,models[,c(1,3,5,6)],by="SAMPLE_NAME")
CCLs_data <- CCLs_data[,c(1:7,273:275,8:272)]

mutations_data <- read.csv("Damaging_Mutations.csv")
P53_mutations_data <- mutations_data[,c(1:5,which(colnames(mutations_data)=="TP53"))]
P53_mutations_data$TP53_DM <- P53_mutations_data$TP53>0
colnames(P53_mutations_data)[2:3] <- c("SAMPLE_NAME","tissue")

###make sample name in uupper case without dashes or dots or whatever (only number and letters)##
CCLs_data$SAMPLE_NAME <- str_to_upper(str_remove_all(CCLs_data$SAMPLE_NAME,"-"))
P53_mutations_data$SAMPLE_NAME <- str_to_upper(str_remove_all(P53_mutations_data$SAMPLE_NAME,"-"))
CCLs_data <- merge(CCLs_data,P53_mutations_data,by="SAMPLE_NAME",all.x=T)

sum(CCLs_data$TP53_DM,na.rm = T)
Nutlin_data <- read.csv("Nutlin_3a_GDSC.csv")
colnames(Nutlin_data)[1] <- "SAMPLE_NAME"
Nutlin_data$SAMPLE_NAME <- str_to_upper(str_remove_all(Nutlin_data$SAMPLE_NAME,"-"))
colnames(Nutlin_data)[5:6] <- c("Nutlin_IC50","Nutlin_AUC")
Nutlin_data$Nutlin_resistant <- Nutlin_data$Nutlin_IC50>24 & Nutlin_data$Nutlin_AUC>0.9375
CCLs_data <- merge(CCLs_data,Nutlin_data[,c(1,5:7)], by="SAMPLE_NAME")
CCLs_data <- CCLs_data[,c(1:10,276,278:284,11:275)]


###### Bleomycin ######
Bleo <- CCLs_data[,c(1:18,which(colnames(CCLs_data)=="Bleomycin"))]
Bleo$Bleomycin <- exp(1)^Bleo$Bleomycin
Bleo <- na.omit(Bleo)
Bleo<- Bleo[Bleo$TUMOUR_SOURCE!="NS" & Bleo$TUMOUR_SOURCE!="recurrent",]

stat=compare_means(Bleomycin~Ploidy_group, data = Bleo, method="t.test",alternative = "greater", p.adjust.method = "BH")
ggbarplot(Bleo, x= "Ploidy_group", y= "Bleomycin",fill="Ploidy_group", add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.adj", y.position = c(95), size = 3,tip.length = 0.001, hide.ns = F)+
  theme_classic()+theme(legend.position = "none")+scale_y_continuous(expand = c(0,0))+
  xlab("Ploidy group")+ylab("Mean IC50")+coord_cartesian(ylim = c(0,100))+
  scale_fill_manual(values = cl)+ggtitle("Mean IC50 for Bleo by ploidy groups")
# ggsave("Barplot of Mean IC50 for Bleo by ploidy groups.pdf",device = "pdf")



###### Cisplatin ######
Cis <- CCLs_data[,c(1:18,which(colnames(CCLs_data)=="Cisplatin"))]
Cis$Cisplatin <- exp(1)^Cis$Cisplatin
Cis <- na.omit(Cis)
Cis<- Cis[Cis$TUMOUR_SOURCE!="NS" & Cis$TUMOUR_SOURCE!="recurrent",]

stat=compare_means(Cisplatin~Ploidy_group, data = Cis, method="t.test", alternative="greater", ref.group = "Near-2n", p.adjust.method = "BH")
ggbarplot(Cis, x= "Ploidy_group", y= "Cisplatin",fill="Ploidy_group", add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.adj", y.position = c(55), size = 3,tip.length = 0.001, hide.ns = F)+
  theme_classic()+theme(legend.position = "none")+scale_y_continuous(expand = c(0,0))+
  xlab("Ploidy group")+ylab("Mean IC50")+coord_cartesian(ylim = c(0,60))+
  scale_fill_manual(values = cl)+ggtitle("Mean IC50 for Cis by ploidy groups")
# ggsave("Barplot of Mean IC50 for Cis by ploidy groups.pdf",device = "pdf")



###### Paclitaxel ######
Pacli <- CCLs_data[,c(1:18,which(colnames(CCLs_data)=="Paclitaxel"))]
Pacli$Paclitaxel <- exp(1)^Pacli$Paclitaxel
Pacli <- na.omit(Pacli)
Pacli<- Pacli[Pacli$TUMOUR_SOURCE!="NS" & Pacli$TUMOUR_SOURCE!="recurrent",]

stat=compare_means(Paclitaxel~Ploidy_group, data = Pacli, method="t.test", alternative="greater",ref.group = "Near-2n", p.adjust.method = "BH")
ggbarplot(Pacli, x= "Ploidy_group", y= "Paclitaxel",fill="Ploidy_group", add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.adj", y.position = c(0.45), size = 3,tip.length = 0.001, hide.ns = F)+
  theme_classic()+theme(legend.position = "none")+scale_y_continuous(expand = c(0,0))+
  xlab("Ploidy group")+ylab("Mean IC50")+coord_cartesian(ylim = c(0,0.5))+
  scale_fill_manual(values = cl)+ggtitle("Mean IC50 for Pacli by ploidy groups")
# ggsave("Barplot of Mean IC50 for Pacli by ploidy groups.pdf",device = "pdf")



###### X5.Fluorouracil ######
FU <- CCLs_data[,c(1:18,which(colnames(CCLs_data)=="X5.Fluorouracil"))]
FU$X5.Fluorouracil <- exp(1)^FU$X5.Fluorouracil
FU <- na.omit(FU)
FU<- FU[FU$TUMOUR_SOURCE!="NS" & FU$TUMOUR_SOURCE!="recurrent",]

stat=compare_means(X5.Fluorouracil~Ploidy_group, data = FU, method="t.test",  alternative="greater",ref.group = "Near-2n", p.adjust.method = "BH")
ggbarplot(FU, x= "Ploidy_group", y= "X5.Fluorouracil",fill="Ploidy_group", add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.adj", y.position = c(185), size = 3,tip.length = 0.001, hide.ns = F)+
  theme_classic()+theme(legend.position = "none")+scale_y_continuous(expand = c(0,0))+
  xlab("Ploidy group")+ylab("Mean IC50")+coord_cartesian(ylim = c(0,200))+
  scale_fill_manual(values = cl)+ggtitle("Mean IC50 for FU by ploidy groups")
# ggsave("Barplot of Mean IC50 for FU by ploidy groups.pdf",device = "pdf")