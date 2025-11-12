library(stringr)
library(dplyr)
library(edgeR)
library(ggplot2)
library(ggfortify)
library(gplots)
library(ggsignif)
library(ggpubr)

setwd("C:/Users/guyha/Desktop/chemotherpy experiment/DAPI (chemo treatments)")

FU_2n_10G <- read.csv("Results 2n 10G 5-FU.csv")
Bleo_2n_10G <- read.csv("Results 2n 10G Bleo.csv")
Cis_2n_10G <- read.csv("Results 2n 10G Cis.csv")
Pacli_2n_10G <- read.csv("Results 2n 10G Pacli.csv")
UT_2n_10G <- read.csv("Results 2n 10G UT.csv")
FU_2n_10R <- read.csv("Results 2n 10R 5-FU.csv")
Bleo_2n_10R <- read.csv("Results 2n 10R Bleo.csv")
Cis_2n_10R <- read.csv("Results 2n 10R Cis.csv")
Pacli_2n_10R <- read.csv("Results 2n 10R Pacli.csv")
UT_2n_10R <- read.csv("Results 2n 10R UT.csv")


FU_H <- read.csv("Results 3n clone H 5-FU.csv")
Cis_H <- read.csv("Results 3n clone H Cis.csv")
Pacli_H <- read.csv("Results 3n clone H Pacli.csv")
Bleo_H <- read.csv("Results 3n clone H Bleo.csv")
UT_H <- read.csv("Results 3n clone H UT.csv")
FU_L <- read.csv("Results 3n clone L 5-FU.csv")
Cis_L <- read.csv("Results 3n clone L Cis.csv")
Pacli_L <- read.csv("Results 3n clone L Pacli.csv")
Bleo_L <- read.csv("Results 3n clone L Bleo.csv")
UT_L <- read.csv("Results 3n clone L UT.csv")
FU_K <- read.csv("Results 3n clone K 5-FU.csv")
Cis_K <- read.csv("Results 3n clone K Cis.csv")
Pacli_K <- read.csv("Results 3n clone K Pacli.csv")
Bleo_K <- read.csv("Results 3n clone K Bleo.csv")
UT_K <- read.csv("Results 3n clone K UT.csv")


Chemo_response <- rbind(FU_2n_10G,Bleo_2n_10G,Cis_2n_10G,
                        Pacli_2n_10G,UT_2n_10G,
                        FU_2n_10R,Bleo_2n_10R,Cis_2n_10R,
                        Pacli_2n_10R,UT_2n_10R,
                        FU_H,Cis_H,Pacli_H,Bleo_H,UT_H,
                        FU_L,Cis_L,Pacli_L,Bleo_L,UT_L,
                        FU_K,Cis_K,Pacli_K,Bleo_K,UT_K)
Chemo_response$X <- c(rep("2n 10G",4030),rep("2n 10R",4173),
                      rep("3n clone H",3398),rep("3n clone L",2075),
                      rep("3n clone K",1661))
Chemo_response$Ploidy <- c(rep("2n",8203),rep("3n",7134))
Chemo_response$Treatment <- c(rep("5-Fluorouracil",672),rep("Cisplatin",613),
                              rep("Paclitaxel",772),rep("Bleomycin",730),
                              rep("Untreated",1243),
                              rep("5-Fluorouracil",828),rep("Cisplatin",507),
                              rep("Paclitaxel",1157),rep("Bleomycin",648),
                              rep("Untreated",1033),
                              rep("5-Fluorouracil",681),rep("Cisplatin",393),
                              rep("Paclitaxel",899),rep("Bleomycin",630),
                              rep("Untreated",795),
                              rep("5-Fluorouracil",231),rep("Cisplatin",275),
                              rep("Paclitaxel",435),rep("Bleomycin",548),
                              rep("Untreated",586),
                              rep("5-Fluorouracil",143),rep("Cisplatin",238),
                              rep("Paclitaxel",437),rep("Bleomycin",352),
                              rep("Untreated",491))


Chemo_response$Relative_Area <- Chemo_response$Area/mean(Chemo_response[Chemo_response$Ploidy=="2n" & Chemo_response$Treatment=="Untreated",2])


Chemo_response$Radius <- (Chemo_response$Area/(4*pi))^0.5
Chemo_response$Volume <- ((Chemo_response$Radius)^3)*pi*4/3

Chemo_response$Relative_Volume_to_2n <- Chemo_response$Volume/mean(Chemo_response[Chemo_response$Ploidy=="2n" & Chemo_response$Treatment=="Untreated",12])

Chemo_response$Relative_Volume_by_ploidy <- rep(0,nrow(Chemo_response))
Chemo_response[Chemo_response$Ploidy=="2n",14] <- Chemo_response[Chemo_response$Ploidy=="2n",12]/mean(Chemo_response[Chemo_response$Ploidy=="2n" & Chemo_response$Treatment=="Untreated",12])
Chemo_response[Chemo_response$Ploidy=="3n",14] <- Chemo_response[Chemo_response$Ploidy=="3n",12]/mean(Chemo_response[Chemo_response$Ploidy=="3n" & Chemo_response$Treatment=="Untreated",12])
stat=compare_means(Volume~Treatment,eps=0, data = Chemo_response, method="t.test",group.by = "Ploidy",ref.group = "Untreated")

ggbarplot(Chemo_response, x= "Ploidy", y= "Relative_Volume_by_ploidy", add="mean_se",position = position_dodge(0.8),fill = "Treatment")+
  stat_pvalue_manual(data = stat, label = "p", y.position = c(2.6),vjust = 1.75,x = "Ploidy", size = 3,tip.length = 0.02, hide.ns = F)+
  theme_classic()+theme(legend.position = "top")+
  xlab("Ploidy")+ylab("Relative Nuclei Volume")+geom_hline(yintercept = 1.25,linetype = 2,colour = "darkgrey")+
  scale_y_continuous(expand = c(0,0),n.breaks = 10)+scale_fill_manual(values = c("#FE9929","#FEC44F","#FEE391","#FFF7BC","lightgrey"))
# ggsave("Barplot of mean nuclear volume relative to each ploidy state and treatment.pdf",device = "pdf")
