library(stringr)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(data.table)
library(zoo)
library(gplots)
library(rstatix)
library(reshape2)
library(FSA)
library(ggsignif)
library(ggpubr)

##### Chemo VS Ploidy #####
setwd(choose.dir())

Bleo <- read.csv("Bleomycin_summary_30.4.23.csv")
Bleo_relative <- Bleo
for (i in unique(Bleo_relative$Sample)){
  for (j in colnames(Bleo_relative)[2:4]){
    Bleo_relative[Bleo_relative$Sample==i,j] <- Bleo_relative[Bleo_relative$Sample==i,j]/
      mean(Bleo_relative[Bleo_relative$Sample==i,j][1:3])
  }
}
Bleo_relative <- data.frame("Sample"=rep(Bleo_relative$Sample,3),
                            "Conc."=rep(Bleo_relative$Conc.,3),
                            "timepoint"=rep(c("Day 2","Day 3","Day 5"),each=390),
                            "Relative_intensity"=c(Bleo_relative$Day.2,Bleo_relative$Day.3,Bleo_relative$Day.5))


Bleo_relative <- data.frame(Bleo_relative,"Ploidy"=rep(c(rep("1n",120),rep("2n",120),rep("3n",150)),3))



FU <- read.csv("5-FU_summary_30.4.23.csv")
FU_relative <- FU
for (i in unique(FU_relative$Sample)){
  for (j in colnames(FU_relative)[2:4]){
    FU_relative[FU_relative$Sample==i,j] <- FU_relative[FU_relative$Sample==i,j]/
      mean(FU_relative[FU_relative$Sample==i,j][1:3])
  }
}
FU_relative <- data.frame("Sample"=rep(FU_relative$Sample,3),
                            "Conc."=rep(FU_relative$Conc.,3),
                            "timepoint"=rep(c("Day 2","Day 3","Day 5"),each=390),
                            "Relative_intensity"=c(FU_relative$Day.2,FU_relative$Day.3,FU_relative$Day.5))


FU_relative <- data.frame(FU_relative,"Ploidy"=rep(c(rep("1n",120),rep("2n",120),rep("3n",150)),3))



Cisplatin <- read.csv("Cisplatin_summary_30.4.23.csv")
Cisplatin_relative <- Cisplatin
for (i in unique(Cisplatin_relative$Sample)){
  for (j in colnames(Cisplatin_relative)[2:4]){
    Cisplatin_relative[Cisplatin_relative$Sample==i,j] <- Cisplatin_relative[Cisplatin_relative$Sample==i,j]/
      mean(Cisplatin_relative[Cisplatin_relative$Sample==i,j][1:3])
  }
}
Cisplatin_relative <- data.frame("Sample"=rep(Cisplatin_relative$Sample,3),
                          "Conc."=rep(Cisplatin_relative$Conc.,3),
                          "timepoint"=rep(c("Day 2","Day 3","Day 5"),each=390),
                          "Relative_intensity"=c(Cisplatin_relative$Day.2,Cisplatin_relative$Day.3,Cisplatin_relative$Day.5))


Cisplatin_relative <- data.frame(Cisplatin_relative,"Ploidy"=rep(c(rep("1n",120),rep("2n",120),rep("3n",150)),3))



Paclitaxel <- read.csv("Paclitaxel_summary_30.4.23.csv")
Paclitaxel_relative <- Paclitaxel
for (i in unique(Paclitaxel_relative$Sample)){
  for (j in colnames(Paclitaxel_relative)[2:4]){
    Paclitaxel_relative[Paclitaxel_relative$Sample==i,j] <- Paclitaxel_relative[Paclitaxel_relative$Sample==i,j]/
      mean(Paclitaxel_relative[Paclitaxel_relative$Sample==i,j][1:3])
  }
}
Paclitaxel_relative <- data.frame("Sample"=rep(Paclitaxel_relative$Sample,3),
                                 "Conc."=rep(Paclitaxel_relative$Conc.,3),
                                 "timepoint"=rep(c("Day 2","Day 3","Day 5"),each=390),
                                 "Relative_intensity"=c(Paclitaxel_relative$Day.2,Paclitaxel_relative$Day.3,Paclitaxel_relative$Day.5))


Paclitaxel_relative <- data.frame(Paclitaxel_relative,"Ploidy"=rep(c(rep("1n",120),rep("2n",120),rep("3n",150)),3))



###### comparison between samples ######

## Bleomycin ##
Bleo_day_2 <- Bleo_relative[Bleo_relative$timepoint=="Day 2",]


Bleo_day_3 <- Bleo_relative[Bleo_relative$timepoint=="Day 3",]

Bleo_day_5 <- Bleo_relative[Bleo_relative$timepoint=="Day 5",]


## 5-FU ##
FU_day_2 <- FU_relative[FU_relative$timepoint=="Day 2",]


FU_day_3 <- FU_relative[FU_relative$timepoint=="Day 3",]

FU_day_5 <- FU_relative[FU_relative$timepoint=="Day 5",]


## Cisplatin ##
Cis_day_2 <- Cisplatin_relative[Cisplatin_relative$timepoint=="Day 2",]


Cis_day_3 <- Cisplatin_relative[Cisplatin_relative$timepoint=="Day 3",]

Cis_day_5 <- Cisplatin_relative[Cisplatin_relative$timepoint=="Day 5",]


## Paclitaxel ##
Pacli_day_2 <- Paclitaxel_relative[Paclitaxel_relative$timepoint=="Day 2",]


Pacli_day_3 <- Paclitaxel_relative[Paclitaxel_relative$timepoint=="Day 3",]

Pacli_day_5 <- Paclitaxel_relative[Paclitaxel_relative$timepoint=="Day 5",]


###### comparison between ploidies (mean) ######

## Bleomycin ##
Bleo_day_2_RI <- data.frame("Sample"=rep(unique(Bleo_day_2$Sample),each=10),
                            "Conc."=rep(unique(Bleo_day_2$Conc.),13),
                            "Mean_relative_intensity"=rep(0,130),
                            "Ploidy"=c(rep("1n",40),rep("2n",40), rep("3n",50)))
Bleo_day_2$Conc.[which(Bleo_day_2$Conc.==0)] <- 0.0006
Bleo_day_2_RI$Conc.[which(Bleo_day_2_RI$Conc.==0)] <- 0.0006
for (i in unique(Bleo_day_2_RI$Sample)){
  for (j in unique(Bleo_day_2_RI$Conc.)){
    Bleo_day_2_RI[Bleo_day_2_RI$Sample==i & Bleo_day_2_RI$Conc.==j,3] <- mean(Bleo_day_2[Bleo_day_2$Sample==i & Bleo_day_2$Conc.==j,4])
  }
}
Bleo_day_2_RI <- Bleo_day_2_RI[Bleo_day_2_RI$Sample!="2n 10R p40",]
stat <- Bleo_day_2_RI[Bleo_day_2_RI$Conc.>0.001,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(Bleo_day_2_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=1.3) +
  stat_summary(fun=mean,geom='point') + ggtitle("Bleomycin day 2")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.1,1.15),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(Bleo_day_2_RI$Conc.))[2:10]),breaks = unique(as.numeric(Bleo_day_2_RI$Conc.)))
# ggsave("WT- Bleo day 2 by ploidy.pdf",device = "pdf")

Bleo_day_3_RI <- data.frame("Sample"=rep(unique(Bleo_day_3$Sample),each=10),
                            "Conc."=rep(unique(Bleo_day_3$Conc.),13),
                            "Mean_relative_intensity"=rep(0,130),
                            "Ploidy"=c(rep("1n",40),rep("2n",40), rep("3n",50)))
Bleo_day_3$Conc.[which(Bleo_day_3$Conc.==0)] <- 0.0006
Bleo_day_3_RI$Conc.[which(Bleo_day_3_RI$Conc.==0)] <- 0.0006
for (i in unique(Bleo_day_3$Sample)){
  for (j in unique(Bleo_day_3$Conc.)){
    Bleo_day_3_RI[Bleo_day_3_RI$Sample==i & Bleo_day_3_RI$Conc.==j,3] <- mean(Bleo_day_3[Bleo_day_3$Sample==i & Bleo_day_3$Conc.==j,4])
  }
}
Bleo_day_3_RI <- Bleo_day_3_RI[Bleo_day_3_RI$Sample!="2n 10R p40",]
stat <- Bleo_day_3_RI[Bleo_day_3_RI$Conc.>0.001,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(Bleo_day_3_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=1.3) +
  stat_summary(fun=mean,geom='point') + ggtitle("Bleomycin day 3")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.1,1.15),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(Bleo_day_3_RI$Conc.))[2:10]),breaks = unique(as.numeric(Bleo_day_3_RI$Conc.)))
# ggsave("WT- Bleo day 3 by ploidy.pdf",device = "pdf")

Bleo_day_5_RI <- data.frame("Sample"=rep(unique(Bleo_day_5$Sample),each=10),
                            "Conc."=rep(unique(Bleo_day_5$Conc.),13),
                            "Mean_relative_intensity"=rep(0,130),
                            "Ploidy"=c(rep("1n",40),rep("2n",40), rep("3n",50)))
Bleo_day_5$Conc.[which(Bleo_day_5$Conc.==0)] <- 0.0006        ## assigning a false & close to zero value just to display the concentration of zero in the log2-transformed axis
Bleo_day_5_RI$Conc.[which(Bleo_day_5_RI$Conc.==0)] <- 0.0006
for (i in unique(Bleo_day_5$Sample)){
  for (j in unique(Bleo_day_5$Conc.)){
    Bleo_day_5_RI[Bleo_day_5_RI$Sample==i & Bleo_day_5_RI$Conc.==j,3] <- mean(Bleo_day_5[Bleo_day_5$Sample==i & Bleo_day_5$Conc.==j,4])
  }
}
Bleo_day_5_RI <- Bleo_day_5_RI[Bleo_day_5_RI$Sample!="2n 10R p40",]

stat <- Bleo_day_5_RI[Bleo_day_5_RI$Conc.>0.001,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(Bleo_day_5_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=1.3) +
  stat_summary(fun=mean,geom='point') + ggtitle("Bleomycin day 5")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.3,1.35),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(Bleo_day_5_RI$Conc.))[2:10]),breaks = unique(as.numeric(Bleo_day_5_RI$Conc.)))
# ggsave("WT- Bleo day 5 by ploidy.pdf",device = "pdf")

### calculate and compare IC50 and AUC between ploidies ###
library(bayestestR)
formula = y ~ poly(x, 2, raw = TRUE)

Bleo_D2=ggplot(Bleo_day_2_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","green","forestgreen","darkgreen","deepskyblue","dodgerblue","darkblue","lightcoral","tomato","red","firebrick","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 0.2,formula = formula)+
  ggtitle("Bleomycin - day 2")+scale_x_continuous(limits = c(0,0.4))
build=ggplot_build(Bleo_D2)
Bleo_D2+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")
eq2=str_remove_all(str_remove_all(str_remove_all(str_remove_all(build[["data"]][[3]][["eq.label"]],pattern = "italic\\("),"\\)"),"~"),"`")
match_names <- c("1n 10R p43","1n 10G p43","1n 10G p42","1n 10R p40",
                 "2n 10G p45","2n 10R p48","2n 10G p47",
                 "3n clone H f+8","3n clone E f+6","3n clone C f+7","3n clone B f+7","3n clone J f+7")
names(eq2) <- match_names
coeff <- str_replace_all(str_replace_all(str_remove(str_remove(str_remove(str_remove_all(str_replace(str_remove(string = eq2,"y="),pattern = "- ",replacement = "-"),"\\*x"),"\\+ "),"\\*"),"\\^2"),"%%","*"),"- ","-")
coeff<- as.numeric(str_replace_all(str_remove_all(unlist(strsplit(coeff," ")),"\\*"),pattern = "10\\^",replacement = "e"))
ic50D2 <- c()
for (i in c(1,4,7,10,13,16,19,22,25,28,31,34)){
  ic50D2 <- c(ic50D2,Re(polyroot(c((coeff[i]-0.5),coeff[i+1],coeff[i+2]))))
}
ic50D2 <- data.frame("Sample"=match_names,"IC50.1"=ic50D2[c(1,3,5,7,9,11,13,15,17,19,21,23)],"IC50.2"=ic50D2[c(2,4,6,8,10,12,14,16,18,20,22,24)])

aucD2 <- c()
for (i in unique(Bleo_day_2_RI$Sample)){
  aucD2 <- c(aucD2,auc(Bleo_day_2_RI[Bleo_day_2_RI$Sample==i,2],Bleo_day_2_RI[Bleo_day_2_RI$Sample==i,3]))
}
names(aucD2) <- unique(Bleo_day_2_RI$Sample)


Bleo_D3=ggplot(Bleo_day_3_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","green","forestgreen","darkgreen","deepskyblue","dodgerblue","darkblue","lightcoral","tomato","red","firebrick","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 0.2,formula = formula)+
  ggtitle("Bleomycin - day 3")+scale_x_continuous(limits = c(0,0.4))
build=ggplot_build(Bleo_D3)
Bleo_D3+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")
eq3=str_remove_all(str_remove_all(str_remove_all(str_remove_all(build[["data"]][[3]][["eq.label"]],pattern = "italic\\("),"\\)"),"~"),"`")
match_names <- c("1n 10R p43","1n 10G p43","1n 10G p42","1n 10R p40",
                 "2n 10G p45","2n 10R p48","2n 10G p47",
                 "3n clone H f+8","3n clone E f+6","3n clone C f+7","3n clone B f+7","3n clone J f+7")
names(eq3) <- match_names
coeff <- str_replace_all(str_replace_all(str_remove(str_remove(str_remove(str_remove_all(str_replace(str_remove(string = eq3,"y="),pattern = "- ",replacement = "-"),"\\*x"),"\\+ "),"\\*"),"\\^2"),"%%","*"),"- ","-")
coeff<- as.numeric(str_replace_all(str_remove_all(unlist(strsplit(coeff," ")),"\\*"),pattern = "10\\^",replacement = "e"))
ic50D3 <- c()
for (i in c(1,4,7,10,13,16,19,22,25,28,31,34)){
  ic50D3 <- c(ic50D3,Re(polyroot(c((coeff[i]-0.5),coeff[i+1],coeff[i+2]))))
}
ic50D3 <- data.frame("Sample"=match_names,"IC50.1"=ic50D3[c(1,3,5,7,9,11,13,15,17,19,21,23)],"IC50.2"=ic50D3[c(2,4,6,8,10,12,14,16,18,20,22,24)])

aucD3 <- c()
for (i in unique(Bleo_day_3_RI$Sample)){
  aucD3 <- c(aucD3,auc(Bleo_day_3_RI[Bleo_day_3_RI$Sample==i,2],Bleo_day_3_RI[Bleo_day_3_RI$Sample==i,3]))
}
names(aucD3) <- unique(Bleo_day_3_RI$Sample)

Bleo_D5=ggplot(Bleo_day_5_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","green","forestgreen","darkgreen","deepskyblue","dodgerblue","darkblue","lightcoral","tomato","red","firebrick","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 0.2,formula = formula)+
  ggtitle("Bleomycin - day 5")+scale_x_continuous(limits = c(0,0.4))
build=ggplot_build(Bleo_D5)
Bleo_D5+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")
eq5=str_remove_all(str_remove_all(str_remove_all(str_remove_all(build[["data"]][[3]][["eq.label"]],pattern = "italic\\("),"\\)"),"~"),"`")
match_names <- c("1n 10R p43","1n 10G p43","1n 10G p42","1n 10R p40",
                 "2n 10G p45","2n 10R p48","2n 10G p47",
                 "3n clone H f+8","3n clone E f+6","3n clone C f+7","3n clone B f+7","3n clone J f+7")
names(eq5) <- match_names
coeff <- str_replace_all(str_replace_all(str_remove(str_remove(str_remove(str_remove_all(str_replace(str_remove(string = eq5,"y="),pattern = "- ",replacement = "-"),"\\*x"),"\\+ "),"\\*"),"\\^2"),"%%","*"),"- ","-")
coeff<- as.numeric(str_replace_all(str_remove_all(unlist(strsplit(coeff," ")),"\\*"),pattern = "10\\^",replacement = "e"))
ic50D5 <- c()
for (i in c(1,4,7,10,13,16,19,22,25,28,31,34)){
  ic50D5 <- c(ic50D5,Re(polyroot(c((coeff[i]-0.5),coeff[i+1],coeff[i+2]))))
}
ic50D5 <- data.frame("Sample"=match_names,"IC50.1"=ic50D5[c(1,3,5,7,9,11,13,15,17,19,21,23)],"IC50.2"=ic50D5[c(2,4,6,8,10,12,14,16,18,20,22,24)])

aucD5 <- c()
for (i in unique(Bleo_day_5_RI$Sample)){
  aucD5 <- c(aucD5,auc(Bleo_day_5_RI[Bleo_day_5_RI$Sample==i,2],Bleo_day_5_RI[Bleo_day_5_RI$Sample==i,3]))
}
names(aucD5) <- unique(Bleo_day_5_RI$Sample)


Bleo_IC50 <- rbind(ic50D2,ic50D3,ic50D5)
IC50 <- c()
for (i in 1:36){
  if (Bleo_IC50[i,2]>0){
    IC50[i] <- Bleo_IC50[i,2]
  }
  else{
    IC50[i] <- Bleo_IC50[i,3]
  }
  
}
Bleo_IC50 <- data.frame("Sample"=Bleo_IC50$Sample,IC50)
ploidy <- rep(0,36)
ploidy[grepl("1n",Bleo_IC50$Sample)] <- "1n"
ploidy[grepl("2n",Bleo_IC50$Sample)] <- "2n"
ploidy[grepl("3n",Bleo_IC50$Sample)] <- "3n"
Timepoint <- rep(c("Day 2","Day 3","Day 5"),each=12)
Bleo_IC50 <- data.frame(Bleo_IC50,"Ploidy"=ploidy,"Timepoint"=Timepoint)
write.csv(Bleo_IC50,"Bleo_IC50 - WT.csv")
Bleo_AUC <- data.frame("Sample"=c(names(aucD2),names(aucD3),names(aucD5)),"AUC"=c(aucD2,aucD3,aucD5),"Ploidy"=c("1n","1n","1n","1n","2n","2n","2n","3n","3n","3n","3n","3n"),"Timepoint"=rep(c("Day 2","Day 3","Day 5"),each=12))
write.csv(Bleo_AUC,"Bleo_AUC - WT.csv")


#mean IC50 (day 5):
stat <- Bleo_IC50[Bleo_IC50$Timepoint=="Day 5",] %>%
  t_test(IC50 ~ Ploidy, alternative = "greater", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Ploidy", dodge = 0.8)
stat$p.adj <- format(stat$p.adj,scientific=T)
ggbarplot(Bleo_IC50[Bleo_IC50$Timepoint=="Day 5",], x= "Ploidy", y= "IC50",fill="Ploidy", position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat,label = "p.adj",y.position = c(0.45,0.5,0.3), size = 3,tip.length = 0.01, hide.ns = F)+
  scale_fill_manual(values=c("forestgreen","royalblue1","red1"))+ggtitle("Bleomycin IC50 Vs. Ploidy")+
  theme(legend.position = "none")+coord_cartesian(ylim = c(0,0.55),xlim = c(0.5,3.5),expand = F)
# ggsave("barplot of Bleomycin mean IC50 Vs. ploidy of WT hESCs - day 5.pdf",device = "pdf")


#mean AUC:
stat <- Bleo_AUC %>%
  t_test(AUC ~ Ploidy, p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Ploidy", dodge = 0.8)
stat$p.adj <- format(stat$p.adj,scientific=T)
ggbarplot(Bleo_AUC, x= "Ploidy", y= "AUC",fill="Ploidy", position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat,label = "p.adj",y.position = c(0.26,0.28,0.17), size = 3,tip.length = 0.01, hide.ns = F)+
  scale_fill_manual(values=c("forestgreen","royalblue1","red1"))+ggtitle("Bleomycin AUC Vs. Ploidy")+scale_y_continuous(expand = c(0,0))
# ggsave("barplot of Bleomycin mean AUC Vs. ploidy of WT hESCs.pdf",device = "pdf")



## 5-FU ##
FU_day_2_RI <- data.frame("Sample"=rep(unique(FU_day_2$Sample),each=10),
                          "Conc."=rep(unique(FU_day_2$Conc.),13),
                          "Mean_relative_intensity"=rep(0,130),
                          "Ploidy"=c(rep("1n",40),rep("2n",40), rep("3n",50)))
FU_day_2$Conc.[which(FU_day_2$Conc.==0)] <- 0.0195
FU_day_2_RI$Conc.[which(FU_day_2_RI$Conc.==0)] <- 0.0195
for (i in unique(FU_day_2_RI$Sample)){
  for (j in unique(FU_day_2_RI$Conc.)){
    FU_day_2_RI[FU_day_2_RI$Sample==i & FU_day_2_RI$Conc.==j,3] <- mean(FU_day_2[FU_day_2$Sample==i & FU_day_2$Conc.==j,4])
  }
}
FU_day_2_RI <- FU_day_2_RI[FU_day_2_RI$Sample!="2n 10R p40",]
stat <- FU_day_2_RI[FU_day_2_RI$Conc.>0.02,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(FU_day_2_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=1.3) +
  stat_summary(fun=mean,geom='point') + ggtitle("5-FU day 2")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.1,1.15),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(FU_day_2_RI$Conc.))[2:10]),breaks = unique(as.numeric(FU_day_2_RI$Conc.)))
# ggsave("WT- FU day 2 by ploidy.pdf",device = "pdf")

FU_day_3_RI <- data.frame("Sample"=rep(unique(FU_day_3$Sample),each=10),
                          "Conc."=rep(unique(FU_day_3$Conc.),13),
                          "Mean_relative_intensity"=rep(0,130),
                          "Ploidy"=c(rep("1n",40),rep("2n",40), rep("3n",50)))
FU_day_3$Conc.[which(FU_day_3$Conc.==0)] <- 0.0195
FU_day_3_RI$Conc.[which(FU_day_3_RI$Conc.==0)] <- 0.0195
for (i in unique(FU_day_3$Sample)){
  for (j in unique(FU_day_3$Conc.)){
    FU_day_3_RI[FU_day_3_RI$Sample==i & FU_day_3_RI$Conc.==j,3] <- mean(FU_day_3[FU_day_3$Sample==i & FU_day_3$Conc.==j,4])
  }
}
FU_day_3_RI <- FU_day_3_RI[FU_day_3_RI$Sample!="2n 10R p40",]
stat <- FU_day_3_RI[FU_day_3_RI$Conc.>0.02,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(FU_day_3_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=1.3) +
  stat_summary(fun=mean,geom='point') + ggtitle("5-FU day 3")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.1,1.15),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(FU_day_3_RI$Conc.))[2:10]),breaks = unique(as.numeric(FU_day_3_RI$Conc.)))
# ggsave("WT- FU day 3 by ploidy.pdf",device = "pdf")

FU_day_5_RI <- data.frame("Sample"=rep(unique(FU_day_5$Sample),each=10),
                          "Conc."=rep(unique(FU_day_5$Conc.),13),
                          "Mean_relative_intensity"=rep(0,130),
                          "Ploidy"=c(rep("1n",40),rep("2n",40), rep("3n",50)))
FU_day_5$Conc.[which(FU_day_5$Conc.==0)] <- 0.0195
FU_day_5_RI$Conc.[which(FU_day_5_RI$Conc.==0)] <- 0.0195
for (i in unique(FU_day_5$Sample)){
  for (j in unique(FU_day_5$Conc.)){
    FU_day_5_RI[FU_day_5_RI$Sample==i & FU_day_5_RI$Conc.==j,3] <- mean(FU_day_5[FU_day_5$Sample==i & FU_day_5$Conc.==j,4])
  }
}
FU_day_5_RI <- FU_day_5_RI[FU_day_5_RI$Sample!="2n 10R p40",]

stat <- FU_day_5_RI[FU_day_5_RI$Conc.>0.02,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(FU_day_5_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=1.3) +
  stat_summary(fun=mean,geom='point') + ggtitle("5-FU day 5")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.3,1.35),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(FU_day_5_RI$Conc.))[2:10]),breaks = unique(as.numeric(FU_day_5_RI$Conc.)))
# ggsave("WT- FU day 5 by ploidy.pdf",device = "pdf")

### calculate and compare IC50 and AUC between ploidies ###
library(bayestestR)
formula = y ~ poly(x, 2, raw = TRUE)

FU_D2=ggplot(FU_day_2_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","green","forestgreen","darkgreen","deepskyblue","dodgerblue","darkblue","lightcoral","tomato","red","firebrick","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 7,formula = formula)+
  ggtitle("5-FU - day 2")+scale_x_continuous(limits = c(0,12))
build=ggplot_build(FU_D2)
FU_D2+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")
eq2=str_remove_all(str_remove_all(str_remove_all(str_remove_all(build[["data"]][[3]][["eq.label"]],pattern = "italic\\("),"\\)"),"~"),"`")
match_names <- c("1n 10R p43","1n 10G p43","1n 10G p42","1n 10R p40",
                 "2n 10G p45","2n 10R p48","2n 10G p47",
                 "3n clone H f+8","3n clone E f+6","3n clone C f+7","3n clone B f+7","3n clone J f+7")
names(eq2) <- match_names
coeff <- str_replace_all(str_replace_all(str_remove(str_remove(str_remove(str_remove_all(str_replace(str_remove(string = eq2,"y="),pattern = "- ",replacement = "-"),"\\*x"),"\\+ "),"\\*"),"\\^2"),"%%","*"),"- ","-")
coeff<- as.numeric(str_replace_all(str_remove_all(unlist(strsplit(coeff," ")),"\\*"),pattern = "10\\^",replacement = "e"))
ic50D2 <- c()
for (i in c(1,4,7,10,13,16,19,22,25,28,31,34)){
  ic50D2 <- c(ic50D2,Re(polyroot(c((coeff[i]-0.5),coeff[i+1],coeff[i+2]))))
}
ic50D2 <- data.frame("Sample"=match_names,"IC50.1"=ic50D2[c(1,3,5,7,9,11,13,15,17,19,21,23)],"IC50.2"=ic50D2[c(2,4,6,8,10,12,14,16,18,20,22,24)])

aucD2 <- c()
for (i in unique(FU_day_2_RI$Sample)){
  aucD2 <- c(aucD2,auc(FU_day_2_RI[FU_day_2_RI$Sample==i,2],FU_day_2_RI[FU_day_2_RI$Sample==i,3]))
}
names(aucD2) <- unique(FU_day_2_RI$Sample)


FU_D3=ggplot(FU_day_3_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","green","forestgreen","darkgreen","deepskyblue","dodgerblue","darkblue","lightcoral","tomato","red","firebrick","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 7,formula = formula)+
  ggtitle("5-FU - day 3")+scale_x_continuous(limits = c(0,13))
build=ggplot_build(FU_D3)
FU_D3+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")
eq3=str_remove_all(str_remove_all(str_remove_all(str_remove_all(build[["data"]][[3]][["eq.label"]],pattern = "italic\\("),"\\)"),"~"),"`")
match_names <- c("1n 10R p43","1n 10G p43","1n 10G p42","1n 10R p40",
                 "2n 10G p45","2n 10R p48","2n 10G p47",
                 "3n clone H f+8","3n clone E f+6","3n clone C f+7","3n clone B f+7","3n clone J f+7")
names(eq3) <- match_names
coeff <- str_replace_all(str_replace_all(str_remove(str_remove(str_remove(str_remove_all(str_replace(str_remove(string = eq3,"y="),pattern = "- ",replacement = "-"),"\\*x"),"\\+ "),"\\*"),"\\^2"),"%%","*"),"- ","-")
coeff<- as.numeric(str_replace_all(str_remove_all(unlist(strsplit(coeff," ")),"\\*"),pattern = "10\\^",replacement = "e"))
ic50D3 <- c()
for (i in c(1,4,7,10,13,16,19,22,25,28,31,34)){
  ic50D3 <- c(ic50D3,Re(polyroot(c((coeff[i]-0.5),coeff[i+1],coeff[i+2]))))
}
ic50D3 <- data.frame("Sample"=match_names,"IC50.1"=ic50D3[c(1,3,5,7,9,11,13,15,17,19,21,23)],"IC50.2"=ic50D3[c(2,4,6,8,10,12,14,16,18,20,22,24)])

aucD3 <- c()
for (i in unique(FU_day_3_RI$Sample)){
  aucD3 <- c(aucD3,auc(FU_day_3_RI[FU_day_3_RI$Sample==i,2],FU_day_3_RI[FU_day_3_RI$Sample==i,3]))
}
names(aucD3) <- unique(FU_day_3_RI$Sample)

FU_D5=ggplot(FU_day_5_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","green","forestgreen","darkgreen","deepskyblue","dodgerblue","darkblue","lightcoral","tomato","red","firebrick","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 7,formula = formula)+
  ggtitle("5-FU - day 5")+scale_x_continuous(limits = c(0,13))
build=ggplot_build(FU_D5)
FU_D5+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")
eq5=str_remove_all(str_remove_all(str_remove_all(str_remove_all(build[["data"]][[3]][["eq.label"]],pattern = "italic\\("),"\\)"),"~"),"`")
match_names <- c("1n 10R p43","1n 10G p43","1n 10G p42","1n 10R p40",
                 "2n 10G p45","2n 10R p48","2n 10G p47",
                 "3n clone H f+8","3n clone E f+6","3n clone C f+7","3n clone B f+7","3n clone J f+7")
names(eq5) <- match_names
coeff <- str_replace_all(str_replace_all(str_remove(str_remove(str_remove(str_remove_all(str_replace(str_remove(string = eq5,"y="),pattern = "- ",replacement = "-"),"\\*x"),"\\+ "),"\\*"),"\\^2"),"%%","*"),"- ","-")
coeff<- as.numeric(str_replace_all(str_remove_all(unlist(strsplit(coeff," ")),"\\*"),pattern = "10\\^",replacement = "e"))
ic50D5 <- c()
for (i in c(1,4,7,10,13,16,19,22,25,28,31,34)){
  ic50D5 <- c(ic50D5,Re(polyroot(c((coeff[i]-0.5),coeff[i+1],coeff[i+2]))))
}
ic50D5 <- data.frame("Sample"=match_names,"IC50.1"=ic50D5[c(1,3,5,7,9,11,13,15,17,19,21,23)],"IC50.2"=ic50D5[c(2,4,6,8,10,12,14,16,18,20,22,24)])

aucD5 <- c()
for (i in unique(FU_day_5_RI$Sample)){
  aucD5 <- c(aucD5,auc(FU_day_5_RI[FU_day_5_RI$Sample==i,2],FU_day_5_RI[FU_day_5_RI$Sample==i,3]))
}
names(aucD5) <- unique(FU_day_5_RI$Sample)


FU_IC50 <- rbind(ic50D2,ic50D3,ic50D5)
FU_IC50 <- FU_IC50[,1:2]  #all the second IC50 values are beyond the max concentration we used
ploidy <- rep(0,36)
ploidy[grepl("1n",FU_IC50$Sample)] <- "1n"
ploidy[grepl("2n",FU_IC50$Sample)] <- "2n"
ploidy[grepl("3n",FU_IC50$Sample)] <- "3n"
Timepoint <- rep(c("Day 2","Day 3","Day 5"),each=12)
FU_IC50 <- data.frame(FU_IC50,"Ploidy"=ploidy,"Timepoint"=Timepoint)
colnames(FU_IC50)[2] <- "IC50"
write.csv(FU_IC50,"FU_IC50 - WT.csv")
FU_AUC <- data.frame("Sample"=c(names(aucD2),names(aucD3),names(aucD5)),"AUC"=c(aucD2,aucD3,aucD5),"Ploidy"=c("1n","1n","1n","1n","2n","2n","2n","3n","3n","3n","3n","3n"),"Timepoint"=rep(c("Day 2","Day 3","Day 5"),each=12))
write.csv(FU_AUC,"FU_AUC - WT.csv")


#mean IC50 (day 2):
stat <- FU_IC50[FU_IC50$Timepoint=="Day 2",] %>%
  t_test(IC50 ~ Ploidy, alternative = "greater", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Ploidy", dodge = 0.8)
stat$p.adj <- format(stat$p.adj,scientific=T)
ggbarplot(FU_IC50[FU_IC50$Timepoint=="Day 2",], x= "Ploidy", y= "IC50",fill="Ploidy", position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat,label = "p.adj",y.position = c(7.6,8,7.2), size = 3,tip.length = 0.01, hide.ns = F)+
  scale_fill_manual(values=c("forestgreen","royalblue1","red1"))+ggtitle("5-FU IC50 Vs. Ploidy")+
  theme(legend.position = "none")+coord_cartesian(ylim = c(0,8.5),xlim = c(0.5,3.5),expand = F)
# ggsave("barplot of 5-FU mean IC50 Vs. ploidy of WT hESCs - day 2.pdf",device = "pdf")


#mean AUC:
stat <- FU_AUC %>%
  t_test(AUC ~ Ploidy, alternative = "greater", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Ploidy", dodge = 0.8)
stat$p.adj <- format(stat$p.adj,scientific=T)
ggbarplot(FU_AUC, x= "Ploidy", y= "AUC",fill="Ploidy", position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat,label = "p.adj",y.position = c(5.5,5.8,5.2), size = 3,tip.length = 0.01, hide.ns = F)+
  scale_fill_manual(values=c("forestgreen","royalblue1","red1"))+ggtitle("5-FU AUC Vs. Ploidy")+coord_cartesian(ylim = c(0,6),xlim = c(0.5,3.5),expand = F)
# ggsave("barplot of 5-FU mean AUC Vs. ploidy of WT hESCs.pdf",device = "pdf")



## Cisplatin ##
Cis_day_2_RI <- data.frame("Sample"=rep(unique(Cis_day_2$Sample),each=10),
                           "Conc."=rep(unique(Cis_day_2$Conc.),13),
                           "Mean_relative_intensity"=rep(0,130),
                           "Ploidy"=c(rep("1n",40),rep("2n",40), rep("3n",50)))
Cis_day_2$Conc.[which(Cis_day_2$Conc.==0)] <- 0.0125
Cis_day_2_RI$Conc.[which(Cis_day_2_RI$Conc.==0)] <- 0.0125
for (i in unique(Cis_day_2_RI$Sample)){
  for (j in unique(Cis_day_2_RI$Conc.)){
    Cis_day_2_RI[Cis_day_2_RI$Sample==i & Cis_day_2_RI$Conc.==j,3] <- mean(Cis_day_2[Cis_day_2$Sample==i & Cis_day_2$Conc.==j,4])
  }
}
Cis_day_2_RI <- Cis_day_2_RI[Cis_day_2_RI$Sample!="2n 10R p40",]
stat <- Cis_day_2_RI[Cis_day_2_RI$Conc.>0.02,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(Cis_day_2_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=1.3) +
  stat_summary(fun=mean,geom='point') + ggtitle("Cisplatin day 2")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.1,1.15),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(Cis_day_2_RI$Conc.))[2:10]),breaks = unique(as.numeric(Cis_day_2_RI$Conc.)))
# ggsave("WT- Cis day 2 by ploidy.pdf",device = "pdf")

Cis_day_3_RI <- data.frame("Sample"=rep(unique(Cis_day_3$Sample),each=10),
                           "Conc."=rep(unique(Cis_day_3$Conc.),13),
                           "Mean_relative_intensity"=rep(0,130),
                           "Ploidy"=c(rep("1n",40),rep("2n",40), rep("3n",50)))
Cis_day_3$Conc.[which(Cis_day_3$Conc.==0)] <- 0.0125
Cis_day_3_RI$Conc.[which(Cis_day_3_RI$Conc.==0)] <- 0.0125
for (i in unique(Cis_day_3$Sample)){
  for (j in unique(Cis_day_3$Conc.)){
    Cis_day_3_RI[Cis_day_3_RI$Sample==i & Cis_day_3_RI$Conc.==j,3] <- mean(Cis_day_3[Cis_day_3$Sample==i & Cis_day_3$Conc.==j,4])
  }
}
Cis_day_3_RI <- Cis_day_3_RI[Cis_day_3_RI$Sample!="2n 10R p40",]
stat <- Cis_day_3_RI[Cis_day_3_RI$Conc.>0.02,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(Cis_day_3_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=1.3) +
  stat_summary(fun=mean,geom='point') + ggtitle("Cisplatin day 3")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.1,1.15),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(Cis_day_3_RI$Conc.))[2:10]),breaks = unique(as.numeric(Cis_day_3_RI$Conc.)))
# ggsave("WT- Cis day 3 by ploidy.pdf",device = "pdf")

Cis_day_5_RI <- data.frame("Sample"=rep(unique(Cis_day_5$Sample),each=10),
                           "Conc."=rep(unique(Cis_day_5$Conc.),13),
                           "Mean_relative_intensity"=rep(0,130),
                           "Ploidy"=c(rep("1n",40),rep("2n",40), rep("3n",50)))
Cis_day_5$Conc.[which(Cis_day_5$Conc.==0)] <- 0.0125
Cis_day_5_RI$Conc.[which(Cis_day_5_RI$Conc.==0)] <- 0.0125
for (i in unique(Cis_day_5$Sample)){
  for (j in unique(Cis_day_5$Conc.)){
    Cis_day_5_RI[Cis_day_5_RI$Sample==i & Cis_day_5_RI$Conc.==j,3] <- mean(Cis_day_5[Cis_day_5$Sample==i & Cis_day_5$Conc.==j,4])
  }
}
Cis_day_5_RI <- Cis_day_5_RI[Cis_day_5_RI$Sample!="2n 10R p40",]

stat <- Cis_day_5_RI[Cis_day_5_RI$Conc.>0.02,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(Cis_day_5_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=1.3) +
  stat_summary(fun=mean,geom='point') + ggtitle("Cisplatin day 5")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.3,1.35),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(Cis_day_5_RI$Conc.))[2:10]),breaks = unique(as.numeric(Cis_day_5_RI$Conc.)))
# ggsave("WT- Cis day 5 by ploidy.pdf",device = "pdf")

### calculate and compare IC50 and AUC between ploidies ###
library(bayestestR)
formula = y ~ poly(x, 2, raw = TRUE)

Cis_D2=ggplot(Cis_day_2_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","green","forestgreen","darkgreen","deepskyblue","dodgerblue","darkblue","lightcoral","tomato","red","firebrick","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 0.5,formula = formula)+
  ggtitle("Cisplatin - day 2")+scale_x_continuous(limits = c(0,1))
build=ggplot_build(Cis_D2)
Cis_D2+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")
eq2=str_remove_all(str_remove_all(str_remove_all(str_remove_all(build[["data"]][[3]][["eq.label"]],pattern = "italic\\("),"\\)"),"~"),"`")
match_names <- c("1n 10R p43","1n 10G p43","1n 10G p42","1n 10R p40",
                 "2n 10G p45","2n 10R p48","2n 10G p47",
                 "3n clone H f+8","3n clone E f+6","3n clone C f+7","3n clone B f+7","3n clone J f+7")
names(eq2) <- match_names
coeff <- str_replace_all(str_replace_all(str_remove(str_remove(str_remove(str_remove_all(str_replace(str_remove(string = eq2,"y="),pattern = "- ",replacement = "-"),"\\*x"),"\\+ "),"\\*"),"\\^2"),"%%","*"),"- ","-")
coeff<- as.numeric(str_replace_all(str_remove_all(unlist(strsplit(coeff," ")),"\\*"),pattern = "10\\^",replacement = "e"))
ic50D2 <- c()
for (i in c(1,4,7,10,13,16,19,22,25,28,31,34)){
  ic50D2 <- c(ic50D2,Re(polyroot(c((coeff[i]-0.5),coeff[i+1],coeff[i+2]))))
}
ic50D2 <- data.frame("Sample"=match_names,"IC50.1"=ic50D2[c(1,3,5,7,9,11,13,15,17,19,21,23)],"IC50.2"=ic50D2[c(2,4,6,8,10,12,14,16,18,20,22,24)])

aucD2 <- c()
for (i in unique(Cis_day_2_RI$Sample)){
  aucD2 <- c(aucD2,auc(Cis_day_2_RI[Cis_day_2_RI$Sample==i,2],Cis_day_2_RI[Cis_day_2_RI$Sample==i,3]))
}
names(aucD2) <- unique(Cis_day_2_RI$Sample)


Cis_D3=ggplot(Cis_day_3_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","green","forestgreen","darkgreen","deepskyblue","dodgerblue","darkblue","lightcoral","tomato","red","firebrick","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 0.5,formula = formula)+
  ggtitle("Cisplatin - day 3")+scale_x_continuous(limits = c(0,1))
build=ggplot_build(Cis_D3)
Cis_D3+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")
eq3=str_remove_all(str_remove_all(str_remove_all(str_remove_all(build[["data"]][[3]][["eq.label"]],pattern = "italic\\("),"\\)"),"~"),"`")
match_names <- c("1n 10R p43","1n 10G p43","1n 10G p42","1n 10R p40",
                 "2n 10G p45","2n 10R p48","2n 10G p47",
                 "3n clone H f+8","3n clone E f+6","3n clone C f+7","3n clone B f+7","3n clone J f+7")
names(eq3) <- match_names
coeff <- str_replace_all(str_replace_all(str_remove(str_remove(str_remove(str_remove_all(str_replace(str_remove(string = eq3,"y="),pattern = "- ",replacement = "-"),"\\*x"),"\\+ "),"\\*"),"\\^2"),"%%","*"),"- ","-")
coeff<- as.numeric(str_replace_all(str_remove_all(unlist(strsplit(coeff," ")),"\\*"),pattern = "10\\^",replacement = "e"))
ic50D3 <- c()
for (i in c(1,4,7,10,13,16,19,22,25,28,31,34)){
  ic50D3 <- c(ic50D3,Re(polyroot(c((coeff[i]-0.5),coeff[i+1],coeff[i+2]))))
}
ic50D3 <- data.frame("Sample"=match_names,"IC50.1"=ic50D3[c(1,3,5,7,9,11,13,15,17,19,21,23)],"IC50.2"=ic50D3[c(2,4,6,8,10,12,14,16,18,20,22,24)])

aucD3 <- c()
for (i in unique(Cis_day_3_RI$Sample)){
  aucD3 <- c(aucD3,auc(Cis_day_3_RI[Cis_day_3_RI$Sample==i,2],Cis_day_3_RI[Cis_day_3_RI$Sample==i,3]))
}
names(aucD3) <- unique(Cis_day_3_RI$Sample)

Cis_D5=ggplot(Cis_day_5_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","green","forestgreen","darkgreen","deepskyblue","dodgerblue","darkblue","lightcoral","tomato","red","firebrick","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 0.5,formula = formula)+
  ggtitle("Cisplatin - day 5")+scale_x_continuous(limits = c(0,1))
build=ggplot_build(Cis_D5)
Cis_D5+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")
eq5=str_remove_all(str_remove_all(str_remove_all(str_remove_all(build[["data"]][[3]][["eq.label"]],pattern = "italic\\("),"\\)"),"~"),"`")
match_names <- c("1n 10R p43","1n 10G p43","1n 10G p42","1n 10R p40",
                 "2n 10G p45","2n 10R p48","2n 10G p47",
                 "3n clone H f+8","3n clone E f+6","3n clone C f+7","3n clone B f+7","3n clone J f+7")
names(eq5) <- match_names
coeff <- str_replace_all(str_replace_all(str_remove(str_remove(str_remove(str_remove_all(str_replace(str_remove(string = eq5,"y="),pattern = "- ",replacement = "-"),"\\*x"),"\\+ "),"\\*"),"\\^2"),"%%","*"),"- ","-")
coeff<- as.numeric(str_replace_all(str_remove_all(unlist(strsplit(coeff," ")),"\\*"),pattern = "10\\^",replacement = "e"))
ic50D5 <- c()
for (i in c(1,4,7,10,13,16,19,22,25,28,31,34)){
  ic50D5 <- c(ic50D5,Re(polyroot(c((coeff[i]-0.5),coeff[i+1],coeff[i+2]))))
}
ic50D5 <- data.frame("Sample"=match_names,"IC50.1"=ic50D5[c(1,3,5,7,9,11,13,15,17,19,21,23)],"IC50.2"=ic50D5[c(2,4,6,8,10,12,14,16,18,20,22,24)])

aucD5 <- c()
for (i in unique(Cis_day_5_RI$Sample)){
  aucD5 <- c(aucD5,auc(Cis_day_5_RI[Cis_day_5_RI$Sample==i,2],Cis_day_5_RI[Cis_day_5_RI$Sample==i,3]))
}
names(aucD5) <- unique(Cis_day_5_RI$Sample)


Cis_IC50 <- rbind(ic50D2,ic50D3,ic50D5)
Cis_IC50 <- Cis_IC50[,1:2]  #all the second IC50 values are beyond the max concentration we used
ploidy <- rep(0,36)
ploidy[grepl("1n",Cis_IC50$Sample)] <- "1n"
ploidy[grepl("2n",Cis_IC50$Sample)] <- "2n"
ploidy[grepl("3n",Cis_IC50$Sample)] <- "3n"
Timepoint <- rep(c("Day 2","Day 3","Day 5"),each=12)
Cis_IC50 <- data.frame(Cis_IC50,"Ploidy"=ploidy,"Timepoint"=Timepoint)
colnames(Cis_IC50)[2] <- "IC50"
write.csv(Cis_IC50,"Cis_IC50 - WT.csv")
Cis_AUC <- data.frame("Sample"=c(names(aucD2),names(aucD3),names(aucD5)),"AUC"=c(aucD2,aucD3,aucD5),"Ploidy"=c("1n","1n","1n","1n","2n","2n","2n","3n","3n","3n","3n","3n"),"Timepoint"=rep(c("Day 2","Day 3","Day 5"),each=12))
write.csv(Cis_AUC,"Cis_AUC - WT.csv")


#mean IC50 (day 2):
stat <- Cis_IC50[Cis_IC50$Timepoint=="Day 2",] %>%
  t_test(IC50 ~ Ploidy, alternative = "greater", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Ploidy", dodge = 0.8)
stat$p.adj <- format(stat$p.adj,scientific=T)
ggbarplot(Cis_IC50[Cis_IC50$Timepoint=="Day 2",], x= "Ploidy", y= "IC50",fill="Ploidy", position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat,label = "p.adj",y.position = c(0.45,0.48,0.3), size = 3,tip.length = 0.01, hide.ns = F)+
  scale_fill_manual(values=c("forestgreen","royalblue1","red1"))+ggtitle("Cisplatin IC50 Vs. Ploidy")+
  theme(legend.position = "none")+coord_cartesian(ylim = c(0,0.5),xlim = c(0.5,3.5),expand = F)
# ggsave("barplot of Cisplatin mean IC50 Vs. ploidy of WT hESCs - day 2.pdf",device = "pdf")


#mean AUC:
stat <- Cis_AUC %>%
  t_test(AUC ~ Ploidy, p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Ploidy", dodge = 0.8)
stat$p.adj <- format(stat$p.adj,scientific=T)
ggbarplot(Cis_AUC, x= "Ploidy", y= "AUC",fill="Ploidy", position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat,label = "p.adj",y.position = c(0.8,0.85,0.6), size = 3,tip.length = 0.01, hide.ns = F)+
  scale_fill_manual(values=c("forestgreen","royalblue1","red1"))+ggtitle("Cisplatin AUC Vs. Ploidy")+scale_y_continuous(expand = c(0,0))
# ggsave("barplot of Cisplatin mean AUC Vs. ploidy of WT hESCs.pdf",device = "pdf")



## Paclitaxel ##
Pacli_day_2_RI <- data.frame("Sample"=rep(unique(Pacli_day_2$Sample),each=10),
                             "Conc."=rep(unique(Pacli_day_2$Conc.),13),
                             "Mean_relative_intensity"=rep(0,130),
                             "Ploidy"=c(rep("1n",40),rep("2n",40), rep("3n",50)))
Pacli_day_2$Conc.[which(Pacli_day_2$Conc.==0)] <- 0.0125
Pacli_day_2_RI$Conc.[which(Pacli_day_2_RI$Conc.==0)] <- 0.0125
for (i in unique(Pacli_day_2_RI$Sample)){
  for (j in unique(Pacli_day_2_RI$Conc.)){
    Pacli_day_2_RI[Pacli_day_2_RI$Sample==i & Pacli_day_2_RI$Conc.==j,3] <- mean(Pacli_day_2[Pacli_day_2$Sample==i & Pacli_day_2$Conc.==j,4])
  }
}
Pacli_day_2_RI <- Pacli_day_2_RI[Pacli_day_2_RI$Sample!="2n 10R p40",]
stat <- Pacli_day_2_RI[Pacli_day_2_RI$Conc.>0.02,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(Pacli_day_2_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=1.3) +
  stat_summary(fun=mean,geom='point') + ggtitle("Paclitaxel day 2")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.1,1.15),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(Pacli_day_2_RI$Conc.))[2:10]),breaks = unique(as.numeric(Pacli_day_2_RI$Conc.)))
# ggsave("WT- Pacli day 2 by ploidy.pdf",device = "pdf")

Pacli_day_3_RI <- data.frame("Sample"=rep(unique(Pacli_day_3$Sample),each=10),
                             "Conc."=rep(unique(Pacli_day_3$Conc.),13),
                             "Mean_relative_intensity"=rep(0,130),
                             "Ploidy"=c(rep("1n",40),rep("2n",40), rep("3n",50)))
Pacli_day_3$Conc.[which(Pacli_day_3$Conc.==0)] <- 0.0125
Pacli_day_3_RI$Conc.[which(Pacli_day_3_RI$Conc.==0)] <- 0.0125
for (i in unique(Pacli_day_3$Sample)){
  for (j in unique(Pacli_day_3$Conc.)){
    Pacli_day_3_RI[Pacli_day_3_RI$Sample==i & Pacli_day_3_RI$Conc.==j,3] <- mean(Pacli_day_3[Pacli_day_3$Sample==i & Pacli_day_3$Conc.==j,4])
  }
}
Pacli_day_3_RI <- Pacli_day_3_RI[Pacli_day_3_RI$Sample!="2n 10R p40",]
stat <- Pacli_day_3_RI[Pacli_day_3_RI$Conc.>0.02,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(Pacli_day_3_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=1.3) +
  stat_summary(fun=mean,geom='point') + ggtitle("Paclitaxel day 3")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.1,1.15),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(Pacli_day_3_RI$Conc.))[2:10]),breaks = unique(as.numeric(Pacli_day_3_RI$Conc.)))
# ggsave("WT- Pacli day 3 by ploidy.pdf",device = "pdf")

Pacli_day_5_RI <- data.frame("Sample"=rep(unique(Pacli_day_5$Sample),each=10),
                             "Conc."=rep(unique(Pacli_day_5$Conc.),13),
                             "Mean_relative_intensity"=rep(0,130),
                             "Ploidy"=c(rep("1n",40),rep("2n",40), rep("3n",50)))
Pacli_day_5$Conc.[which(Pacli_day_5$Conc.==0)] <- 0.0125
Pacli_day_5_RI$Conc.[which(Pacli_day_5_RI$Conc.==0)] <- 0.0125
for (i in unique(Pacli_day_5$Sample)){
  for (j in unique(Pacli_day_5$Conc.)){
    Pacli_day_5_RI[Pacli_day_5_RI$Sample==i & Pacli_day_5_RI$Conc.==j,3] <- mean(Pacli_day_5[Pacli_day_5$Sample==i & Pacli_day_5$Conc.==j,4])
  }
}
Pacli_day_5_RI <- Pacli_day_5_RI[Pacli_day_5_RI$Sample!="2n 10R p40",]

stat <- Pacli_day_5_RI[Pacli_day_5_RI$Conc.>0.02,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(Pacli_day_5_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=1.3) +
  stat_summary(fun=mean,geom='point') + ggtitle("Paclitaxel day 5")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.3,1.35),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(Pacli_day_5_RI$Conc.))[2:10]),breaks = unique(as.numeric(Pacli_day_5_RI$Conc.)))
# ggsave("WT- Pacli day 5 by ploidy.pdf",device = "pdf")

### calculate and compare IC50 and AUC between ploidies ###
library(bayestestR)
formula = y ~ poly(x, 2, raw = TRUE)

Pacli_D2=ggplot(Pacli_day_2_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","green","forestgreen","darkgreen","deepskyblue","dodgerblue","darkblue","lightcoral","tomato","red","firebrick","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 6,formula = formula)+
  ggtitle("Paclitaxel - day 2")+scale_x_continuous(limits = c(0,10))
build=ggplot_build(Pacli_D2)
Pacli_D2+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")
eq2=str_remove_all(str_remove_all(str_remove_all(str_remove_all(build[["data"]][[3]][["eq.label"]],pattern = "italic\\("),"\\)"),"~"),"`")
match_names <- c("1n 10R p43","1n 10G p43","1n 10G p42","1n 10R p40",
                 "2n 10G p45","2n 10R p48","2n 10G p47",
                 "3n clone H f+8","3n clone E f+6","3n clone C f+7","3n clone B f+7","3n clone J f+7")
names(eq2) <- match_names
coeff <- str_replace_all(str_replace_all(str_remove(str_remove(str_remove(str_remove_all(str_replace(str_remove(string = eq2,"y="),pattern = "- ",replacement = "-"),"\\*x"),"\\+ "),"\\*"),"\\^2"),"%%","*"),"- ","-")
coeff<- as.numeric(str_replace_all(str_replace_all(str_remove_all(unlist(strsplit(coeff," ")),"\\*"),pattern = "10\\^",replacement = "e"),pattern = "x","1"))
ic50D2 <- c()
for (i in c(1,4,7,10,13,16,19,22,25,28,31,34)){
  ic50D2 <- c(ic50D2,Re(polyroot(c((coeff[i]-0.5),coeff[i+1],coeff[i+2]))))
}
ic50D2 <- data.frame("Sample"=match_names,"IC50.1"=ic50D2[c(1,3,5,7,9,11,13,15,17,19,21,23)],"IC50.2"=ic50D2[c(2,4,6,8,10,12,14,16,18,20,22,24)])

aucD2 <- c()
for (i in unique(Pacli_day_2_RI$Sample)){
  aucD2 <- c(aucD2,auc(Pacli_day_2_RI[Pacli_day_2_RI$Sample==i,2],Pacli_day_2_RI[Pacli_day_2_RI$Sample==i,3]))
}
names(aucD2) <- unique(Pacli_day_2_RI$Sample)


Pacli_D3=ggplot(Pacli_day_3_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","green","forestgreen","darkgreen","deepskyblue","dodgerblue","darkblue","lightcoral","tomato","red","firebrick","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 6,formula = formula)+
  ggtitle("Paclitaxel - day 3")+scale_x_continuous(limits = c(0,10))
build=ggplot_build(Pacli_D3)
Pacli_D3+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")
eq3=str_remove_all(str_remove_all(str_remove_all(str_remove_all(build[["data"]][[3]][["eq.label"]],pattern = "italic\\("),"\\)"),"~"),"`")
match_names <- c("1n 10R p43","1n 10G p43","1n 10G p42","1n 10R p40",
                 "2n 10G p45","2n 10R p48","2n 10G p47",
                 "3n clone H f+8","3n clone E f+6","3n clone C f+7","3n clone B f+7","3n clone J f+7")
names(eq3) <- match_names
coeff <- str_replace_all(str_replace_all(str_remove(str_remove(str_remove(str_remove_all(str_replace(str_remove(string = eq3,"y="),pattern = "- ",replacement = "-"),"\\*x"),"\\+ "),"\\*"),"\\^2"),"%%","*"),"- ","-")
coeff<- as.numeric(str_replace_all(str_remove_all(unlist(strsplit(coeff," ")),"\\*"),pattern = "10\\^",replacement = "e"))
ic50D3 <- c()
for (i in c(1,4,7,10,13,16,19,22,25,28,31,34)){
  ic50D3 <- c(ic50D3,Re(polyroot(c((coeff[i]-0.5),coeff[i+1],coeff[i+2]))))
}
ic50D3 <- data.frame("Sample"=match_names,"IC50.1"=ic50D3[c(1,3,5,7,9,11,13,15,17,19,21,23)],"IC50.2"=ic50D3[c(2,4,6,8,10,12,14,16,18,20,22,24)])

aucD3 <- c()
for (i in unique(Pacli_day_3_RI$Sample)){
  aucD3 <- c(aucD3,auc(Pacli_day_3_RI[Pacli_day_3_RI$Sample==i,2],Pacli_day_3_RI[Pacli_day_3_RI$Sample==i,3]))
}
names(aucD3) <- unique(Pacli_day_3_RI$Sample)

Pacli_D5=ggplot(Pacli_day_5_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","green","forestgreen","darkgreen","deepskyblue","dodgerblue","darkblue","lightcoral","tomato","red","firebrick","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 6,formula = formula)+
  ggtitle("Paclitaxel - day 5")+scale_x_continuous(limits = c(0,10))
build=ggplot_build(Pacli_D5)
Pacli_D5+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")
eq5=str_remove_all(str_remove_all(str_remove_all(str_remove_all(build[["data"]][[3]][["eq.label"]],pattern = "italic\\("),"\\)"),"~"),"`")
match_names <- c("1n 10R p43","1n 10G p43","1n 10G p42","1n 10R p40",
                 "2n 10G p45","2n 10R p48","2n 10G p47",
                 "3n clone H f+8","3n clone E f+6","3n clone C f+7","3n clone B f+7","3n clone J f+7")
names(eq5) <- match_names
coeff <- str_replace_all(str_replace_all(str_remove(str_remove(str_remove(str_remove_all(str_replace(str_remove(string = eq5,"y="),pattern = "- ",replacement = "-"),"\\*x"),"\\+ "),"\\*"),"\\^2"),"%%","*"),"- ","-")
coeff<- as.numeric(str_replace_all(str_remove_all(unlist(strsplit(coeff," ")),"\\*"),pattern = "10\\^",replacement = "e"))
ic50D5 <- c()
for (i in c(1,4,7,10,13,16,19,22,25,28,31,34)){
  ic50D5 <- c(ic50D5,Re(polyroot(c((coeff[i]-0.5),coeff[i+1],coeff[i+2]))))
}
ic50D5 <- data.frame("Sample"=match_names,"IC50.1"=ic50D5[c(1,3,5,7,9,11,13,15,17,19,21,23)],"IC50.2"=ic50D5[c(2,4,6,8,10,12,14,16,18,20,22,24)])

aucD5 <- c()
for (i in unique(Pacli_day_5_RI$Sample)){
  aucD5 <- c(aucD5,auc(Pacli_day_5_RI[Pacli_day_5_RI$Sample==i,2],Pacli_day_5_RI[Pacli_day_5_RI$Sample==i,3]))
}
names(aucD5) <- unique(Pacli_day_5_RI$Sample)


Pacli_IC50 <- rbind(ic50D2,ic50D3,ic50D5)
IC50 <- c()
for (i in 1:36){
  if (Pacli_IC50[i,2]>0){
    IC50[i] <- Pacli_IC50[i,2]
  }
  else{
    IC50[i] <- Pacli_IC50[i,3]
  }
  
}
Pacli_IC50 <- data.frame("Sample"=Pacli_IC50$Sample,IC50)
ploidy <- rep(0,36)
ploidy[grepl("1n",Pacli_IC50$Sample)] <- "1n"
ploidy[grepl("2n",Pacli_IC50$Sample)] <- "2n"
ploidy[grepl("3n",Pacli_IC50$Sample)] <- "3n"
Timepoint <- rep(c("Day 2","Day 3","Day 5"),each=12)
Pacli_IC50 <- data.frame(Pacli_IC50,"Ploidy"=ploidy,"Timepoint"=Timepoint)
colnames(Pacli_IC50)[2] <- "IC50"
write.csv(Pacli_IC50,"Pacli_IC50 - WT.csv")
Pacli_AUC <- data.frame("Sample"=c(names(aucD2),names(aucD3),names(aucD5)),"AUC"=c(aucD2,aucD3,aucD5),"Ploidy"=c("1n","1n","1n","1n","2n","2n","2n","3n","3n","3n","3n","3n"),"Timepoint"=rep(c("Day 2","Day 3","Day 5"),each=12))
write.csv(Pacli_AUC,"Pacli_AUC - WT.csv")


#mean IC50 (day 3):
stat <- Pacli_IC50[Pacli_IC50$Timepoint=="Day 3",] %>%
  t_test(IC50 ~ Ploidy, alternative = "greater", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Ploidy", dodge = 0.8)
stat$p.adj <- format(stat$p.adj,scientific=T)
ggbarplot(Pacli_IC50[Pacli_IC50$Timepoint=="Day 3",], x= "Ploidy", y= "IC50",fill="Ploidy", position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat,label = "p.adj",y.position = c(6,6.3,5), size = 3,tip.length = 0.01, hide.ns = F)+
  scale_fill_manual(values=c("forestgreen","royalblue1","red1"))+ggtitle("Paclitaxel IC50 Vs. Ploidy")+
  theme(legend.position = "none")+coord_cartesian(ylim = c(0,7),xlim = c(0.5,3.5),expand = F)
# ggsave("barplot of Paclitaxel mean IC50 Vs. ploidy of WT hESCs - day 3.pdf",device = "pdf")


#mean AUC:
stat <- Pacli_AUC %>%
  t_test(AUC ~ Ploidy, p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Ploidy", dodge = 0.8)
stat$p.adj <- format(stat$p.adj,scientific=T)
ggbarplot(Pacli_AUC, x= "Ploidy", y= "AUC",fill="Ploidy", position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat,label = "p.adj",y.position = c(6,6.5,5), size = 3,tip.length = 0.01, hide.ns = F)+
  scale_fill_manual(values=c("forestgreen","royalblue1","red1"))+ggtitle("Paclitaxel AUC Vs. Ploidy")+coord_cartesian(ylim = c(0,7),xlim = c(0.5,3.5),expand = F)
# ggsave("barplot of Paclitaxel mean AUC Vs. ploidy of WT hESCs.pdf",device = "pdf")
