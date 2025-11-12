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

# AUC calculation + Two methods for IC50 calculation: full curve linear regression vs local linear regression


# Function to calculate x at y=0.5 using full curve linear regression
calculate_x_full_regression <- function(x_vals, y_vals) {
  # Fit linear model to entire curve
  model <- lm(y_vals ~ x_vals)
  
  # Extract coefficients
  intercept <- coef(model)[1]
  slope <- coef(model)[2]
  
  # Calculate x when y = 0.5: y = intercept + slope * x
  # 0.5 = intercept + slope * x
  # x = (0.5 - intercept) / slope
  x_at_y_05 <- (0.5 - intercept) / slope
  
  return(x_at_y_05)
}

# Function to calculate x at y=0.5 using local linear regression
calculate_x_local_regression <- function(x_vals, y_vals) {
  # Find the two closest points where y=0.5 falls between them
  
  # Find points above and below y=0.5
  below_05 <- which(y_vals <= 0.5)
  above_05 <- which(y_vals >= 0.5)
  
  if (length(below_05) == 0 || length(above_05) == 0) {
    # If all points are above or below 0.5, return NA
    return(NA)
  }
  
  # Find the closest points on either side of y=0.5
  max_below_idx <- below_05[which.max(y_vals[below_05])]
  min_above_idx <- above_05[which.min(y_vals[above_05])]
  
  # Use these two points for linear regression
  local_x <- c(x_vals[max_below_idx], x_vals[min_above_idx])
  local_y <- c(y_vals[max_below_idx], y_vals[min_above_idx])
  
  # Fit linear model to these two points
  model <- lm(local_y ~ local_x)
  
  # Extract coefficients
  intercept <- coef(model)[1]
  slope <- coef(model)[2]
  
  # Calculate x when y = 0.5
  x_at_y_05 <- (0.5 - intercept) / slope
  
  return(x_at_y_05)
}

# Function to process multiple samples
calculate_IC50 <- function(data_frame) {
  # Assume data_frame has columns: sample_name, x1, x2, ..., x10, y1, y2, ..., y10
  # Or alternatively: sample_name, point_number, x_value, y_value
  
  # Get unique sample names
  sample_names <- unique(data_frame$Sample)
  
  # Initialize results data frame
  results <- data.frame(
    Sample = character(),
    x_at_y05_full_regression = numeric(),
    x_at_y05_local_regression = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process each sample
  for (sample in sample_names) {
    # Extract data for current sample
    sample_data <- data_frame[data_frame$Sample == sample, ]
    
    # Extract x and y values (assuming columns x_value and y_value exist)
    x_vals <- sample_data$Conc.
    y_vals <- sample_data$Mean_relative_intensity
    
    # Calculate x at y=0.5 using both methods
    x_full <- calculate_x_full_regression(x_vals, y_vals)
    x_local <- calculate_x_local_regression(x_vals, y_vals)
    
    # Add to results
    results <- rbind(results, data.frame(
      Sample = sample,
      x_at_y05_full_regression = x_full,
      x_at_y05_local_regression = x_local
    ))
  }
  # Create result data frame with modified name
  result_df_name <- paste0(results, "_res")
  assign(result_df_name, results, envir = .GlobalEnv)
  return(results)
}


##### Chemo VS Ploidy #####
setwd(choose.dir())
# In this analysis, we expected that the more genomic copies the p53-KO hESCs will have, 
# the higher their IC50 and AUC values will be, 
# so we use the 'alternative = "less"' argument to specify we want to check 
# whether the haploids' IC50 and AUC is lower than the diploids IC50 and AUC,
# and whether the diploids' IC50 and AUC is lower than the triploids IC50 and AUC


## Bleo ##
Bleo <- read.csv("Bleomycin_P53KO_summary (18.02.24).csv")
Bleo_relative <- na.omit(Bleo[grepl(Bleo$Sample,pattern = "Nutlin-"),])
for (i in unique(Bleo_relative$Sample)){
  for (j in colnames(Bleo_relative)[2:4]){
    Bleo_relative[Bleo_relative$Sample==i,j] <- Bleo_relative[Bleo_relative$Sample==i,j]/
      mean(Bleo_relative[Bleo_relative$Sample==i,j][1:3])
  }
}
Bleo_relative <- data.frame("Sample"=rep(Bleo_relative$Sample,3),
                            "Conc."=rep(Bleo_relative$Conc.,3),
                            "timepoint"=rep(c("Day 2","Day 3","Day 5"),each=270),
                            "Relative_intensity"=c(Bleo_relative$Day.2,Bleo_relative$Day.3,Bleo_relative$Day.5))

Bleo_relative <- data.frame(Bleo_relative,"Ploidy"=rep(c(rep("1n",90),rep("2n",90),rep("3n",90)),3))

Bleo_day_2 <- Bleo_relative[Bleo_relative$timepoint=="Day 2",]
Bleo_day_3 <- Bleo_relative[Bleo_relative$timepoint=="Day 3",]
Bleo_day_5 <- Bleo_relative[Bleo_relative$timepoint=="Day 5",]

Bleo_day_2_RI <- data.frame("Sample"=rep(unique(Bleo_day_2$Sample),each=10),
                            "Conc."=rep(unique(Bleo_day_2$Conc.),9),
                            "Mean_relative_intensity"=rep(0,90),
                            "Ploidy"=c(rep("1n",30),rep("2n",30), rep("3n",30)))
Bleo_day_2$Conc.[which(Bleo_day_2$Conc.==0)] <- 0.15625
Bleo_day_2_RI$Conc.[which(Bleo_day_2_RI$Conc.==0)] <- 0.15625
for (i in unique(Bleo_day_2$Sample)){
  for (j in unique(Bleo_day_2$Conc.)){
    Bleo_day_2_RI[Bleo_day_2_RI$Sample==i & Bleo_day_2_RI$Conc.==j,3] <- mean(Bleo_day_2[Bleo_day_2$Sample==i & Bleo_day_2$Conc.==j,4])
  }
}
shapiro_pvals <- c()
for (i in unique(Bleo_day_2_RI$Conc.)[-1]){
  for (j in unique(Bleo_day_2_RI$Ploidy)){
    shapiro_pvals <- c(shapiro_pvals,
                       shapiro.test(Bleo_day_2_RI[Bleo_day_2_RI$Conc.==i & Bleo_day_2_RI$Ploidy==j,3])$p.value)
  }
}

stat <- Bleo_day_2_RI[Bleo_day_2_RI$Conc.>0.2,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", alternative = "less", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(Bleo_day_2_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=0.75) +
  stat_summary(fun=mean,geom='point') + ggtitle("Bleomycin day 2")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.2,1.25),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(Bleo_day_2_RI$Conc.))[2:10]),breaks = unique(as.numeric(Bleo_day_2_RI$Conc.)))
# ggsave("Bleo - P53KO day 2 by ploidy.pdf",device = "pdf")

Bleo_day_3_RI <- data.frame("Sample"=rep(unique(Bleo_day_3$Sample),each=10),
                            "Conc."=rep(unique(Bleo_day_3$Conc.),9),
                            "Mean_relative_intensity"=rep(0,90),
                            "Ploidy"=c(rep("1n",30),rep("2n",30), rep("3n",30)))
Bleo_day_3$Conc.[which(Bleo_day_3$Conc.==0)] <- 0.15625
Bleo_day_3_RI$Conc.[which(Bleo_day_3_RI$Conc.==0)] <- 0.15625
for (i in unique(Bleo_day_3$Sample)){
  for (j in unique(Bleo_day_3$Conc.)){
    Bleo_day_3_RI[Bleo_day_3_RI$Sample==i & Bleo_day_3_RI$Conc.==j,3] <- mean(Bleo_day_3[Bleo_day_3$Sample==i & Bleo_day_3$Conc.==j,4])
  }
}
shapiro_pvals <- c()
for (i in unique(Bleo_day_3_RI$Conc.)[-1]){
  for (j in unique(Bleo_day_3_RI$Ploidy)){
    shapiro_pvals <- c(shapiro_pvals,
                       shapiro.test(Bleo_day_3_RI[Bleo_day_3_RI$Conc.==i & Bleo_day_3_RI$Ploidy==j,3])$p.value)
  }
}
stat <- Bleo_day_3_RI[Bleo_day_3_RI$Conc.>0.2,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", alternative = "less", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(Bleo_day_3_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=0.75) +
  stat_summary(fun=mean,geom='point') + ggtitle("Bleomycin day 3")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.1,1.15),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(Bleo_day_3_RI$Conc.))[2:10]),breaks = unique(as.numeric(Bleo_day_3_RI$Conc.)))
# ggsave("Bleo - P53KO day 3 by ploidy.pdf",device = "pdf")


Bleo_day_5_RI <- data.frame("Sample"=rep(unique(Bleo_day_5$Sample),each=10),
                            "Conc."=rep(unique(Bleo_day_5$Conc.),9),
                            "Mean_relative_intensity"=rep(0,90),
                            "Ploidy"=c(rep("1n",30),rep("2n",30), rep("3n",30)))
Bleo_day_5$Conc.[which(Bleo_day_5$Conc.==0)] <- 0.15625
Bleo_day_5_RI$Conc.[which(Bleo_day_5_RI$Conc.==0)] <- 0.15625
for (i in unique(Bleo_day_5$Sample)){
  for (j in unique(Bleo_day_5$Conc.)){
    Bleo_day_5_RI[Bleo_day_5_RI$Sample==i & Bleo_day_5_RI$Conc.==j,3] <- mean(Bleo_day_5[Bleo_day_5$Sample==i & Bleo_day_5$Conc.==j,4])
  }
}
shapiro_pvals <- c()
for (i in unique(Bleo_day_5_RI$Conc.)[-1]){
  for (j in unique(Bleo_day_5_RI$Ploidy)){
    shapiro_pvals <- c(shapiro_pvals,
                       shapiro.test(Bleo_day_5_RI[Bleo_day_5_RI$Conc.==i & Bleo_day_5_RI$Ploidy==j,3])$p.value)
  }
}
stat <- Bleo_day_5_RI[Bleo_day_5_RI$Conc.>0.2,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", alternative = "less", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(Bleo_day_5_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=0.75) +
  stat_summary(fun=mean,geom='point') + ggtitle("Bleomycin day 5")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.05,1.1),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(Bleo_day_5_RI$Conc.))[2:10]),breaks = unique(as.numeric(Bleo_day_5_RI$Conc.)))
# ggsave("Bleo - P53KO day 5 by ploidy.pdf",device = "pdf")

### calculate and compare IC50 and AUC between ploidies ###
library(bayestestR)
formula = y ~ poly(x, 2, raw = TRUE)

Bleo_D2=ggplot(Bleo_day_2_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","deepskyblue","dodgerblue","green","darkgreen","darkblue","lightcoral","red","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 50,formula = formula)+
  ggtitle("Bleomycin - day 2")+scale_x_continuous(limits = c(0,120))
build=ggplot_build(Bleo_D2)
Bleo_D2+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")

aucD2 <- c()
for (i in unique(Bleo_day_2_RI$Sample)){
  aucD2 <- c(aucD2,auc(Bleo_day_2_RI[Bleo_day_2_RI$Sample==i,2],Bleo_day_2_RI[Bleo_day_2_RI$Sample==i,3]))
}
names(aucD2) <- unique(Bleo_day_2_RI$Sample)


Bleo_D3=ggplot(Bleo_day_3_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","deepskyblue","dodgerblue","green","darkgreen","darkblue","lightcoral","red","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 50,formula = formula)+
  ggtitle("Bleomycin - day 3")+scale_x_continuous(limits = c(0,120))
build=ggplot_build(Bleo_D3)
Bleo_D3+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")

aucD3 <- c()
for (i in unique(Bleo_day_3_RI$Sample)){
  aucD3 <- c(aucD3,auc(Bleo_day_3_RI[Bleo_day_3_RI$Sample==i,2],Bleo_day_3_RI[Bleo_day_3_RI$Sample==i,3]))
}
names(aucD3) <- unique(Bleo_day_3_RI$Sample)

Bleo_D5=ggplot(Bleo_day_5_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","deepskyblue","dodgerblue","green","darkgreen","darkblue","lightcoral","red","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 50,formula = formula)+
  ggtitle("Bleomycin - day 5")+scale_x_continuous(limits = c(0,120))
build=ggplot_build(Bleo_D5)
Bleo_D5+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")

aucD5 <- c()
for (i in unique(Bleo_day_5_RI$Sample)){
  aucD5 <- c(aucD5,auc(Bleo_day_5_RI[Bleo_day_5_RI$Sample==i,2],Bleo_day_5_RI[Bleo_day_5_RI$Sample==i,3]))
}
names(aucD5) <- unique(Bleo_day_5_RI$Sample)

Bleo_AUC <- data.frame("Sample"=c(names(aucD2),names(aucD3),names(aucD5)),"AUC"=c(aucD2,aucD3,aucD5),"Ploidy"=c("1n","1n","1n","2n","2n","2n","3n","3n","3n"),"Timepoint"=rep(c("Day 2","Day 3","Day 5"),each=9))
write.csv(Bleo_AUC,"Bleo_AUC - P53KO.csv")

#Calculate IC50:
Bleo_day_2_RI_res <- calculate_IC50(Bleo_day_2_RI)
Bleo_day_3_RI_res <- calculate_IC50(Bleo_day_3_RI)
Bleo_day_5_RI_res <- calculate_IC50(Bleo_day_5_RI)
Bleo_IC50 <- rbind(Bleo_day_2_RI_res,Bleo_day_3_RI_res,Bleo_day_5_RI_res)
ploidy <- rep(rep(c("1n","2n","3n"),each=3),3)
Timepoint <- rep(c("Day 2","Day 3","Day 5"),each=9)
Bleo_IC50 <- data.frame(Bleo_IC50,"Ploidy"=ploidy,"Timepoint"=Timepoint)
colnames(Bleo_IC50)[2:3] <- c("IC50_full_regression","IC50_local_regression")


#mean IC50 (day 3):
stat <- Bleo_IC50[Bleo_IC50$Timepoint=="Day 3",] %>%
  t_test(IC50_local_regression ~ Ploidy, alternative = "less", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Ploidy", dodge = 0.8)
stat$p.adj <- format(stat$p.adj,scientific=T)
ggbarplot(Bleo_IC50[Bleo_IC50$Timepoint=="Day 3",], x= "Ploidy", y= "IC50_local_regression",fill="Ploidy", position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat,label = "p.adj",y.position = c(20,38,35), size = 3,tip.length = 0.01, hide.ns = F)+
  scale_fill_manual(values=c("forestgreen","royalblue1","red1"))+ggtitle("Bleomycin IC50 Vs. Ploidy")+
  theme(legend.position = "none")
# ggsave("barplot of Bleomycin mean IC50 Vs. ploidy of P53KO hESCs.pdf",device = "pdf")


#mean AUC (all days):
stat <- Bleo_AUC %>%
  t_test(AUC ~ Ploidy, alternative = "less", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Ploidy", dodge = 0.8)
stat$p.adj <- format(stat$p.adj,scientific=T)
ggbarplot(Bleo_AUC, x= "Ploidy", y= "AUC",fill="Ploidy", position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat,label = "p.adj",y.position = c(24,38,35), size = 3,tip.length = 0.01, hide.ns = F)+
  scale_fill_manual(values=c("forestgreen","royalblue1","red1"))+ggtitle("Bleomycin AUC Vs. Ploidy")+
  theme(legend.position = "none")
# ggsave("barplot of Bleomycin mean combined AUC Vs. ploidy of P53KO hESCs.pdf",device = "pdf")



## 5-FU ##
FU <- read.csv("5-FU_P53KO_summary (18.02.24).csv")
FU_relative <- FU[grepl(FU$Sample,pattern = "Nutlin-"),]
for (i in unique(FU_relative$Sample)){
  for (j in colnames(FU_relative)[2:4]){
    FU_relative[FU_relative$Sample==i,j] <- FU_relative[FU_relative$Sample==i,j]/
      mean(FU_relative[FU_relative$Sample==i,j][1:3],na.rm=TRUE)
  }
}
FU_relative <- data.frame("Sample"=rep(FU_relative$Sample,3),
                          "Conc."=rep(FU_relative$Conc.,3),
                          "timepoint"=rep(c("Day 2","Day 3","Day 5"),each=270),
                          "Relative_intensity"=c(FU_relative$Day.2,FU_relative$Day.3,FU_relative$Day.5))

FU_relative <- data.frame(FU_relative,"Ploidy"=rep(c(rep("1n",90),rep("2n",90),rep("3n",90)),3))

FU_day_2 <- FU_relative[FU_relative$timepoint=="Day 2",]
FU_day_3 <- FU_relative[FU_relative$timepoint=="Day 3",]
FU_day_5 <- FU_relative[FU_relative$timepoint=="Day 5",]

FU_day_2_RI <- data.frame("Sample"=rep(unique(FU_day_2$Sample),each=10),
                          "Conc."=rep(unique(FU_day_2$Conc.),9),
                          "Mean_relative_intensity"=rep(0,90),
                          "Ploidy"=c(rep("1n",30),rep("2n",30), rep("3n",30)))
FU_day_2$Conc.[which(FU_day_2$Conc.==0)] <- 3.125
FU_day_2_RI$Conc.[which(FU_day_2_RI$Conc.==0)] <- 3.125
for (i in unique(FU_day_2$Sample)){
  for (j in unique(FU_day_2$Conc.)){
    FU_day_2_RI[FU_day_2_RI$Sample==i & FU_day_2_RI$Conc.==j,3] <- mean(FU_day_2[FU_day_2$Sample==i & FU_day_2$Conc.==j,4])
  }
}
shapiro_pvals <- c()
for (i in unique(FU_day_2_RI$Conc.)[-1]){
  for (j in unique(FU_day_2_RI$Ploidy)){
    shapiro_pvals <- c(shapiro_pvals,
                       shapiro.test(FU_day_2_RI[FU_day_2_RI$Conc.==i & FU_day_2_RI$Ploidy==j,3])$p.value)
  }
}

stat <- FU_day_2_RI[FU_day_2_RI$Conc.>4,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", alternative = "less", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(FU_day_2_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=0.75) +
  stat_summary(fun=mean,geom='point') + ggtitle("5-FU day 2")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.22,1.27),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(FU_day_2_RI$Conc.))[2:10]),breaks = unique(as.numeric(FU_day_2_RI$Conc.)))
# ggsave("FU - P53KO day 2 by ploidy.pdf",device = "pdf")

FU_day_3_RI <- data.frame("Sample"=rep(unique(FU_day_3$Sample),each=10),
                          "Conc."=rep(unique(FU_day_3$Conc.),9),
                          "Mean_relative_intensity"=rep(0,90),
                          "Ploidy"=c(rep("1n",30),rep("2n",30), rep("3n",30)))
FU_day_3$Conc.[which(FU_day_3$Conc.==0)] <- 3.125
FU_day_3_RI$Conc.[which(FU_day_3_RI$Conc.==0)] <- 3.125
for (i in unique(FU_day_3$Sample)){
  for (j in unique(FU_day_3$Conc.)){
    FU_day_3_RI[FU_day_3_RI$Sample==i & FU_day_3_RI$Conc.==j,3] <- mean(FU_day_3[FU_day_3$Sample==i & FU_day_3$Conc.==j,4])
  }
}
shapiro_pvals <- c()
for (i in unique(FU_day_3_RI$Conc.)[-1]){
  for (j in unique(FU_day_3_RI$Ploidy)){
    shapiro_pvals <- c(shapiro_pvals,
                       shapiro.test(FU_day_3_RI[FU_day_3_RI$Conc.==i & FU_day_3_RI$Ploidy==j,3])$p.value)
  }
}
stat <- FU_day_3_RI[FU_day_3_RI$Conc.>4,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", alternative = "less", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(FU_day_3_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=0.75) +
  stat_summary(fun=mean,geom='point') + ggtitle("5-FU day 3")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.28,1.33),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(FU_day_3_RI$Conc.))[2:10]),breaks = unique(as.numeric(FU_day_3_RI$Conc.)))
# ggsave("FU - P53KO day 3 by ploidy.pdf",device = "pdf")


FU_day_5_RI <- data.frame("Sample"=rep(unique(FU_day_5$Sample),each=10),
                          "Conc."=rep(unique(FU_day_5$Conc.),9),
                          "Mean_relative_intensity"=rep(0,90),
                          "Ploidy"=c(rep("1n",30),rep("2n",30), rep("3n",30)))
FU_day_5$Conc.[which(FU_day_5$Conc.==0)] <- 3.125
FU_day_5_RI$Conc.[which(FU_day_5_RI$Conc.==0)] <- 3.125
for (i in unique(FU_day_5$Sample)){
  for (j in unique(FU_day_5$Conc.)){
    FU_day_5_RI[FU_day_5_RI$Sample==i & FU_day_5_RI$Conc.==j,3] <- mean(FU_day_5[FU_day_5$Sample==i & FU_day_5$Conc.==j,4])
  }
}
shapiro_pvals <- c()
for (i in unique(FU_day_5_RI$Conc.)[-1]){
  for (j in unique(FU_day_5_RI$Ploidy)){
    shapiro_pvals <- c(shapiro_pvals,
                       shapiro.test(FU_day_5_RI[FU_day_5_RI$Conc.==i & FU_day_5_RI$Ploidy==j,3])$p.value)
  }
}
stat <- FU_day_5_RI[FU_day_5_RI$Conc.>4,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", alternative = "less", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(FU_day_5_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=0.75) +
  stat_summary(fun=mean,geom='point') + ggtitle("5-FU day 5")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.5,1.55),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(FU_day_5_RI$Conc.))[2:10]),breaks = unique(as.numeric(FU_day_5_RI$Conc.)))
# ggsave("FU - P53KO day 5 by ploidy.pdf",device = "pdf")

### calculate and compare IC50 and AUC between ploidies ###
library(bayestestR)
formula = y ~ poly(x, 2, raw = TRUE)

FU_D2=ggplot(FU_day_2_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","deepskyblue","dodgerblue","green","darkgreen","darkblue","lightcoral","red","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 1600,formula = formula)+
  ggtitle("5-FU - day 2")+scale_x_continuous(limits = c(0,4000))
build=ggplot_build(FU_D2)
FU_D2+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")

aucD2 <- c()
for (i in unique(FU_day_2_RI$Sample)){
  aucD2 <- c(aucD2,auc(FU_day_2_RI[FU_day_2_RI$Sample==i,2],FU_day_2_RI[FU_day_2_RI$Sample==i,3]))
}
names(aucD2) <- unique(FU_day_2_RI$Sample)


FU_D3=ggplot(FU_day_3_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","deepskyblue","dodgerblue","green","darkgreen","darkblue","lightcoral","red","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 1600,formula = formula)+
  ggtitle("5-FU - day 3")+scale_x_continuous(limits = c(0,4000))
build=ggplot_build(FU_D3)
FU_D3+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")

aucD3 <- c()
for (i in unique(FU_day_3_RI$Sample)){
  aucD3 <- c(aucD3,auc(FU_day_3_RI[FU_day_3_RI$Sample==i,2],FU_day_3_RI[FU_day_3_RI$Sample==i,3]))
}
names(aucD3) <- unique(FU_day_3_RI$Sample)

FU_D5=ggplot(FU_day_5_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","deepskyblue","dodgerblue","green","darkgreen","darkblue","lightcoral","red","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 1600,formula = formula)+
  ggtitle("5-FU - day 5")+scale_x_continuous(limits = c(0,4000))
build=ggplot_build(FU_D5)
FU_D5+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")

aucD5 <- c()
for (i in unique(FU_day_5_RI$Sample)){
  aucD5 <- c(aucD5,auc(FU_day_5_RI[FU_day_5_RI$Sample==i,2],FU_day_5_RI[FU_day_5_RI$Sample==i,3]))
}
names(aucD5) <- unique(FU_day_5_RI$Sample)

FU_AUC <- data.frame("Sample"=c(names(aucD2),names(aucD3),names(aucD5)),"AUC"=c(aucD2,aucD3,aucD5),"Ploidy"=c("1n","1n","1n","2n","2n","2n","3n","3n","3n"),"Timepoint"=rep(c("Day 2","Day 3","Day 5"),each=9))
write.csv(FU_AUC,"FU_AUC - P53KO.csv")

#Calculate IC50:
FU_day_2_RI_res <- calculate_IC50(FU_day_2_RI)
FU_day_3_RI_res <- calculate_IC50(FU_day_3_RI)
FU_day_5_RI_res <- calculate_IC50(FU_day_5_RI)
FU_IC50 <- rbind(FU_day_2_RI_res,FU_day_3_RI_res,FU_day_5_RI_res)
ploidy <- rep(rep(c("1n","2n","3n"),each=3),3)
Timepoint <- rep(c("Day 2","Day 3","Day 5"),each=9)
FU_IC50 <- data.frame(FU_IC50,"Ploidy"=ploidy,"Timepoint"=Timepoint)
colnames(FU_IC50)[2:3] <- c("IC50_full_regression","IC50_local_regression")

# assigning full regression IC50 values where local IC50 values cannot be calculated:
FU_IC50[is.na(FU_IC50$IC50_local_regression),3] <- FU_IC50[is.na(FU_IC50$IC50_local_regression),2]


#mean IC50 (day 5):
stat <- FU_IC50[FU_IC50$Timepoint=="Day 5",] %>%
  t_test(IC50_local_regression ~ Ploidy, alternative = "less", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Ploidy", dodge = 0.8)
stat$p.adj <- format(stat$p.adj,scientific=T)
ggbarplot(FU_IC50[FU_IC50$Timepoint=="Day 5",], x= "Ploidy", y= "IC50_local_regression",fill="Ploidy", position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat,label = "p.adj",y.position = c(325,800,700), size = 3,tip.length = 0.01, hide.ns = F)+
  scale_fill_manual(values=c("forestgreen","royalblue1","red1"))+ggtitle("5-FU IC50 Vs. Ploidy")+
  theme(legend.position = "none")
# ggsave("barplot of 5-FU mean IC50 Vs. ploidy of P53KO hESCs.pdf",device = "pdf")


#mean AUC (all days):
stat <- FU_AUC %>%
  t_test(AUC ~ Ploidy, alternative = "less", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Ploidy", dodge = 0.8)
stat$p.adj <- format(stat$p.adj,scientific=T)
ggbarplot(FU_AUC, x= "Ploidy", y= "AUC",fill="Ploidy", position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat,label = "p.adj",y.position = c(1050,1350,1250), size = 3,tip.length = 0.01, hide.ns = F)+
  scale_fill_manual(values=c("forestgreen","royalblue1","red1"))+ggtitle("5-FU AUC Vs. Ploidy")+
  theme(legend.position = "none")
# ggsave("barplot of 5-FU mean combined AUC Vs. ploidy of P53KO hESCs.pdf",device = "pdf")



## Cis ##
Cis <- read.csv("Cisplatin_P53KO_summary (18.02.24).csv")
Cis_relative <- na.omit(Cis[grepl(Cis$Sample,pattern = "Nutlin-"),])
for (i in unique(Cis_relative$Sample)){
  for (j in colnames(Cis_relative)[2:4]){
    Cis_relative[Cis_relative$Sample==i,j] <- Cis_relative[Cis_relative$Sample==i,j]/
      mean(Cis_relative[Cis_relative$Sample==i,j][1:3])
  }
}
Cis_relative <- data.frame("Sample"=rep(Cis_relative$Sample,3),
                           "Conc."=rep(Cis_relative$Conc.,3),
                           "timepoint"=rep(c("Day 2","Day 3","Day 5"),each=270),
                           "Relative_intensity"=c(Cis_relative$Day.2,Cis_relative$Day.3,Cis_relative$Day.5))

Cis_relative <- data.frame(Cis_relative,"Ploidy"=rep(c(rep("1n",90),rep("2n",90),rep("3n",90)),3))

Cis_day_2 <- Cis_relative[Cis_relative$timepoint=="Day 2",]
Cis_day_3 <- Cis_relative[Cis_relative$timepoint=="Day 3",]
Cis_day_5 <- Cis_relative[Cis_relative$timepoint=="Day 5",]

Cis_day_2_RI <- data.frame("Sample"=rep(unique(Cis_day_2$Sample),each=10),
                           "Conc."=rep(unique(Cis_day_2$Conc.),9),
                           "Mean_relative_intensity"=rep(0,90),
                           "Ploidy"=c(rep("1n",30),rep("2n",30), rep("3n",30)))
Cis_day_2$Conc.[which(Cis_day_2$Conc.==0)] <- 0.1953125
Cis_day_2_RI$Conc.[which(Cis_day_2_RI$Conc.==0)] <- 0.1953125
for (i in unique(Cis_day_2$Sample)){
  for (j in unique(Cis_day_2$Conc.)){
    Cis_day_2_RI[Cis_day_2_RI$Sample==i & Cis_day_2_RI$Conc.==j,3] <- mean(Cis_day_2[Cis_day_2$Sample==i & Cis_day_2$Conc.==j,4])
  }
}
shapiro_pvals <- c()
for (i in unique(Cis_day_2_RI$Conc.)[-1]){
  for (j in unique(Cis_day_2_RI$Ploidy)){
    shapiro_pvals <- c(shapiro_pvals,
                       shapiro.test(Cis_day_2_RI[Cis_day_2_RI$Conc.==i & Cis_day_2_RI$Ploidy==j,3])$p.value)
  }
}

stat <- Cis_day_2_RI[Cis_day_2_RI$Conc.>0.2,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", alternative = "less", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(Cis_day_2_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=0.75) +
  stat_summary(fun=mean,geom='point') + ggtitle("Cisplatin day 2")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.35,1.4),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(Cis_day_2_RI$Conc.))[2:10]),breaks = unique(as.numeric(Cis_day_2_RI$Conc.)))
# ggsave("Cis - P53KO day 2 by ploidy.pdf",device = "pdf")

Cis_day_3_RI <- data.frame("Sample"=rep(unique(Cis_day_3$Sample),each=10),
                           "Conc."=rep(unique(Cis_day_3$Conc.),9),
                           "Mean_relative_intensity"=rep(0,90),
                           "Ploidy"=c(rep("1n",30),rep("2n",30), rep("3n",30)))
Cis_day_3$Conc.[which(Cis_day_3$Conc.==0)] <- 0.1953125
Cis_day_3_RI$Conc.[which(Cis_day_3_RI$Conc.==0)] <- 0.1953125
for (i in unique(Cis_day_3$Sample)){
  for (j in unique(Cis_day_3$Conc.)){
    Cis_day_3_RI[Cis_day_3_RI$Sample==i & Cis_day_3_RI$Conc.==j,3] <- mean(Cis_day_3[Cis_day_3$Sample==i & Cis_day_3$Conc.==j,4])
  }
}
shapiro_pvals <- c()
for (i in unique(Cis_day_3_RI$Conc.)[-1]){
  for (j in unique(Cis_day_3_RI$Ploidy)){
    shapiro_pvals <- c(shapiro_pvals,
                       shapiro.test(Cis_day_3_RI[Cis_day_3_RI$Conc.==i & Cis_day_3_RI$Ploidy==j,3])$p.value)
  }
}
stat <- Cis_day_3_RI[Cis_day_3_RI$Conc.>0.2,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", alternative = "less", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(Cis_day_3_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=0.75) +
  stat_summary(fun=mean,geom='point') + ggtitle("Cisplatin day 3")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.35,1.4),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(Cis_day_3_RI$Conc.))[2:10]),breaks = unique(as.numeric(Cis_day_3_RI$Conc.)))
# ggsave("Cis - P53KO day 3 by ploidy.pdf",device = "pdf")


Cis_day_5_RI <- data.frame("Sample"=rep(unique(Cis_day_5$Sample),each=10),
                           "Conc."=rep(unique(Cis_day_5$Conc.),9),
                           "Mean_relative_intensity"=rep(0,90),
                           "Ploidy"=c(rep("1n",30),rep("2n",30), rep("3n",30)))
Cis_day_5$Conc.[which(Cis_day_5$Conc.==0)] <- 0.1953125
Cis_day_5_RI$Conc.[which(Cis_day_5_RI$Conc.==0)] <- 0.1953125
for (i in unique(Cis_day_5$Sample)){
  for (j in unique(Cis_day_5$Conc.)){
    Cis_day_5_RI[Cis_day_5_RI$Sample==i & Cis_day_5_RI$Conc.==j,3] <- mean(Cis_day_5[Cis_day_5$Sample==i & Cis_day_5$Conc.==j,4])
  }
}
shapiro_pvals <- c()
for (i in unique(Cis_day_5_RI$Conc.)[-1]){
  for (j in unique(Cis_day_5_RI$Ploidy)){
    shapiro_pvals <- c(shapiro_pvals,
                       shapiro.test(Cis_day_5_RI[Cis_day_5_RI$Conc.==i & Cis_day_5_RI$Ploidy==j,3])$p.value)
  }
}
stat <- Cis_day_5_RI[Cis_day_5_RI$Conc.>0.2,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", alternative = "less", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(Cis_day_5_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=0.75) +
  stat_summary(fun=mean,geom='point') + ggtitle("Cisplatin day 5")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.25,1.3),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(Cis_day_5_RI$Conc.))[2:10]),breaks = unique(as.numeric(Cis_day_5_RI$Conc.)))
# ggsave("Cis - P53KO day 5 by ploidy.pdf",device = "pdf")

### calculate and compare IC50 and AUC between ploidies ###
library(bayestestR)
formula = y ~ poly(x, 2, raw = TRUE)

Cis_D2=ggplot(Cis_day_2_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","deepskyblue","dodgerblue","green","darkgreen","darkblue","lightcoral","red","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 50,formula = formula)+
  ggtitle("Cisplatin - day 2")+scale_x_continuous(limits = c(0,120))
build=ggplot_build(Cis_D2)
Cis_D2+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")

aucD2 <- c()
for (i in unique(Cis_day_2_RI$Sample)){
  aucD2 <- c(aucD2,auc(Cis_day_2_RI[Cis_day_2_RI$Sample==i,2],Cis_day_2_RI[Cis_day_2_RI$Sample==i,3]))
}
names(aucD2) <- unique(Cis_day_2_RI$Sample)


Cis_D3=ggplot(Cis_day_3_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","deepskyblue","dodgerblue","green","darkgreen","darkblue","lightcoral","red","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 50,formula = formula)+
  ggtitle("Cisplatin - day 3")+scale_x_continuous(limits = c(0,120))
build=ggplot_build(Cis_D3)
Cis_D3+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")

aucD3 <- c()
for (i in unique(Cis_day_3_RI$Sample)){
  aucD3 <- c(aucD3,auc(Cis_day_3_RI[Cis_day_3_RI$Sample==i,2],Cis_day_3_RI[Cis_day_3_RI$Sample==i,3]))
}
names(aucD3) <- unique(Cis_day_3_RI$Sample)

Cis_D5=ggplot(Cis_day_5_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","deepskyblue","dodgerblue","green","darkgreen","darkblue","lightcoral","red","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 50,formula = formula)+
  ggtitle("Cisplatin - day 5")+scale_x_continuous(limits = c(0,120))
build=ggplot_build(Cis_D5)
Cis_D5+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")

aucD5 <- c()
for (i in unique(Cis_day_5_RI$Sample)){
  aucD5 <- c(aucD5,auc(Cis_day_5_RI[Cis_day_5_RI$Sample==i,2],Cis_day_5_RI[Cis_day_5_RI$Sample==i,3]))
}
names(aucD5) <- unique(Cis_day_5_RI$Sample)

Cis_AUC <- data.frame("Sample"=c(names(aucD2),names(aucD3),names(aucD5)),"AUC"=c(aucD2,aucD3,aucD5),"Ploidy"=c("1n","1n","1n","2n","2n","2n","3n","3n","3n"),"Timepoint"=rep(c("Day 2","Day 3","Day 5"),each=9))
write.csv(Cis_AUC,"Cis_AUC - P53KO.csv")

#Calculate IC50:
Cis_day_2_RI_res <- calculate_IC50(Cis_day_2_RI)
Cis_day_3_RI_res <- calculate_IC50(Cis_day_3_RI)
Cis_day_5_RI_res <- calculate_IC50(Cis_day_5_RI)
Cis_IC50 <- rbind(Cis_day_2_RI_res,Cis_day_3_RI_res,Cis_day_5_RI_res)
ploidy <- rep(rep(c("1n","2n","3n"),each=3),3)
Timepoint <- rep(c("Day 2","Day 3","Day 5"),each=9)
Cis_IC50 <- data.frame(Cis_IC50,"Ploidy"=ploidy,"Timepoint"=Timepoint)
colnames(Cis_IC50)[2:3] <- c("IC50_full_regression","IC50_local_regression")


#mean IC50 (day 5):
stat <- Cis_IC50[Cis_IC50$Timepoint=="Day 5",] %>%
  t_test(IC50_local_regression ~ Ploidy, alternative = "less", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Ploidy", dodge = 0.8)
stat$p.adj <- format(stat$p.adj,scientific=T)
ggbarplot(Cis_IC50[Cis_IC50$Timepoint=="Day 5",], x= "Ploidy", y= "IC50_local_regression",fill="Ploidy", position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat,label = "p.adj",y.position = c(7,13,11), size = 3,tip.length = 0.01, hide.ns = F)+
  scale_fill_manual(values=c("forestgreen","royalblue1","red1"))+ggtitle("Cisplatin IC50 Vs. Ploidy")+
  theme(legend.position = "none")
# ggsave("barplot of Cisplatin mean IC50 Vs. ploidy of P53KO hESCs.pdf",device = "pdf")


#mean AUC (all days):
stat <- Cis_AUC %>%
  t_test(AUC ~ Ploidy, alternative = "less", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Ploidy", dodge = 0.8)
stat$p.adj <- format(stat$p.adj,scientific=T)
ggbarplot(Cis_AUC, x= "Ploidy", y= "AUC",fill="Ploidy", position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat,label = "p.adj",y.position = c(14,24,22), size = 3,tip.length = 0.01, hide.ns = F)+
  scale_fill_manual(values=c("forestgreen","royalblue1","red1"))+ggtitle("Cisplatin AUC Vs. Ploidy")+
  theme(legend.position = "none")
# ggsave("barplot of Cisplatin combined mean AUC Vs. ploidy of P53KO hESCs.pdf",device = "pdf")



## Pacli ##
Pacli <- read.csv("Paclitaxel_P53KO_summary (18.02.24).csv")
Pacli_relative <- na.omit(Pacli[grepl(Pacli$Sample,pattern = "Nutlin-"),])
for (i in unique(Pacli_relative$Sample)){
  for (j in colnames(Pacli_relative)[2:4]){
    Pacli_relative[Pacli_relative$Sample==i,j] <- Pacli_relative[Pacli_relative$Sample==i,j]/
      mean(Pacli_relative[Pacli_relative$Sample==i,j][1:3])
  }
}
Pacli_relative <- data.frame("Sample"=rep(Pacli_relative$Sample,3),
                             "Conc."=rep(Pacli_relative$Conc.,3),
                             "timepoint"=rep(c("Day 2","Day 3","Day 5"),each=270),
                             "Relative_intensity"=c(Pacli_relative$Day.2,Pacli_relative$Day.3,Pacli_relative$Day.5))

Pacli_relative <- data.frame(Pacli_relative,"Ploidy"=rep(c(rep("1n",90),rep("2n",90),rep("3n",90)),3))

Pacli_day_2 <- Pacli_relative[Pacli_relative$timepoint=="Day 2",]
Pacli_day_3 <- Pacli_relative[Pacli_relative$timepoint=="Day 3",]
Pacli_day_5 <- Pacli_relative[Pacli_relative$timepoint=="Day 5",]

Pacli_day_2_RI <- data.frame("Sample"=rep(unique(Pacli_day_2$Sample),each=10),
                             "Conc."=rep(unique(Pacli_day_2$Conc.),9),
                             "Mean_relative_intensity"=rep(0,90),
                             "Ploidy"=c(rep("1n",30),rep("2n",30), rep("3n",30)))
Pacli_day_2$Conc.[which(Pacli_day_2$Conc.==0)] <- 0.09765625
Pacli_day_2_RI$Conc.[which(Pacli_day_2_RI$Conc.==0)] <- 0.09765625
for (i in unique(Pacli_day_2$Sample)){
  for (j in unique(Pacli_day_2$Conc.)){
    Pacli_day_2_RI[Pacli_day_2_RI$Sample==i & Pacli_day_2_RI$Conc.==j,3] <- mean(Pacli_day_2[Pacli_day_2$Sample==i & Pacli_day_2$Conc.==j,4])
  }
}
shapiro_pvals <- c()
for (i in unique(Pacli_day_2_RI$Conc.)[-1]){
  for (j in unique(Pacli_day_2_RI$Ploidy)){
    shapiro_pvals <- c(shapiro_pvals,
                       shapiro.test(Pacli_day_2_RI[Pacli_day_2_RI$Conc.==i & Pacli_day_2_RI$Ploidy==j,3])$p.value)
  }
}

stat <- Pacli_day_2_RI[Pacli_day_2_RI$Conc.>0.1,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", alternative = "less", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(Pacli_day_2_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=0.75) +
  stat_summary(fun=mean,geom='point') + ggtitle("Paclitaxel day 2")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.2,1.25),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(Pacli_day_2_RI$Conc.))[2:10]),breaks = unique(as.numeric(Pacli_day_2_RI$Conc.)))
# ggsave("Pacli - P53KO day 2 by ploidy.pdf",device = "pdf")

Pacli_day_3_RI <- data.frame("Sample"=rep(unique(Pacli_day_3$Sample),each=10),
                             "Conc."=rep(unique(Pacli_day_3$Conc.),9),
                             "Mean_relative_intensity"=rep(0,90),
                             "Ploidy"=c(rep("1n",30),rep("2n",30), rep("3n",30)))
Pacli_day_3$Conc.[which(Pacli_day_3$Conc.==0)] <- 0.09765625
Pacli_day_3_RI$Conc.[which(Pacli_day_3_RI$Conc.==0)] <- 0.09765625
for (i in unique(Pacli_day_3$Sample)){
  for (j in unique(Pacli_day_3$Conc.)){
    Pacli_day_3_RI[Pacli_day_3_RI$Sample==i & Pacli_day_3_RI$Conc.==j,3] <- mean(Pacli_day_3[Pacli_day_3$Sample==i & Pacli_day_3$Conc.==j,4])
  }
}
shapiro_pvals <- c()
for (i in unique(Pacli_day_3_RI$Conc.)[-1]){
  for (j in unique(Pacli_day_3_RI$Ploidy)){
    shapiro_pvals <- c(shapiro_pvals,
                       shapiro.test(Pacli_day_3_RI[Pacli_day_3_RI$Conc.==i & Pacli_day_3_RI$Ploidy==j,3])$p.value)
  }
}
stat <- Pacli_day_3_RI[Pacli_day_3_RI$Conc.>0.1,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", alternative = "less", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(Pacli_day_3_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=0.75) +
  stat_summary(fun=mean,geom='point') + ggtitle("Paclitaxel day 3")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.1,1.15),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(Pacli_day_3_RI$Conc.))[2:10]),breaks = unique(as.numeric(Pacli_day_3_RI$Conc.)))
# ggsave("Pacli - P53KO day 3 by ploidy.pdf",device = "pdf")


Pacli_day_5_RI <- data.frame("Sample"=rep(unique(Pacli_day_5$Sample),each=10),
                             "Conc."=rep(unique(Pacli_day_5$Conc.),9),
                             "Mean_relative_intensity"=rep(0,90),
                             "Ploidy"=c(rep("1n",30),rep("2n",30), rep("3n",30)))
Pacli_day_5$Conc.[which(Pacli_day_5$Conc.==0)] <- 0.09765625
Pacli_day_5_RI$Conc.[which(Pacli_day_5_RI$Conc.==0)] <- 0.09765625
for (i in unique(Pacli_day_5$Sample)){
  for (j in unique(Pacli_day_5$Conc.)){
    Pacli_day_5_RI[Pacli_day_5_RI$Sample==i & Pacli_day_5_RI$Conc.==j,3] <- mean(Pacli_day_5[Pacli_day_5$Sample==i & Pacli_day_5$Conc.==j,4])
  }
}
shapiro_pvals <- c()
for (i in unique(Pacli_day_5_RI$Conc.)[-1]){
  for (j in unique(Pacli_day_5_RI$Ploidy)){
    shapiro_pvals <- c(shapiro_pvals,
                       shapiro.test(Pacli_day_5_RI[Pacli_day_5_RI$Conc.==i & Pacli_day_5_RI$Ploidy==j,3])$p.value)
  }
}
stat <- Pacli_day_5_RI[Pacli_day_5_RI$Conc.>0.1,] %>%
  group_by(Conc.) %>%
  t_test(Mean_relative_intensity~Ploidy,ref.group = "2n", alternative = "less", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Conc.")
stat[,14:16] <- stat[,1]

ggplot(Pacli_day_5_RI, aes(x=Conc., y=Mean_relative_intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=0.75) +
  stat_summary(fun=mean,geom='point') + ggtitle("Paclitaxel day 5")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',aes(color=Ploidy),width=0.1) +
  theme_classic()+scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),expand = c(0,0.005))+  
  scale_colour_manual(values=c("forestgreen","royalblue1","red1"))+
  # stat_pvalue_manual(data = stat,label.size = 3,remove.bracket = T,color = "group2",y.position = rep(c(1.2,1.25),9), label = "p.adj.signif", hide.ns = F)+
  scale_x_continuous(trans='log2',labels = c(0,unique(as.numeric(Pacli_day_5_RI$Conc.))[2:10]),breaks = unique(as.numeric(Pacli_day_5_RI$Conc.)))
# ggsave("Pacli - P53KO day 5 by ploidy.pdf",device = "pdf")

### calculate and compare IC50 and AUC between ploidies ###
library(bayestestR)
formula = y ~ poly(x, 2, raw = TRUE)

Pacli_D2=ggplot(Pacli_day_2_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","deepskyblue","dodgerblue","green","darkgreen","darkblue","lightcoral","red","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 50,formula = formula)+
  ggtitle("Paclitaxel - day 2")+scale_x_continuous(limits = c(0,120))
build=ggplot_build(Pacli_D2)
Pacli_D2+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")

aucD2 <- c()
for (i in unique(Pacli_day_2_RI$Sample)){
  aucD2 <- c(aucD2,auc(Pacli_day_2_RI[Pacli_day_2_RI$Sample==i,2],Pacli_day_2_RI[Pacli_day_2_RI$Sample==i,3]))
}
names(aucD2) <- unique(Pacli_day_2_RI$Sample)


Pacli_D3=ggplot(Pacli_day_3_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","deepskyblue","dodgerblue","green","darkgreen","darkblue","lightcoral","red","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 50,formula = formula)+
  ggtitle("Paclitaxel - day 3")+scale_x_continuous(limits = c(0,120))
build=ggplot_build(Pacli_D3)
Pacli_D3+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")

aucD3 <- c()
for (i in unique(Pacli_day_3_RI$Sample)){
  aucD3 <- c(aucD3,auc(Pacli_day_3_RI[Pacli_day_3_RI$Sample==i,2],Pacli_day_3_RI[Pacli_day_3_RI$Sample==i,3]))
}
names(aucD3) <- unique(Pacli_day_3_RI$Sample)

Pacli_D5=ggplot(Pacli_day_5_RI, aes(x=Conc., y=Mean_relative_intensity)) +geom_point(aes(fill = Sample, color = Sample))+
  theme_bw()+stat_smooth(aes(fill = Sample, color = Sample), method = "lm",se = F, formula = formula)+
  scale_colour_manual(values=c("lightgreen","deepskyblue","dodgerblue","green","darkgreen","darkblue","lightcoral","red","darkred"))+
  stat_regline_equation(aes(color = Sample,label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),label.y.npc = 1,label.x = 50,formula = formula)+
  ggtitle("Paclitaxel - day 5")+scale_x_continuous(limits = c(0,120))
build=ggplot_build(Pacli_D5)
Pacli_D5+geom_hline(yintercept = 0.5,lty=2,color="darkgrey")

aucD5 <- c()
for (i in unique(Pacli_day_5_RI$Sample)){
  aucD5 <- c(aucD5,auc(Pacli_day_5_RI[Pacli_day_5_RI$Sample==i,2],Pacli_day_5_RI[Pacli_day_5_RI$Sample==i,3]))
}
names(aucD5) <- unique(Pacli_day_5_RI$Sample)

Pacli_AUC <- data.frame("Sample"=c(names(aucD2),names(aucD3),names(aucD5)),"AUC"=c(aucD2,aucD3,aucD5),"Ploidy"=c("1n","1n","1n","2n","2n","2n","3n","3n","3n"),"Timepoint"=rep(c("Day 2","Day 3","Day 5"),each=9))
write.csv(Pacli_AUC,"Pacli_AUC - P53KO.csv")

#Calculate IC50:
Pacli_day_2_RI_res <- calculate_IC50(Pacli_day_2_RI)
Pacli_day_3_RI_res <- calculate_IC50(Pacli_day_3_RI)
Pacli_day_5_RI_res <- calculate_IC50(Pacli_day_5_RI)
Pacli_IC50 <- rbind(Pacli_day_2_RI_res,Pacli_day_3_RI_res,Pacli_day_5_RI_res)
ploidy <- rep(rep(c("1n","2n","3n"),each=3),3)
Timepoint <- rep(c("Day 2","Day 3","Day 5"),each=9)
Pacli_IC50 <- data.frame(Pacli_IC50,"Ploidy"=ploidy,"Timepoint"=Timepoint)
colnames(Pacli_IC50)[2:3] <- c("IC50_full_regression","IC50_local_regression")


#mean IC50 (day 2):
stat <- Pacli_IC50[Pacli_IC50$Timepoint=="Day 2",] %>%
  t_test(IC50_local_regression ~ Ploidy, alternative = "less", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Ploidy", dodge = 0.8)
stat$p.adj <- format(stat$p.adj,scientific=T)
ggbarplot(Pacli_IC50[Pacli_IC50$Timepoint=="Day 2",], x= "Ploidy", y= "IC50_local_regression",fill="Ploidy", position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat,label = "p.adj",y.position = c(15,24,22), size = 3,tip.length = 0.01, hide.ns = F)+
  scale_fill_manual(values=c("forestgreen","royalblue1","red1"))+ggtitle("Paclitaxel IC50 Vs. Ploidy")+
  theme(legend.position = "none")
# ggsave("barplot of Paclitaxel mean IC50 Vs. ploidy of P53KO hESCs.pdf",device = "pdf")


#mean AUC (all days):
stat <- Pacli_AUC %>%
  t_test(AUC ~ Ploidy, alternative = "less", p.adjust.method = "BH")
stat <- stat %>% add_xy_position(x = "Ploidy", dodge = 0.8)
stat$p.adj <- format(stat$p.adj,scientific=T)
ggbarplot(Pacli_AUC, x= "Ploidy", y= "AUC",fill="Ploidy", position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat,label = "p.adj",y.position = c(20,30,27), size = 3,tip.length = 0.01, hide.ns = F)+
  scale_fill_manual(values=c("forestgreen","royalblue1","red1"))+ggtitle("Paclitaxel AUC Vs. Ploidy")+
  theme(legend.position = "none")
# ggsave("barplot of Paclitaxel mean combined AUC Vs. ploidy of P53KO hESCs.pdf",device = "pdf")
