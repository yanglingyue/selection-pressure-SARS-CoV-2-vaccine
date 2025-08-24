setwd("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine")
# Load packages (assumed already installed)
rm(list = ls())
library(data.table);library(corrplot);library(ggcorrplot);library(RColorBrewer)
library(ggpubr);

dat1<-read.csv("01_data/subsample_adaptation_vaccine_low.csv")

source("load_package.R")
source("plot_theme.R")
source("load_data_clean.R")

dir.create("04_extended_analysis/plot_output/Figure_S5")
outpath <-  "04_extended_analysis/plot_output/Figure_S5/"

custom_colors <- c("#00468BFF","#ED0000FF",'#F79F1F',"#42B540B2","#7570b3",
                   "#1b9e77",'#b33939','#ffb142',"#ED0000FF","#4DBBD5FF")
                                         
####################################fig.xx##############################################
#####Scatter plot of adjusted vaccine coverage and selection pressure on S regions######
###Before the emergence of Omicron variants and during Omicron variants were dominant###
dat4 <- dat1
dat4$type <- NA
dat4[which(dat4$month=="2021-11"|dat4$month=="2021-12"|dat4$month=="2022-01"|dat4$month=="2022-02"|
             dat4$month=="2022-03"|dat4$month=="2022-04"|dat4$month=="2022-05"|dat4$month=="2022-06"|
             dat4$month=="2022-07"|dat4$month=="2022-08"|
             dat4$month=="2022-09"),"type"] <- "From  November 2021 to September 2022"

dat4[which(is.na(dat4$type)),"type"] <- "From  March 2020 to October 2021"

dat4_S1= subset(dat4,dat4$gene== "S1")

p1_dnds <- ggplot(dat4_S1) +
  geom_point(aes(x = FW_BW_P, y = dnds_mean_focal,color=type), size=2,shape=19 ,alpha=0.5) +
  xlab("Ajusted vaccine coverage") +
  ylab("Nonsynonymous to synonymous\n divergence ratio") +
  stat_smooth(aes(x=FW_BW_P,y=dnds_mean_focal,color=type,fill=type),method=lm,se=T,formula = y ~ x)+
  scale_x_continuous(limits=c(0,80),breaks = seq(0,75,25))+
  scale_y_continuous(limits=c(0,9),breaks = seq(0,9,3))+
  stat_cor(aes(x=FW_BW_P,y=dnds_mean_focal,color=type),size=6,method = "spearman",size=6)+
  ggtitle(paste0(""))+
  mytheme+
  scale_fill_manual(name = 'Time period',values = custom_colors)+
  scale_color_manual(name = 'Time period',values = custom_colors)+
  theme(plot.margin=unit(c(0.25,1.25,1.25,1.25),'lines'))

p1_dnds

p1_nonsyn <- ggplot(dat4_S1) +
  geom_point(aes(x = FW_BW_P, y = muts_nonsyn_focal,color=type), size=2,shape=19 ,alpha=0.5) +
  xlab("Ajusted vaccine coverage") +
  ylab("Number of nonsynonymous\n mutations") +
  stat_smooth(aes(x=FW_BW_P,y=muts_nonsyn_focal,color=type,fill=type),method=lm,se=T,formula = y ~ x)+
  scale_x_continuous(limits=c(0,80),breaks = seq(0,75,25))+
  scale_y_continuous(limits=c(0,40),breaks = seq(0,40,10))+
  stat_cor(aes(x=FW_BW_P,y=muts_nonsyn_focal,color=type),size=6,method = "spearman",size=6)+
  ggtitle(paste0(""))+
  mytheme+
  scale_fill_manual(name = 'Time period',values = custom_colors)+
  scale_color_manual(name = 'Time period',values = custom_colors)+
  theme(plot.margin=unit(c(0.25,1.25,1.25,1.25),'lines'))

p1_nonsyn

p1_fitness <- ggplot(dat4_S1) +
  geom_point(aes(x = FW_BW_P, y = fitness_focal,color=type), size=2,shape=19 ,alpha=0.5) +
  xlab("Ajusted vaccine coverage") +
  ylab("Mutational fitness") +
  stat_smooth(aes(x=FW_BW_P,y=fitness_focal,color=type,fill=type),method=lm,se=T,formula = y ~ x)+
  scale_x_continuous(limits=c(0,80),breaks = seq(0,75,25))+
  scale_y_continuous(limits=c(-0.1,4),breaks = seq(0,4,1))+
  stat_cor(aes(x=FW_BW_P,y=fitness_focal,color=type),size=6,method = "spearman",size=6)+
  ggtitle(paste0(""))+
  mytheme+
  scale_fill_manual(name = 'Time period',values = custom_colors)+
  scale_color_manual(name = 'Time period',values = custom_colors)+
  theme(plot.margin=unit(c(0.25,1.25,1.25,1.25),'lines'))
p1_fitness


figs8 <-ggarrange(p1_dnds, p1_nonsyn,p1_fitness,ncol =3, nrow = 1,heights = c(1,1,1),
                  labels = c("A","B","C"),  hjust=-2, vjust=0.8,align = "hv",
                  legend = "bottom", common.legend = T,
                  font.label = list(size = 24, face = "bold"))


ggsave(figs8,filename = file.path(outpath,"figS5_AVC.pdf"),
       width = 18,height =6)




p2_dnds <- ggplot(dat4_S1) +
  geom_point(aes(x = F_Vero, y = dnds_mean_focal,color=type), size=2,shape=19 ,alpha=0.5) +
  xlab("Full-dose vaccine coverage (%)") +
  ylab("Nonsynonymous to synonymous\n divergence ratio") +
  stat_smooth(aes(x=F_Vero,y=dnds_mean_focal,color=type,fill=type),method=lm,se=T,formula = y ~ x)+
  scale_x_continuous(limits=c(0,100),breaks = seq(0,100,25))+
  scale_y_continuous(limits=c(0,9),breaks = seq(0,9,3))+
  stat_cor(aes(x=F_Vero,y=dnds_mean_focal,color=type),size=6,method = "spearman")+
  ggtitle(paste0(""))+
  mytheme+
  scale_fill_manual(name = 'Time period',values = custom_colors)+
  scale_color_manual(name = 'Time period',values = custom_colors)+
  theme(plot.margin=unit(c(0.25,1.25,1.25,1.25),'lines'))

p2_dnds

p2_nonsyn <- ggplot(dat4_S1) +
  geom_point(aes(x = F_Vero, y = muts_nonsyn_focal,color=type), size=2,shape=19 ,alpha=0.5) +
  xlab("Full-dose vaccine coverage (%)") +
  ylab("Number of nonsynonymous\n mutations") +
  stat_smooth(aes(x=F_Vero,y=muts_nonsyn_focal,color=type,fill=type),method=lm,se=T,formula = y ~ x)+
  scale_x_continuous(limits=c(0,100),breaks = seq(0,100,25))+
  scale_y_continuous(limits=c(0,40),breaks = seq(0,40,10))+
  stat_cor(aes(x=F_Vero,y=muts_nonsyn_focal,color=type),size=6,method = "spearman")+
  ggtitle(paste0(""))+
  mytheme+
  scale_fill_manual(name = 'Time period',values = custom_colors)+
  scale_color_manual(name = 'Time period',values = custom_colors)+
  theme(plot.margin=unit(c(0.25,1.25,1.25,1.25),'lines'))

p2_nonsyn

p2_fitness <- ggplot(dat4_S1) +
  geom_point(aes(x = F_Vero, y = fitness_focal,color=type), size=2,shape=19 ,alpha=0.5) +
  xlab("Full-dose vaccine coverage (%)") +
  ylab("Mutational fitness") +
  stat_smooth(aes(x=F_Vero,y=fitness_focal,color=type,fill=type),method=lm,se=T,formula = y ~ x)+
  scale_x_continuous(limits=c(0,100),breaks = seq(0,100,25))+
  scale_y_continuous(limits=c(-0.1,4),breaks = seq(0,4,1))+
  stat_cor(aes(x=F_Vero,y=fitness_focal,color=type),size=6,method = "spearman")+
  ggtitle(paste0(""))+
  mytheme+
  scale_fill_manual(name = 'Time period',values = custom_colors)+
  scale_color_manual(name = 'Time period',values = custom_colors)+
  theme(plot.margin=unit(c(0.25,1.25,1.25,1.25),'lines'))
p2_fitness


figs8_F_Vero <-ggarrange(p2_dnds, p2_nonsyn,p2_fitness,ncol =3, nrow = 1,heights = c(1,1,1),
                  labels = c("A","B","C"),  hjust=-2, vjust=0.8,align = "hv",
                  legend = "bottom", common.legend = T,
                  font.label = list(size = 24, face = "bold"))


ggsave(figs8_F_Vero,filename = file.path(outpath,"figS5_F_Vero.pdf"),
       width = 18,height =6)

