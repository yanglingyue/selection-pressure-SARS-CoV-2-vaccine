setwd("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine")
# Load packages (assumed already installed)
rm(list = ls())
library(data.table);library(corrplot);library(ggcorrplot);library(RColorBrewer)
library(ggpubr);

dat<-read.csv("01_data/subsample_dnds_vaccine_diversity_low_death.csv")
dat1 <- dat

source("load_package.R")
source("plot_theme.R")
source("load_data_clean.R")

dir.create("04_extended_analysis/plot_output/Figure_S7")
outpath <-  "04_extended_analysis/plot_output/Figure_S7/"

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

dat4_S= subset(dat4,dat4$gene== "S")
dat4_S1= subset(dat4,dat4$gene== "S1")
dat4_S2= subset(dat4,dat4$gene== "S2")

p1_S <- ggplot(dat4_S) +
  geom_point(aes(x = FW_BW_P, y = dn_ds_mean,color=type), size=1,shape=19 ,alpha=0.5) +
  xlab("") +
  ylab("Ratio of nonsynonymous \nto synonymous divergence") +
  stat_smooth(aes(x=FW_BW_P,y=dn_ds_mean,color=type,fill=type),method=lm,se=T,formula = y ~ x)+
  scale_x_continuous(limits=c(0,80),breaks = seq(0,75,25))+
  scale_y_continuous(limits=c(0,6.3),breaks = seq(0,6,2))+
  stat_cor(aes(x=FW_BW_P,y=dn_ds_mean,color=type),method = "spearman",size=6)+
  ggtitle(paste0("S"))+
  mytheme+
  scale_fill_manual(name = 'Time period',values = custom_colors)+
  scale_color_manual(name = 'Time period',values = custom_colors)+
  theme(plot.margin=unit(c(0.25,0.5,1.5,0.25),'lines'))

p1_S

p1_S1 <- ggplot(dat4_S1) +
  geom_point(aes(x = FW_BW_P, y = dn_ds_mean,color=type), size=1,shape=19 ,alpha=0.5) +
  xlab("Ajusted vaccine coverage") +
  ylab(" ") +
  stat_smooth(aes(x=FW_BW_P,y=dn_ds_mean,color=type,fill=type),method=lm,se=T,formula = y ~ x)+
  scale_x_continuous(limits=c(0,80),breaks = seq(0,75,25))+
  scale_y_continuous(limits=c(0,10.5),breaks = seq(0,10,2))+
  stat_cor(aes(x=FW_BW_P,y=dn_ds_mean,color=type),method = "spearman",size=6)+
  ggtitle(paste0("S1"))+
  mytheme+
  scale_fill_manual(name = 'Time period',values = custom_colors)+
  scale_color_manual(name = 'Time period',values = custom_colors)+
  theme(plot.margin=unit(c(0.25,0.5,1.5,0.25),'lines'))

p1_S1

p1_S2 <- ggplot(dat4_S2) +
  geom_point(aes(x = FW_BW_P, y = dn_ds_mean,color=type), size=1,shape=19 ,alpha=0.5) +
  xlab("") +
  ylab(" ") +
  stat_smooth(aes(x=FW_BW_P,y=dn_ds_mean,color=type,fill=type),method=lm,se=T,formula = y ~ x)+
  scale_x_continuous(limits=c(0,80),breaks = seq(0,75,25))+
  scale_y_continuous(limits=c(0,1.05),breaks = seq(0,1,0.2))+
  stat_cor(aes(x=FW_BW_P,y=dn_ds_mean,color=type),method = "spearman",size=6)+
  ggtitle(paste0("S2"))+
  mytheme+
  scale_fill_manual(name = 'Time period',values = custom_colors)+
  scale_color_manual(name = 'Time period',values = custom_colors)+
  theme(plot.margin=unit(c(0.25,0.5,1.5,0.25),'lines'))
p1_S2


figs8 <-ggarrange(p1_S, p1_S1,p1_S2,ncol =3, nrow = 1,heights = c(1,1,1),
                  labels = c("A","B","C"),  hjust=-2, vjust=0.8,align = "hv",
                  legend = "bottom", common.legend = T,
                  font.label = list(size = 24, face = "bold"))


ggsave(figs8,filename = file.path(outpath,"figS7_AVC.pdf"),
       width = 15,height =5.5)




p2_S <- ggplot(dat4_S) +
  geom_point(aes(x = F_Vero, y = dn_ds_mean,color=type), size=1,shape=19 ,alpha=0.5) +
  xlab("") +
  ylab("Ratio of nonsynonymous \nto synonymous divergence") +
  stat_smooth(aes(x=F_Vero,y=dn_ds_mean,color=type,fill=type),method=lm,se=T,formula = y ~ x)+
  scale_x_continuous(limits=c(0,100),breaks = seq(0,100,25))+
  scale_y_continuous(limits=c(0,6.3),breaks = seq(0,6,2))+
  stat_cor(aes(x=F_Vero,y=dn_ds_mean,color=type),method = "spearman")+
  ggtitle(paste0("S"))+
  mytheme+
  scale_fill_manual(name = 'Time period',values = custom_colors)+
  scale_color_manual(name = 'Time period',values = custom_colors)+
  theme(plot.margin=unit(c(0.25,0.5,1.5,0.25),'lines'))

p2_S

p2_S1 <- ggplot(dat4_S1) +
  geom_point(aes(x = F_Vero, y = dn_ds_mean,color=type), size=1,shape=19 ,alpha=0.5) +
  xlab("Ajusted vaccine coverage") +
  ylab(" ") +
  stat_smooth(aes(x=F_Vero,y=dn_ds_mean,color=type,fill=type),method=lm,se=T,formula = y ~ x)+
  scale_x_continuous(limits=c(0,100),breaks = seq(0,100,25))+
  scale_y_continuous(limits=c(0,10.5),breaks = seq(0,10,2))+
  stat_cor(aes(x=F_Vero,y=dn_ds_mean,color=type),method = "spearman")+
  ggtitle(paste0("S1"))+
  mytheme+
  scale_fill_manual(name = 'Time period',values = custom_colors)+
  scale_color_manual(name = 'Time period',values = custom_colors)+
  theme(plot.margin=unit(c(0.25,0.5,1.5,0.25),'lines'))

p2_S1

p2_S2 <- ggplot(dat4_S2) +
  geom_point(aes(x = F_Vero, y = dn_ds_mean,color=type), size=1,shape=19 ,alpha=0.5) +
  xlab("") +
  ylab(" ") +
  stat_smooth(aes(x=F_Vero,y=dn_ds_mean,color=type,fill=type),method=lm,se=T,formula = y ~ x)+
  scale_x_continuous(limits=c(0,100),breaks = seq(0,100,25))+
  scale_y_continuous(limits=c(0,1.05),breaks = seq(0,1,0.2))+
  stat_cor(aes(x=F_Vero,y=dn_ds_mean,color=type),method = "spearman")+
  ggtitle(paste0("S2"))+
  mytheme+
  scale_fill_manual(name = 'Time period',values = custom_colors)+
  scale_color_manual(name = 'Time period',values = custom_colors)+
  theme(plot.margin=unit(c(0.25,0.5,1.5,0.25),'lines'))
p2_S2


figs8_F_Vero <-ggarrange(p2_S, p2_S1,p2_S2,ncol =3, nrow = 1,heights = c(1,1,1),
                  labels = c("A","B","C"),  hjust=-2, vjust=0.8,align = "hv",
                  legend = "bottom", common.legend = T,
                  font.label = list(size = 24, face = "bold"))


ggsave(figs8_F_Vero,filename = file.path(outpath,"figS7_F_Vero.pdf"),
       width = 15,height =5.5)

