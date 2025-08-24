rm(list = ls())
root_path <- "C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine"
setwd(root_path)
dat1<-read.csv("01_data/subsample_adaptation_vaccine_low.csv")
source("load_package.R")
source("plot_theme.R")
source("load_data_clean.R")


#outpath <-  "02_model-output/all_vaccine_202003_202110"
outpath0 <- "03_output/figure_2_3"

outpath1 <-  paste0(outpath0,"/","dnds_vaccine_202003_202110",sep="")
outpath2 <-  paste0(outpath0,"/","dnds_vaccine_202201_202209",sep="")

outpath3 <-  paste0(outpath0,"/","muts_nonsyn_vaccine_202003_202110",sep="")
outpath4 <-  paste0(outpath0,"/","muts_nonsyn_vaccine_202201_202209")

outpath5 <-  paste0(outpath0,"/","fitness_vaccine_202003_202110",sep="")
outpath6 <-  paste0(outpath0,"/","fitness_vaccine_202201_202209",sep="")

#Plot

S <- read.csv(paste0(root_path,"/",outpath1,"/","vero_dnds_slice_S1.CSV"))
colnames(S) <- c("S_dat_vero_low","S_dat_vero_fit","S_dat_vero_high","vero_coverage" ,"S")

S_O <- read.csv(paste0(root_path,"/",outpath2,"/","vero_dnds_slice_2022_S1.CSV"))
colnames(S_O) <- c("SO_dat_vero_low","SO_dat_vero_fit","SO_dat_vero_high","vero_coverage" ,"SO")

#custom_colors <- c("#ED0000FF", "#00468BFF")
custom_colors <- c("#5f9723", "#00468BFF")

DF1 <- merge(S,S_O,by=c("vero_coverage"),all.y = T)



p1 <-ggplot(DF1,aes(x=vero_coverage,y=S_dat_vero_fit))+
  geom_rect(aes(xmin=54.6, xmax=76, ymin=-Inf, ymax=Inf),fill='grey90',alpha = 0.01)+
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",linewidth=1)+
  geom_vline(aes(xintercept=25),lty='dashed',colour ="#92C5DE",linewidth=1)+
  geom_vline(aes(xintercept=50),lty='dashed',colour ="#4393C3",linewidth=1)+
  geom_vline(aes(xintercept=75),lty='dashed',colour ="#2166AC",linewidth=1)+
  geom_line(aes(x=vero_coverage,y=S_dat_vero_fit,group=1,color="Pre-Omicron period"),linewidth=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = S_dat_vero_low, ymax = S_dat_vero_high,group=1,fill="Pre-Omicron period"),alpha=0.25)+
  geom_line(aes(x=vero_coverage,y=SO_dat_vero_fit,group=1,color="Omicron period"),linewidth=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = SO_dat_vero_low, ymax = SO_dat_vero_high,group=1,fill="Omicron period"),alpha=0.25)+
  scale_x_continuous(limits=c(0,80),breaks = seq(0,75,25))+
  scale_y_continuous(limits=c(-6,3),breaks = seq(-6,3,3))+
  #labels = seq(0,100,25),expand = c(0, 0))+
  labs(x="Adjusted vaccine coverage", y= "Effect on nonsynonymous to \nsynonymous divergence ratio" )+  
  ggtitle("")+
  #scale_x_discrete(labels=NULL)+
  scale_color_manual(name = '',values = custom_colors)+ 
  scale_fill_manual(name = '',values = custom_colors)+ 
  mytheme+
  theme(plot.margin=unit(c(0.75,1,0.75,1.5),'lines'))+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+ 
  theme(panel.grid=element_line(linetype=2,color='grey'))+
  guides(color=guide_legend(nrow=1),fill=guide_legend(nrow=1))+
  theme(axis.title.y.left = element_text( margin = margin(0,0.1,0,0,'cm')),
        axis.title.x = element_text(margin = margin(0.1,0,0,0,'cm')))+
  theme(legend.position = "bottom")

p1



S1 <- read.csv(paste0(root_path,"/",outpath3,"/","vero_muts_nonsyn_slice_S1.CSV"))
colnames(S1) <- c("S1_dat_vero_low","S1_dat_vero_fit","S1_dat_vero_high","vero_coverage" ,"S1")

S1_O <- read.csv(paste0(root_path,"/",outpath4,"/","vero_muts_nonsyn_slice_2022_S1.CSV"))
colnames(S1_O) <- c("S1O_dat_vero_low","S1O_dat_vero_fit","S1O_dat_vero_high","vero_coverage" ,"S1O")



DF2 <- merge(S1,S1_O,by=c("vero_coverage"),all.y = T)

p2 <-ggplot(DF2,aes(x=vero_coverage,y=S1_dat_vero_fit))+
  geom_rect(aes(xmin=59.8, xmax=76, ymin=-Inf, ymax=Inf),fill='grey90',alpha = 0.01)+
  geom_vline(aes(xintercept=25),lty='dashed',colour ="#92C5DE",linewidth=1)+
  geom_vline(aes(xintercept=50),lty='dashed',colour ="#4393C3",linewidth=1)+
  geom_vline(aes(xintercept=75),lty='dashed',colour ="#2166AC",linewidth=1)+
  
  geom_line(aes(x=vero_coverage,y=S1_dat_vero_fit,group=1,color="Pre-Omicron period"),size=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = S1_dat_vero_low, ymax = S1_dat_vero_high,group=1,fill="Pre-Omicron"),alpha=0.25)+
  geom_line(aes(x=vero_coverage,y=S1O_dat_vero_fit,group=1,color="Omicron period"),size=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = S1O_dat_vero_low, ymax = S1O_dat_vero_high,group=1,fill="Omicron"),alpha=0.25)+
  scale_x_continuous(limits=c(0,80),breaks = seq(0,75,25))+
  scale_y_continuous(limits=c(-18,18),breaks = seq(-18,18,9))+
  #labels = seq(0,100,25),expand = c(0, 0))+
  labs(x="Adjusted vaccine coverage", y= "Effect on number \nof nonsynonymous mutations" )+  
  ggtitle("")+
  #scale_x_discrete(labels=NULL)+
  scale_color_manual(name = '',values = custom_colors)+ 
  scale_fill_manual(name = '',values = custom_colors)+ 
  mytheme+
  theme(panel.grid=element_line(linetype=2,color='grey'))+
  guides(color=guide_legend(nrow=1),fill=guide_legend(nrow=1))+
  theme(legend.position = "bottom")+
  theme(plot.margin=unit(c(0.75,1,0.75,1.5),'lines'))+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+ 
  theme(axis.title.y.left = element_text( margin = margin(0,0.1,0,0,'cm')),
        axis.title.x = element_text(margin = margin(0.1,0,0,0,'cm')))+
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",size=1)
p2



S2 <- read.csv(paste0(root_path,"/",outpath5,"/","vero_fitness_slice_S1.CSV"))
colnames(S2) <- c("S2_dat_vero_low","S2_dat_vero_fit","S2_dat_vero_high","vero_coverage" ,"S2")

S2_O <- read.csv(paste0(root_path,"/",outpath6,"/","vero_fitness_slice_2022_S1.CSV"))
colnames(S2_O) <- c("S2O_dat_vero_low","S2O_dat_vero_fit","S2O_dat_vero_high","vero_coverage" ,"S2O")


DF3 <- merge(S2,S2_O,by=c("vero_coverage"),all.y = T)

p3 <-ggplot(DF3,aes(x=vero_coverage,y=S2_dat_vero_fit))+
  geom_rect(aes(xmin=62.7, xmax=76, ymin=-Inf, ymax=Inf),fill='grey90',alpha = 0.01)+
  geom_vline(aes(xintercept=25),lty='dashed',colour ="#92C5DE",linewidth=1)+
  geom_vline(aes(xintercept=50),lty='dashed',colour ="#4393C3",linewidth=1)+
  geom_vline(aes(xintercept=75),lty='dashed',colour ="#2166AC",linewidth=1)+
  
  geom_line(aes(x=vero_coverage,y=S2_dat_vero_fit,group=1,color="Pre-Omicron period"),size=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = S2_dat_vero_low, ymax = S2_dat_vero_high,group=1,fill="Pre-Omicron"),alpha=0.25)+
  scale_x_continuous(limits=c(0,80),breaks = seq(0,75,25))+
  geom_line(aes(x=vero_coverage,y=S2O_dat_vero_fit,group=1,color="Omicron period"),size=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = S2O_dat_vero_low, ymax = S2O_dat_vero_high,group=1,fill="Omicron"),alpha=0.25)+
  scale_y_continuous(limits=c(-2,2),breaks = seq(-2,2,1))+
  #labels = seq(0,100,25),expand = c(0, 0))+
  labs(x="Adjusted vaccine coverage", y= "Effect on mutational fitness" )+  
  ggtitle("")+
  #scale_x_discrete(labels=NULL)+
  scale_color_manual(name = '',values = custom_colors)+ 
  scale_fill_manual(name = '',values = custom_colors)+ 
  mytheme+
  theme(panel.grid=element_line(linetype=2,color='grey'))+
  guides(color=guide_legend(nrow=1),fill=guide_legend(nrow=1))+
  theme(legend.position = "bottom")+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+ 
  #theme(legend.position = "none")+
  theme(plot.margin=unit(c(0.75,1,0.75,1.5),'lines'))+
  theme(axis.title.y.left = element_text( margin = margin(0,0.1,0,0,'cm')),
        axis.title.x = element_text(margin = margin(0.1,0,0,0,'cm')))+
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",size=1)
p3




p<-ggarrange(p1, p2, p3,ncol =3, nrow = 1,widths = c(1,1,1),heights = c(1,1,1),
             labels = c("A","B","C"), align = "hv",
             legend = "bottom", common.legend = T,
             font.label = list(size = 30, face = "bold"))
p

ggsave(p,filename = file.path(outpath0,"figure_2.pdf"),
       width = 18.2,height =6.2)
