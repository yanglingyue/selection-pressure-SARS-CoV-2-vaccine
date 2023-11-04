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

dir.create("04_extended_analysis/plot_output/Figure_S6")
outpath <-  "04_extended_analysis/plot_output/Figure_S6/"
dat1_SF_202110 <- subset(dat1,dat1$gene== "S2"&dat1$month=="2021-10")


####################################fig.xx##############################################
#####Hot plot of adjusted vaccine coverage from March 2020 to September 2022######
max(dat1$FW_BW_P,na.rm=T)
cols <- brewer.pal(8,'Blues')
#cols <- brewer.pal(3,'Purples')
dat1$country <- factor(dat1$country,
                           levels = dat1_SF_202110$country[order(dat1_SF_202110$F_Vero, decreasing = F)],ordered=TRUE)

p1 <- ggplot(dat1, aes(month,country)) +
  geom_tile(aes(fill=FW_BW_P),color = "white",size=0.3) +
  scale_fill_gradientn(colors=cols, limits = c(0,80),na.value="grey88") +
  #geom_vline(aes(xintercept="2021-01"),lty='solid',colour ="black",size=1)+
  geom_vline(aes(xintercept="2021-11"),lty="dashed",colour ="red",size=2)+
  #Rotating labels
  coord_flip()+
  #theme_bw()+
  #theme(panel.grid=element_blank(),axis.line=element_line(size=0.3,colour="black"))+
  labs(x="",y="") + 
  scale_x_discrete(breaks=c("2020-03","2020-06","2020-09","2020-12","2021-03","2021-06",
                            "2021-09","2021-12","2022-03","2022-06","2022-09"),
                   labels=c("Mar \n2020","Jun \n2020","Sep \n2020","Dec \n2020","Mar \n2021","Jun \n2021",
                            "Sep \n2021","Dec \n2021","Mar \n2022","Jun \n2022","Sep \n2022"))+
  mytheme+
  theme(axis.text.x = element_text(size = 20, color="black",vjust = 0.5,hjust = 1, angle = 90))+
  theme(legend.key.size = unit(0.5,"cm"),legend.key.width = unit(3,"cm"),legend.key.height  = unit(0.35,"cm"))+
  theme(plot.margin=unit(c(1,1,0.25,0.5),'lines'))+
  guides(fill=guide_colorbar("Adjusted vaccine coverage",title.position = "top",
                             title.theme = element_text(size = 20, angle = 0,hjust = 0.5),
                             label.theme= element_text(size = 20),
                             title.vjust = 0.9)) 
p1 



max(dat1$FW_BW_P_jama,na.rm=T)
p2 <- ggplot(dat1, aes(month,country)) +
  geom_tile(aes(fill=FW_BW_P_jama),color = "white",size=0.3) +
  scale_fill_gradientn(colors=cols, limits = c(0,80),na.value="grey88") +
  #geom_vline(aes(xintercept="2021-01"),lty='solid',colour ="black",size=1)+
  geom_vline(aes(xintercept="2021-11"),lty="dashed",colour ="red",size=2)+
  #Rotating labels
  coord_flip()+
  #theme_bw()+
  #theme(panel.grid=element_blank(),axis.line=element_line(size=0.3,colour="black"))+
  labs(x="",y="") + 
  scale_x_discrete(breaks=c("2020-03","2020-06","2020-09","2020-12","2021-03","2021-06",
                            "2021-09","2021-12","2022-03","2022-06","2022-09"),
                   labels=c("Mar \n2020","Jun \n2020","Sep \n2020","Dec \n2020","Mar \n2021","Jun \n2021",
                            "Sep \n2021","Dec \n2021","Mar \n2022","Jun \n2022","Sep \n2022"))+
  mytheme+
  theme(axis.text.x = element_text(size = 20, color="black",vjust = 0.5,hjust = 1, angle = 90))+
  theme(legend.key.size = unit(0.5,"cm"),legend.key.width = unit(3,"cm"),legend.key.height  = unit(0.35,"cm"))+
  theme(plot.margin=unit(c(1,1,0.25,0.5),'lines'))+
  guides(fill=guide_colorbar("Adjusted vaccine coverage",title.position = "top",
                             title.theme = element_text(size = 20, angle = 0,hjust = 0.5),
                             label.theme= element_text(size = 20),
                             title.vjust = 0.9)) 
p2


figs6 <-ggarrange(p1, p2,ncol =1, nrow = 2,heights = c(1,1),
                  labels = c("A","B"),  hjust= -1, vjust=-1,align = "hv",
                  legend = "bottom", common.legend = T,
                  font.label = list(size = 24, face = "bold"))
figs6 

ggsave(figs6,filename = file.path(outpath,"figS6.pdf"),
       width = 20,height =18)

