# # install.packages("lubridate")
rm(list = ls())
setwd("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine")
dat<-read.csv("01_data/subsample_dnds_vaccine_diversity_low_death.csv")
dat1 <- dat

source("load_package.R")
source("plot_theme.R")
source("load_data_clean.R")

dir.create("04_extended_analysis/plot_output/Figure_S1")
outpath <-  "04_extended_analysis/plot_output/Figure_S1/"


dat1_S2= subset(dat1,dat1$gene== "S2")
dat1_S2_202110 <- subset(dat1,dat1$gene== "S2"&dat1$month=="2021-10")
dat1_S2$country <- factor(dat1_S2$country,
                                 levels = dat1_S2_202110$country[order(dat1_S2_202110$F_Vero, decreasing = F)],ordered=TRUE)

pal <- rev(brewer.pal(11, "RdBu"))
levels <- pretty(c(0, 10), 80)
#col2 <- colorRampPalette(pal[6:10])
#col1 <- colorRampPalette(c("#5c4591",'#a786b9',"#d4b9da","#e7e1ef"))
#col2 <- colorRampPalette(brewer.pal(10,"RdYlBu")[6:9])
col1 <-  colorRampPalette(brewer.pal(11,"RdBu")[2:6])
col2 <- colorRampPalette(brewer.pal(10,"RdYlBu")[6:10])

#cols <- c(col1(sum(levels < 1)), col2(sum(levels > 1)))
cols <- c(col1(sum(levels < 1)),  col2(sum(levels > 1)))
#colors <- c("#d4b9da","#e7e1ef",brewer.pal(10,"RdYlBu")[6:10])

ps1 <- ggplot(dat1_S2, aes(country,month)) +
  geom_tile(aes(fill=dn_ds_mean),color = "white",size=0.3) +
  #geom_text(aes(label = values), color = "black",size=3) + 
  scale_fill_gradientn(colors = cols, limits = c(0, 10), na.value = "grey80",
                       breaks = seq(0, 10, by = 1)) + 
  scale_y_discrete(breaks=c("2020-03","2020-09","2021-03",
                            "2021-09","2022-03","2022-09"),
                   labels=c("Mar \n2020","Sep \n2020","Mar \n2021",
                            "Sep \n2021","Mar \n2022","Sep \n2022"),expand = c(0, 0))+
  geom_hline(aes(yintercept="2021-01"),lty='solid',colour ="black",size=0.5)+
  geom_hline(aes(yintercept="2021-11"),lty='solid',colour ="black",size=0.5)+
  labs(x="",y="") + 
  ggtitle("S2")+
  mytheme+
  theme(axis.text.x = element_text(size = 20, color="black",vjust = 0.5,hjust = 1, angle = 90))+
  theme(legend.position = "right",    legend.key.size=unit(2,'cm'),
        legend.key.width = unit(3,"cm"),legend.key.height  = unit(0.5,"cm"))+
  guides(fill = guide_colorbar("dN/dS ratio", title.position = "right",
                               title.theme = element_text(size = 20, angle = 90, hjust = 0.5),
                               label.theme = element_text(size = 20),
                               title.vjust = 0.9,ticks.linewidth = 0.5,lty='solid',  # Set tick line width
                               ticks.colour = "black")) +
  guides(fill = "none",color="none")+
  theme(plot.title = element_text(hjust = 0.5,size=24, face = "bold"))

ps1



dat1_S= subset(dat1,dat1$gene== "S")
dat1_S_202110 <- subset(dat1,dat1$gene== "S"&dat1$month=="2021-10")
dat1_S$country <- factor(dat1_S$country,
                          levels = dat1_S_202110$country[order(dat1_S_202110$F_Vero, decreasing = F)],ordered=TRUE)

ps2 <- ggplot(dat1_S, aes(country,month)) +
  geom_tile(aes(fill=dn_ds_mean),color = "white",size=0.3) +
  #geom_text(aes(label = values), color = "black",size=3) + 
  scale_fill_gradientn(colors = cols, limits = c(0, 10), na.value = "grey80",
                       breaks = seq(0, 10, by = 1)) + 
  scale_y_discrete(breaks=c("2020-03","2020-09","2021-03",
                            "2021-09","2022-03","2022-09"),
                   labels=c("Mar \n2020","Sep \n2020","Mar \n2021",
                            "Sep \n2021","Mar \n2022","Sep \n2022"),expand = c(0, 0))+
  geom_hline(aes(yintercept="2021-01"),lty='solid',colour ="black",size=0.5)+
  geom_hline(aes(yintercept="2021-11"),lty='solid',colour ="black",size=0.5)+
  labs(x="",y="") + 
  ggtitle("S")+
  mytheme+
  theme(axis.text.x = element_text(size = 20, color="black",vjust = 0.5,hjust = 1, angle = 90))+
  theme(legend.position = "bottom",    legend.key.size=unit(2,'cm'),
        legend.key.width = unit(4,"cm"),legend.key.height  = unit(0.5,"cm"))+
  guides(fill = guide_colorbar("Effect on ratio of nonsynonymous to synonymous divergence", title.position = "top",
                               title.theme = element_text(size = 20, angle = 0, hjust = 0.5),
                               label.theme = element_text(size = 20),
                               title.vjust = 0.9,ticks.linewidth = 0.5,lty='solid',  # Set tick line width
                               ticks.colour = "black")) +

  theme(plot.title = element_text(hjust = 0.5,size=24, face = "bold"))

ps2



library(ggpubr)
library(egg)
library(patchwork)


patchwork =  ps1 + ps2 +plot_layout(nrow = 2,heights = c(1, 1))

# Remove title from second subplot
patchwork[[1]] = patchwork[[1]] + theme(axis.text.x = element_blank(),
                                        axis.ticks.x = element_blank(),
                                        axis.title.x = element_blank() )

patchwork[[2]] = patchwork[[2]] 
               

patchwork


ggsave(patchwork,filename = file.path(paste0(outpath,"fig_S1.pdf")),
       width = 20,height =18)







