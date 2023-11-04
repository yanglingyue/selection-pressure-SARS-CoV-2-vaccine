# # install.packages("lubridate")
rm(list = ls())
setwd("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine")
dat<-read.csv("01_data/subsample_dnds_vaccine_diversity_low_death.csv")
dat1 <- dat

source("load_package.R")
source("plot_theme.R")
source("load_data_clean.R")

dir.create("03_output/figure_1")
outpath <-  "03_output/figure_1/"


dat1_S1_202110= subset(dat1,dat1$gene== "S1"&dat1$month=="2021-10")


region <- read.csv("02_sup-data/region_code.csv")
colnames(region) <- c("iso_code","continent","country")
dat1_S1_202110 <- merge(dat1_S1_202110,region,by=c("country"),all.x = T)


#custom_colors <- c("#c25b88","#96c597","#806ca9","#b1dbf3","#425da2","#e8b50f")
custom_colors <- c("#3cb371","#ee1289","#1874cd","#7d26cd","#eec900","#ee7600")
custom_colors <- c("#eb4c1c","#3cb371","#7d26cd","#1874cd","#425da2","#e8b50f")
custom_colors <- c("#3cb371","#c25b88","#1874cd","#7d26cd","#eec900","#ee7600")
custom_colors <- c("#4daf4a","#4DBBD5FF","#7570b3","#ee7600","#eec900","#c25b88")
custom_colors <- c("#eb4c1c","#f7b455","#6cbb73","#2d9596","#425da2","#5c4591")
#custom_colors <- c("#fedb3f","#da2828","#6cbb73","#7aacdf","#a3298e","#5c4591")
custom_colors <- c("#e4b986","#ffdc75","#f29f5f","#d06a6c","#a1799f","#5f1c83")
custom_colors <- c("#e4b986","#ffdc75","#f29f5f","#7eb780","#2f9798","#406175")

library(mgcv)
dat1_S1= subset(dat1,dat1$gene== "S1")
pal <- rev(brewer.pal(11, "RdBu"))
levels <- pretty(c(0, 10), 80)
#col2 <- colorRampPalette(pal[6:10])
#col1 <- colorRampPalette(pal[6:10])

#col1 <- colorRampPalette(c("#5c4591",'#a786b9',"#d4b9da","#e7e1ef"))
col1 <-  colorRampPalette(brewer.pal(11,"RdBu")[2:6])
col2 <- colorRampPalette(brewer.pal(10,"RdYlBu")[6:10])

#cols <- c(col1(sum(levels < 1)), col2(sum(levels > 1)))
cols <- c(col1(sum(levels < 1)),  col2(sum(levels > 1)))

dat1_S1_202209 <- subset(dat1,dat1$gene== "S1"&dat1$month=="2021-10")

dat1_S1_202209 <- merge(dat1_S1_202209,region,by=c("country"),all.x = T)

dat1_S1$country <- factor(dat1_S1$country,
                                 levels = dat1_S1_202209$country[order(dat1_S1_202209$F_Vero, decreasing = F)],ordered=TRUE)

#colors <- c("#d4b9da","#e7e1ef",brewer.pal(10,"RdYlBu")[6:10])

ps1 <- ggplot(dat1_S1, aes(country,month)) +
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
  #scale_fill_stepsn("Difference (%)", labels = seq(0,10,0.5),
  #                  breaks=seq(0,10,0.5), limits = c(0,10),
  #                  colours = cols) +
  #coord_flip()+
  labs(x="",y="") + 
  mytheme+
  theme(axis.text.x = element_text(size = 18, color="black",vjust = 0.5,hjust = 1, angle = 90),
        axis.text.y = element_text(size = 18, color = "black", hjust = 0))+
  theme(legend.position = "bottom",    legend.key.size=unit(2,'cm'),
        legend.key.width = unit(3,"cm"),legend.key.height  = unit(0.3,"cm"))+
  guides(fill = guide_colorbar("Ratio of nonsynonymous to synonymous divergence", title.position = "right",
                               title.theme = element_text(size = 18, angle = 0, hjust = 0.5),
                               label.theme = element_text(size = 18),
                               title.vjust = 0.9,ticks.linewidth = 0.5,lty='solid',  # Set tick line width
                               ticks.colour = "black")) +
  #guides(fill = "none",color="none")+
  theme(plot.title = element_text(hjust = 0.5,size=24, face = "bold"))

ps1


dat1_S1_202110$country <- factor(dat1_S1_202110$country,
                          levels = dat1_S1_202110$country[order(dat1_S1_202110$F_Vero, decreasing = F)],ordered=TRUE)


ps2 <- ggplot(dat1_S1_202110,aes(x=country,y=F_Vero,group=1))+
  geom_bar(aes(x=country,y=F_Vero,fill=continent),stat="identity",alpha=0.75)+
  scale_fill_manual("",values=custom_colors)+
  #scale_y_continuous(limits=c(0,102),breaks = seq(0,100,25))+
  #guides(fill = "none",color="none")+
  guides(fill=guide_legend(nrow=1,byrow = T))+
  xlab("") +
  ylab("Full-dose vaccine\n coverage (%)") +
  geom_hline(aes(yintercept=25),lty='dashed',colour ="grey25",size=0.5)+
  geom_hline(aes(yintercept=50),lty='dashed',colour ="grey25",size=0.5)+
  geom_hline(aes(yintercept=75),lty='dashed',colour ="grey25",size=0.5)+
  mytheme+
  theme(axis.text.x = element_text(size = 18, color="black",vjust = 0.5,hjust = 1, angle = 90))+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())
#ggtitle("Natural immunity")+
ps2


library(ggpubr)
library(egg)
library(patchwork)


patchwork =  ps1 + ps2 +plot_layout(nrow = 2,heights = c(60, 40))

# Remove title from second subplot
patchwork[[1]] = patchwork[[1]] + theme(axis.text.x = element_blank(),
                                        axis.ticks.x = element_blank(),
                                        axis.title.x = element_blank() )

patchwork[[2]] = patchwork[[2]] 
 #+theme(legend.position = ("none"),legend.justification = c("bottom"))


patchwork


ggsave(patchwork,filename = file.path(paste0(outpath,"fig1-a-b.pdf")),
       width = 18,height =11)







