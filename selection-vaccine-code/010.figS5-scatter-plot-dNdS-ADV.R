setwd("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/selection-vaccine-code-v8")
# Load packages (assumed already installed)
rm(list = ls())
library(data.table);library(corrplot);library(ggcorrplot);library(RColorBrewer)
library(ggpubr);

custom_colors <- list()
colors_dutch <- c(
           '#EE5A24','#45aaf2','#A3CB38','#1289A7','#009432',
           '#0652DD',"#c4daee","#E64B35BF",'#F79F1F','#EA2027',
           "#9084b4",'#833471','#1B1464','#5758BB','#f7f1e3'
)
colors_spanish <- c(
            "#B09C85B2","#bcd1a2",'#6F1E51',"#c4daee",'#34ace0',
              '#33d9b2','#2c2c54','#474787','#aaa69d','#227093',
             '#ff5252','#ff793f','#f7f1e3','#ffb142','#ffda79',
             '#EA2027','#cd6133','#84817a','#cc8e35','#ccae62')

custom_colors$discrete <- c(colors_dutch, colors_spanish)
custom_colors$cell_cycle <- c("#00468BFF","#ED0000FF",'#F79F1F',"#42B540B2","#7570b3",
                              "#1b9e77",'#b33939','#ffb142',"#ED0000FF","#4DBBD5FF")

#custom_colors$cell_cycle <- c("#45aaf2",'#ffb142',"#EA2027", "#504190")

mytheme <- theme_minimal()+
  theme(
    #panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    #panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    # plot.title=element_text(size=12, vjust = 0.5, hjust = 0.5),
    legend.position = c("bottom"),
    legend.key.size = unit(0.25,"cm"),
    legend.key.width = unit(1,"cm"),
    legend.key.height  = unit(0.5,"cm"),
    #legend.spacing.x = unit(1, 'cm'),
    legend.margin = margin(0, 0, 0, 0),
    legend.background = element_rect(fill=NA, size=0,color=NA),
    legend.text=element_text(size=18),
    legend.title=element_text(face="bold",size=18),
    axis.line.y=element_line(linetype=1,color='black',size=0.5),
    axis.line.x=element_line(linetype=1,color='black',size=0.5),
    axis.ticks = element_line(linetype=2,color='black'),
    panel.grid=element_line(linetype=1,color='grey'),
    plot.title = element_text(hjust = 0.5,size=20, face = "bold"),
    axis.title.y.left = element_text(size = 18,color="black",vjust=0,face = "bold"),
    axis.title.y.right =element_text(size = 18,color="black",vjust=0,angle=90),
    axis.title.x = element_text(size = 18, color="black",face = "bold", vjust = 0),
    axis.text.y.left = element_text(size = 18,color="black",vjust = 0.5,angle = 0),
    axis.text.y.right = element_text(size = 18,color="black",vjust = 0.5,angle = 0),
    axis.text.x = element_text(size = 18, color="black",vjust = 0.5, angle = 0),
    axis.ticks.length=unit(0.15,"cm"),
    axis.ticks.y=element_line(color="black",size=.5),
    axis.ticks.x=element_line(color="black",size=.5),
    plot.margin=unit(c(0.5,0.5,0.5,0.5),'lines'),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )




####################################fig.xx##############################################
#####Scatter plot of adjusted vaccine coverage and selection pressure on S regions######
###Before the emergence of Omicron variants and during Omicron variants were dominant###
dat1<-  fread("01_data/input_dnds_subsample_low_200_BV_Nextstrain.csv")
dat1<- as.data.frame(dat1)

y <- which(dat1$month=="2020-01"|dat1$month=="2020-02")
dat1 <-dat1[-c(y),]

x <- which(dat1$country=="Puerto Rico"|dat1$country=="Canary Islands"|dat1$country=="French Guiana"|
             dat1$country=="Reunion"|dat1$country=="Curacao"|dat1$country=="North Macedonia"|
             dat1$country=="Slovakia")
dat1 <-dat1[-x,]


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

p8_S <- ggplot(dat4_S) +
  geom_point(aes(x = FW_BW_P, y = dn_ds_mean,color=type), size=1,shape=19 ,alpha=0.5) +
  xlab("") +
  ylab("Ratio of nonsynonymous \nto synonymous divergence") +
  stat_smooth(aes(x=FW_BW_P,y=dn_ds_mean,color=type,fill=type),method=lm,se=T,formula = y ~ x)+
  scale_x_continuous(limits=c(0,80),breaks = seq(0,75,25))+
  scale_y_continuous(limits=c(0,6.3),breaks = seq(0,6,2))+
  stat_cor(aes(x=FW_BW_P,y=dn_ds_mean,color=type),method = "spearman")+
  ggtitle(paste0("S"))+
  mytheme+
  scale_fill_manual(name = 'Time period',values = custom_colors$cell_cycle)+
  scale_color_manual(name = 'Time period',values = custom_colors$cell_cycle)+
  theme(plot.margin=unit(c(0.25,0.5,1.5,0.25),'lines'))

p8_S

p8_S1 <- ggplot(dat4_S1) +
  geom_point(aes(x = FW_BW_P, y = dn_ds_mean,color=type), size=1,shape=19 ,alpha=0.5) +
  xlab("Ajusted vaccine coverage") +
  ylab(" ") +
  stat_smooth(aes(x=FW_BW_P,y=dn_ds_mean,color=type,fill=type),method=lm,se=T,formula = y ~ x)+
  scale_x_continuous(limits=c(0,80),breaks = seq(0,75,25))+
  scale_y_continuous(limits=c(0,10.5),breaks = seq(0,10,2))+
  stat_cor(aes(x=FW_BW_P,y=dn_ds_mean,color=type),method = "spearman")+
  ggtitle(paste0("S1"))+
  mytheme+
  scale_fill_manual(name = 'Time period',values = custom_colors$cell_cycle)+
  scale_color_manual(name = 'Time period',values = custom_colors$cell_cycle)+
  theme(plot.margin=unit(c(0.25,0.5,1.5,0.25),'lines'))

p8_S1

p8_S2 <- ggplot(dat4_S2) +
  geom_point(aes(x = FW_BW_P, y = dn_ds_mean,color=type), size=1,shape=19 ,alpha=0.5) +
  xlab("") +
  ylab(" ") +
  stat_smooth(aes(x=FW_BW_P,y=dn_ds_mean,color=type,fill=type),method=lm,se=T,formula = y ~ x)+
  scale_x_continuous(limits=c(0,80),breaks = seq(0,75,25))+
  scale_y_continuous(limits=c(0,1.05),breaks = seq(0,1,0.2))+
  stat_cor(aes(x=FW_BW_P,y=dn_ds_mean,color=type),method = "spearman")+
  ggtitle(paste0("S2"))+
  mytheme+
  scale_fill_manual(name = 'Time period',values = custom_colors$cell_cycle)+
  scale_color_manual(name = 'Time period',values = custom_colors$cell_cycle)+
  theme(plot.margin=unit(c(0.25,0.5,1.5,0.25),'lines'))
p8_S2


figs5 <-ggarrange(p8_S, p8_S1,p8_S2,ncol =3, nrow = 1,heights = c(1,1,1),
                  labels = c("a","b","c"),  hjust=-2, vjust=0.8,align = "hv",
                  legend = "bottom", common.legend = T,
                  font.label = list(size = 24, face = "bold"))
figs5 
ggsave("03_fig-main/figS5.pdf", figs5,width = 15,height =5.5)


