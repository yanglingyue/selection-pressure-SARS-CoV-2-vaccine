# # install.packages("lubridate")
rm(list = ls())
setwd("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine")

source("load_package.R")
source("plot_theme.R")

outpath0 <-  "04_extended_analysis/plot_output/Figure_S16"
dir.create(outpath0)
outpath <-  "03_output/figure_4/Supplementary_fig4_correlation/cor_202110_"

dat_F<-read.csv("01_data/subsample_dnds_vaccine_diversity_low_death.csv")
dat1_SF_202110 <- subset(dat_F,dat_F$gene== "S2"&dat_F$month=="2021-10")

cor_S <- read.csv(paste0(outpath,"S","_","lag2.csv"))
cor_S1 <- read.csv(paste0(outpath,"S1","_","lag2.csv"))
cor_S2 <- read.csv(paste0(outpath,"S2","_","lag2.csv"))
cor_gene <- rbind(cor_S,cor_S1,cor_S2)

cor_gene<- cor_gene[order(cor_gene$Coverage,cor_gene$country,decreasing=F),]
cor_gene$country <-factor(cor_gene$country,levels=unique(cor_gene$country))

pal <- rev(brewer.pal(11, "RdBu"))
levels <- pretty(c(-1, 1), 40)
col1 <- colorRampPalette(pal[1:6])
col2 <- colorRampPalette(pal[6:10])
cols <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))

cor_gene$country <- factor(cor_gene$country,
                        levels = dat1_SF_202110$country[order(dat1_SF_202110$F_Vero, decreasing = F)],ordered=TRUE)

pc1 <- ggplot(cor_gene, aes(country,gene)) +
  geom_tile(aes(fill=cor),color = "white",size=0.3) +
  #geom_text(aes(label = values), color = "black",size=3) + 
  scale_fill_gradientn(colors=cols, limits = c(-1,1),na.value="grey88") +
  #facet_wrap(~type1,scales = "free_y",labeller = "label_parsed",nrow=1)+
  scale_y_discrete(expand = c(0, 0))+
  #coord_flip()+
  xlab("") + 
  ylab("") +
  theme(panel.grid=element_blank(),
        axis.line=element_line(size=0.3,colour="black"))+
  mytheme+
  theme(axis.text.x = element_text(size = 20, color="black",vjust = 0.5,hjust = 1, angle = 90))+
  theme(legend.position = "bottom",    legend.key.size=unit(0.5,'cm'),
        legend.key.width = unit(3,"cm"),legend.key.height  = unit(0.5,"cm"))+
  guides(fill=guide_colorbar("Spearman' correlation coefficient",title.position = "top",
                             title.theme = element_text(size = 20, angle = 0,hjust = 0.5),
                             label.theme= element_text(size = 20),
                             title.vjust = 0.9)) +
  ggtitle("Correlation between monthly new reported mortality and ratio of nonsynonymous to synonymous divergence")+
  theme(plot.title = element_text(hjust = 0.5,size=24, face = "bold"))

pc1


ggsave(pc1,file=paste0(outpath0,"/FigureS16_lag2.pdf"),
       width = 20,height =10)