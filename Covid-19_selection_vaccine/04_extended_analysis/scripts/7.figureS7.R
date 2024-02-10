setwd("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine")
# Load packages (assumed already installed)
rm(list = ls())
library(data.table);library(corrplot);library(ggcorrplot);library(RColorBrewer)
library(ggpubr);


source("load_package.R")
source("plot_theme.R")

dir.create("04_extended_analysis/plot_output/Figure_S7")
outpath <-  "04_extended_analysis/plot_output/Figure_S7/"

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




df_prop <- read_csv("02_sup-data/lineage_prop_low_200_globe_mean.csv")

df_prop$type <-factor(df_prop$type,levels=c("20I (Alpha, V1)","20H (Beta, V2)","20J (Gamma, V3)",
                                      "21A (Delta)","21I (Delta)","21J (Delta)",
                                      "21M (Omicron)","21K (Omicron)","21L (Omicron)",
                                      "22A (Omicron)","22B (Omicron)","22C (Omicron)",
                                      "22D (Omicron)","22E (Omicron)","22F (Omicron)",
                                      "Other"))

p1 <- ggplot(df_prop,aes(month,prop,fill=type))+
  geom_bar(stat="identity",position = position_stack(),size =0.15,color = "black",alpha=0.75)+
  scale_y_continuous(breaks = seq(0,1,0.25),
                     labels = seq(0,100,25),expand = c(0, 0))+
  labs(x="", y= "Proportion (%)")+
  scale_x_discrete(breaks=c("2020-03","2020-06","2020-09","2020-12","2021-03","2021-06","2021-09","2021-12",
                            "2022-03","2022-06","2022-09"),
                   labels=c("Mar\n2020","Jun\n2020","Sep\n2020","Dec\n2020","Mar\n2021","Jun\n2021","Sep\n2021","Dec\n2021",
                            "Mar\n2022","Jun\n2022","Sep\n2022"))+
  scale_fill_manual(labels =c(" 20I (Alpha, V1)     "," 20H (Beta, V2)     "," 20J (Gamma, V3)     ",
                              " 21A (Delta)     "," 21I (Delta)     "," 21J (Delta)     ",
                              " 21M (Omicron)     "," 21K (Omicron)     "," 21L (Omicron)    ",
                              " 22A (Omicron)     "," 22B (Omicron)     "," 22C (Omicron)     ",
                              " 22D (Omicron)     "," 22E (Omicron)     "," 22F (Omicron)     ",
                              " Other     ") ,
                    values =custom_colors$discrete)+
  labs(fill=NULL)+
  guides(fill=guide_legend(nrow=4,byrow = T))+
  mytheme+
  theme(axis.title.y.left = element_text(size = 18,color="black",vjust=2,face = "bold", margin = margin(0,0.4,0,0,'cm')),
        axis.title.x = element_text(size = 18, color="black",face = "bold",margin = margin(0.4,0,0,0,'cm')))+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+
  theme(axis.ticks.x=element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(legend.spacing.x = unit(0.05,"cm"))
p1

ggsave(p1,filename = file.path(outpath,"figS7.pdf"),
       width = 12,height =6)

