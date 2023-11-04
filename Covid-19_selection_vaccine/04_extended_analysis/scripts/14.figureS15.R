rm(list = ls())
setwd("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine/03_output/figure_4/ccm_2179_death_E2_4/All_result_p_2179")
datapath <- "C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine/03_output/figure_4/ccm_2179_death_E2_4/All_result_p_2179"
datafile <- list.files(datapath, pattern = "*.csv$", full.names = TRUE)
files=list()
read.txt <- function(x) {
  des <- readLines(x)
  return(paste(des, collapse = ""))
}
review <- lapply(datafile, read.csv)
docname <- list.files(datapath, pattern = "*.csv$")
myfiles = do.call(rbind, lapply(docname, function(x) read.csv(x, stringsAsFactors = FALSE)))
data <- myfiles
dat1 <- data

source("../../../../load_package.R")
source("../../../../plot_theme.R")

outpath <- "../../../../04_extended_analysis/plot_output/Figure_S15"
dir.create(outpath)
country_name <- trimws(unique(dat1[,"country"]))
n =length(country_name)
print(country_name)

dat_F<-read.csv("../../../../01_data/subsample_dnds_vaccine_diversity_low_death.csv")
dat1_SF_202110 <- subset(dat_F,dat_F$gene== "S2"&dat_F$month=="2021-10")





type_name <- trimws(unique(dat1[,"Type"]))
#t <- length(type_name)
print(type_name)
dat1$Type <- as.character(dat1$Type)
gene_type <- unique(dat1$gene)
s=length(gene_type)
print(gene_type)

dat1[which(dat1$rho_actual<=0),"rho_actual"] <- 0
dat1[which(dat1$rho_surr<=0),"rho_surr"] <- 0

custom_colors <- c("#3d69ad","#B09C85B2","#fdbe67","#b153a0","#4e69b1","#949698","#b12425","#009e73","#d55e00","#99d594","#77488c","#a86ca0","#229f99","#377eb8","#e41a1c")


#df_gene[which(df_gene$Type=="dnds_xmap_death"),"Type"] <- "Deaths forces dN/dS ratio"
#df_gene[which(df_gene$Type=="death_xmap_dnds"),"Type"] <- "dN/dS ratio forces deaths"

#df_gene[which(df_gene$Type=="dnds_xmap_death_rate"),"Type"] <- "Deaths forces dN/dS ratio"
#df_gene[which(df_gene$Type=="death_rate_xmap_dnds"),"Type"] <- "dN/dS ratio forces deaths"
df_S1 <-  subset(dat1,dat1$gene=="S1")
df_S1[which(df_S1$Type=="dnds_xmap_death_rate"),"Type"] <- "Mortality drives ratio of nonsynonymous to synonymous divergence (S1)"
df_S1[which(df_S1$Type=="death_rate_xmap_dnds"),"Type"] <- "Ratio of nonsynonymous to synonymous divergence drives mortality (S1)"

df_S1$Type <-factor(df_S1$Type,levels= c("Mortality drives ratio of nonsynonymous to synonymous divergence (S1)",
                                         "Ratio of nonsynonymous to synonymous divergence drives mortality (S1)",
                                             ordered=TRUE))  

df_S1 <- df_S1[-which(is.na(df_S1$Type)),]

df_S1$country <- factor(df_S1$country,
                          levels = dat1_SF_202110$country[order(dat1_SF_202110$F_Vero, decreasing = F)],ordered=TRUE)

pc_S1 <- ggplot(df_S1, aes(x=country, y=rho_surr,fill=Type)) + 
  geom_boxplot(notch=F, outlier.shape = NA,outlier.size = 1,
               outlier.fill = "black")+  
  geom_point(aes(x=country,y=rho_actual,shape=factor(sign1),color=factor(sign1)),size=3,stroke = 0.8)+
  #geom_point(position=position_jitterdodge(jitter.width=2, dodge.width = 0), 
  #           pch=21, aes(fill=factor(Type),color=factor(Type)), show.legend = F,size=0.5)+
  #coord_flip(ylim = c(0, 1))+
  facet_wrap(~Type,scales = "free_y",ncol=1)+
  coord_cartesian(ylim = c(0, 1))+
  scale_y_continuous(name=expression(paste("Cross Mapping Skill (",rho,")")), 
                     breaks = seq(0,1,0.5))+
  labs(x=NULL,color="" ,fill="")+
  mytheme+
  scale_fill_manual(name = '',values = custom_colors) +
  #scale_color_manual(name = '',values = custom_colors) +
  scale_color_manual(name = '', values = c("gray", "black")) +  # Set colors for p < 0.05 and p > 0.05
  scale_shape_manual(name = '',values = c(1, 19),label=c("p > 0.05","p < 0.05")) +
  theme(panel.spacing = unit(1, "cm"))+
  theme(strip.text = element_text(face="bold",size=20))+
  #ggtitle(gene_type[g])+
  theme(axis.text.x = element_text(size = 20, color="black",vjust = 0.5,hjust = 1, angle = 90))+
  theme(axis.title.y.left = element_text(size = 24,color="black",vjust=2,face = "bold", margin = margin(0,0.4,0,0,'cm')),
        axis.title.x = element_text(size = 24, color="black",face = "bold",margin = margin(0.4,0,0,0,'cm')))+
  #guides(shape = FALSE,color = guide_legend(override.aes = list(shape = NA, size = 3)))   
  guides(color = guide_legend(override.aes = list(shape = NA, size = 3))) +
  guides(fill = "none",color="none",shape="none")

#guides(color = guide_legend(override.aes = list(shape = c(NA, NA), size = c(3, 3), fill = c("gray", "black"))))   

pc_S1 





df_S2 <-  subset(dat1,dat1$gene=="S2")
df_S2[which(df_S2$Type=="dnds_xmap_death_rate"),"Type"] <- "Mortality drives ratio of nonsynonymous to synonymous divergence (S2)"
df_S2[which(df_S2$Type=="death_rate_xmap_dnds"),"Type"] <- "Ratio of nonsynonymous to synonymous divergence drives mortality (S2)"

df_S2$Type <-factor(df_S2$Type,levels= c("Mortality drives ratio of nonsynonymous to synonymous divergence (S2)",
                                         "Ratio of nonsynonymous to synonymous divergence drives mortality (S2)",
                                             ordered=TRUE))  

df_S2 <- df_S2[-which(is.na(df_S2$Type)),]

df_S2$country <- factor(df_S2$country,
                        levels = dat1_SF_202110$country[order(dat1_SF_202110$F_Vero, decreasing = F)],ordered=TRUE)

pc_S2 <- ggplot(df_S2, aes(x=country, y=rho_surr,fill=Type)) + 
  geom_boxplot(notch=F, outlier.shape = NA,outlier.size = 1,
               outlier.fill = "black")+  
  geom_point(aes(x=country,y=rho_actual,shape=factor(sign1),color=factor(sign1)),size=3,stroke = 0.8)+
  #geom_point(position=position_jitterdodge(jitter.width=2, dodge.width = 0), 
  #           pch=21, aes(fill=factor(Type),color=factor(Type)), show.legend = F,size=0.5)+
  #coord_flip(ylim = c(0, 1))+
  facet_wrap(~Type,scales = "free_y",ncol=1)+
  coord_cartesian(ylim = c(0, 1))+
  scale_y_continuous(name=expression(paste("Cross Mapping Skill (",rho,")")), 
                     breaks = seq(0,1,0.5))+
  labs(x=NULL,color="" ,fill="")+
  mytheme+
  scale_fill_manual(name = '',values = custom_colors) +
  #scale_color_manual(name = '',values = custom_colors) +
  scale_color_manual(name = '', values = c("gray", "black")) +  # Set colors for p < 0.05 and p > 0.05
  scale_shape_manual(name = '',values = c(1, 19),label=c("p > 0.05","p < 0.05")) +
  theme(panel.spacing = unit(1, "cm"))+
  theme(strip.text = element_text(face="bold",size=20))+
  #ggtitle(gene_type[g])+
  theme(axis.text.x = element_text(size = 20, color="black",vjust = 0.5,hjust = 1, angle = 90))+
  theme(axis.title.y.left = element_text(size = 24,color="black",vjust=2,face = "bold", margin = margin(0,0.4,0,0,'cm')),
        axis.title.x = element_text(size = 24, color="black",face = "bold",margin = margin(0.4,0,0,0,'cm')))+
  #guides(shape = FALSE,color = guide_legend(override.aes = list(shape = NA, size = 3)))   
  guides(color = guide_legend(override.aes = list(shape = NA, size = 3)))   
#guides(color = guide_legend(override.aes = list(shape = c(NA, NA), size = c(3, 3), fill = c("gray", "black"))))   

pc_S2




df_S <-  subset(dat1,dat1$gene=="S")
df_S[which(df_S$Type=="dnds_xmap_death_rate"),"Type"] <- "Mortality drives ratio of nonsynonymous to synonymous divergence (S)"
df_S[which(df_S$Type=="death_rate_xmap_dnds"),"Type"] <- "Ratio of nonsynonymous to synonymous divergence drives mortality (S)"

df_S$Type <-factor(df_S$Type,levels= c("Mortality drives ratio of nonsynonymous to synonymous divergence (S)",
                                       "Ratio of nonsynonymous to synonymous divergence drives mortality (S)",
                                         ordered=TRUE))  

df_S <- df_S[-which(is.na(df_S$Type)),]
df_S$country <- factor(df_S$country,
                        levels = dat1_SF_202110$country[order(dat1_SF_202110$F_Vero, decreasing = F)],ordered=TRUE)

pc_S <- ggplot(df_S, aes(x=country, y=rho_surr,fill=Type)) + 
  geom_boxplot(notch=F, outlier.shape = NA,outlier.size = 1,
               outlier.fill = "black")+  
  geom_point(aes(x=country,y=rho_actual,shape=factor(sign1),color=factor(sign1)),size=3,stroke = 0.8)+
  #geom_point(position=position_jitterdodge(jitter.width=2, dodge.width = 0), 
  #           pch=21, aes(fill=factor(Type),color=factor(Type)), show.legend = F,size=0.5)+
  #coord_flip(ylim = c(0, 1))+
  facet_wrap(~Type,scales = "free_y",ncol=1)+
  coord_cartesian(ylim = c(0, 1))+
  scale_y_continuous(name=expression(paste("Cross Mapping Skill (",rho,")")), 
                     breaks = seq(0,1,0.5))+
  labs(x=NULL,color="" ,fill="")+
  mytheme+
  scale_fill_manual(name = '',values = custom_colors) +
  #scale_color_manual(name = '',values = custom_colors) +
  scale_color_manual(name = '', values = c("gray", "black")) +  # Set colors for p < 0.05 and p > 0.05
  scale_shape_manual(name = '',values = c(1, 19),label=c("p > 0.05","p < 0.05")) +
  theme(panel.spacing = unit(1, "cm"))+
  theme(strip.text = element_text(face="bold",size=20))+
  #ggtitle(gene_type[g])+
  theme(axis.text.x = element_text(size = 20, color="black",vjust = 0.5,hjust = 1, angle = 90))+
  theme(axis.title.y.left = element_text(size = 24,color="black",vjust=2,face = "bold", margin = margin(0,0.4,0,0,'cm')),
        axis.title.x = element_text(size = 24, color="black",face = "bold",margin = margin(0.4,0,0,0,'cm')))+
  #guides(shape = FALSE,color = guide_legend(override.aes = list(shape = NA, size = 3)))   
  guides(color = guide_legend(override.aes = list(shape = NA, size = 3)))+
  guides(fill = "none",color="none",shape="none")
  
#guides(color = guide_legend(override.aes = list(shape = c(NA, NA), size = c(3, 3), fill = c("gray", "black"))))   

pc_S



p<-ggarrange(pc_S, pc_S1,pc_S2,ncol =1, nrow = 3,widths = c(1,1,1),heights = c(1,1,1),
             labels = c("A","B","C"),  vjust=1.1,align = "hv",
             common.legend=TRUE,legend = "bottom",
             font.label = list(size = 24, face = "bold"))
p

ggsave(p,filename = file.path(outpath,"/","figure_S15.pdf"),
       width = 20,height =18)

