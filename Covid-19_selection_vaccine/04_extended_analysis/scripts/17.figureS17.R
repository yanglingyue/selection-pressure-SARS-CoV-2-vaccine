rm(list = ls())
setwd("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine/03_output/figure_4/figure_4_ccm/All_result_p_2179")
datapath <- "C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine/03_output/figure_4/figure_4_ccm/All_result_p_2179"
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

outpath <- "../../../../04_extended_analysis/plot_output/Figure_S17"
dir.create(outpath)
country_name <- trimws(unique(dat1[,"country"]))
n =length(country_name)
print(country_name)

dat_F<-read.csv("../../../../01_data/subsample_adaptation_vaccine_low.csv")
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

custom_colors <- c("#3d69ad","#88b5d6","#77488c","#d4b9da","#2f9798","#7eb780")


df_S1 <-  subset(dat1,dat1$gene=="S1")
df_S1_dnds <- df_S1[which(df_S1$Type=="dnds_xmap_death_rate"|df_S1$Type=="death_rate_xmap_dnds"),]

df_S1_dnds[which(df_S1_dnds$Type=="dnds_xmap_death_rate"),"Type"] <- "Deaths/million drives nonsynonymous to synonymous divergence ratio"
df_S1_dnds[which(df_S1_dnds$Type=="death_rate_xmap_dnds"),"Type"] <- "Nonsynonymous to synonymous divergence ratio drives deaths/million"

df_S1_dnds$Type <-factor(df_S1_dnds$Type,levels= c("Deaths/million drives nonsynonymous to synonymous divergence ratio",
                                         "Nonsynonymous to synonymous divergence ratio drives deaths/million",
                                             ordered=TRUE))  

df_S1_dnds$country <- factor(df_S1_dnds$country,
                          levels = dat1_SF_202110$country[order(dat1_SF_202110$F_Vero, decreasing = F)],ordered=TRUE)

pc_S1 <- ggplot(df_S1_dnds, aes(x=country, y=rho_surr,fill=Type)) + 
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
  scale_fill_manual(name = '',values = custom_colors[1:2]) +
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



df_S1_nonsyn <- df_S1[which(df_S1$Type=="muts_nonsyn_xmap_death_rate"|df_S1$Type=="death_rate_xmap_muts_nonsyn"),]


df_S1_nonsyn[which(df_S1_nonsyn$Type=="muts_nonsyn_xmap_death_rate"),"Type"] <- "Deaths/million drives number of nonsynonymous mutations"
df_S1_nonsyn[which(df_S1_nonsyn$Type=="death_rate_xmap_muts_nonsyn"),"Type"] <- "Number of nonsynonymous mutations drives deaths/million"

df_S1_nonsyn$Type <-factor(df_S1_nonsyn$Type,levels= c("Deaths/million drives number of nonsynonymous mutations",
                                         "Number of nonsynonymous mutations drives deaths/million",
                                             ordered=TRUE))  


df_S1_nonsyn$country <- factor(df_S1_nonsyn$country,
                        levels = dat1_SF_202110$country[order(dat1_SF_202110$F_Vero, decreasing = F)],ordered=TRUE)

pc_S2 <- ggplot(df_S1_nonsyn, aes(x=country, y=rho_surr,fill=Type)) + 
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
  scale_fill_manual(name = '',values = custom_colors[3:4]) +
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



df_S1_fitness <- df_S1[which(df_S1$Type=="fitness_xmap_death_rate"|df_S1$Type=="death_rate_xmap_fitness"),]

df_S1_fitness[which(df_S1_fitness$Type=="fitness_xmap_death_rate"),"Type"] <- "Deaths/million drives mutational fitness"
df_S1_fitness[which(df_S1_fitness$Type=="death_rate_xmap_fitness"),"Type"] <- "Mutational fitness forces deaths/million"

df_S1_fitness$Type <-factor(df_S1_fitness$Type,levels= c("Deaths/million drives mutational fitness",
                                       "Mutational fitness forces deaths/million",
                                         ordered=TRUE))  

df_S1_fitness$country <- factor(df_S1_fitness$country,
                        levels = dat1_SF_202110$country[order(dat1_SF_202110$F_Vero, decreasing = F)],ordered=TRUE)

pc_S3 <- ggplot(df_S1_fitness, aes(x=country, y=rho_surr,fill=Type)) + 
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
  scale_fill_manual(name = '',values = custom_colors[5:6]) +
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

pc_S3



p<-ggarrange(pc_S1, pc_S2,pc_S3,ncol =1, nrow = 3,widths = c(1,1,1),heights = c(1,1,1),
             labels = c("A","B","C"),  vjust=1.1,align = "hv",
             common.legend=TRUE,legend = "bottom",
             font.label = list(size = 24, face = "bold"))
p

ggsave(p,filename = file.path(outpath,"/","figure_S17.pdf"),
       width = 20,height =18)

