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

outpath <- "../../"

country_name <- trimws(unique(dat1[,"country"]))
n =length(country_name)
print(country_name)

type_name <- trimws(unique(dat1[,"Type"]))
#t <- length(type_name)
print(type_name)
dat1$Type <- as.character(dat1$Type)
gene_type <- unique(dat1$gene)
s=length(gene_type)
print(gene_type)

res_type <- matrix(NA,  nrow = length(type_name), ncol=ncol(dat1))
res_type <- as.data.frame(res_type)

res_gene <- matrix(NA,  nrow = length(type_name)*length(gene_type), ncol=ncol(dat1))
res_gene <- as.data.frame(res_gene)

res_country <- matrix(NA,  nrow = length(type_name)*length(gene_type)*length(country_name), ncol=ncol(dat1))
res_country <- as.data.frame(res_country)

for(co1 in 1:length(country_name))
{
  print(country_name[co1])
  df_country <-   subset(dat1,dat1$country== country_name[co1])
  for(g in 1:s)
  {
    print(gene_type[g])
    df_gene= subset(df_country,df_country$gene== gene_type[g])
  for(t in 1:length(type_name))
    {
    print(type_name[t])
    df_type <-   subset(df_gene,df_gene$Type== type_name[t])
    res_type[t,] <- df_type[1,]
    }
    res_gene[(c(length(type_name)*g-c(length(type_name)-1)):(length(type_name)*g)),] <- res_type[,]
  }
  res_country[(c(length(type_name)*s*co1-c(length(type_name)*s-1)):(length(type_name)*s*co1)),] <- res_gene[,]
}

colnames(res_country) <- colnames(dat1)
res_country1 <- res_country[,c("Type","rho_actual","kend_tau","kend_p","sign","sign1","gene" ,"country","seed")]

custom_colors <- c("#3d69ad","#B09C85B2","#377eb8","#e41a1c","#009e73","#d55e00","#99d594","#77488c","#a86ca0","#229f99")

#res_country1 <- res_country
res_country1[which(res_country1$rho_actual<=0),"rho_actual"] <- 0

res_country1$type1 <- res_country1$Type

res_country1[which(res_country1$Type=="dnds_xmap_death_rate"),"Type"] <- "Mortality forces dN/dS ratio"
res_country1[which(res_country1$Type=="death_rate_xmap_dnds"),"Type"] <- "dN/dS ratio forces mortality"


#res_country1[which(res_country1$Type=="dnds_xmap_death_rate"),"Type"] <- "Death rate forces dN/dS ratio"
#res_country1[which(res_country1$Type=="death_rate_xmap_dnds"),"Type"] <- "dN/dS ratio forces death rate"

res_country1$Type <-factor(res_country1$Type,
                           levels= c("Mortality forces dN/dS ratio", "dN/dS ratio forces mortality",ordered=TRUE)) 
 
res_country2 <- res_country1[-which(is.na(res_country1$Type)),]
res_country2$gene <-factor(res_country2$gene,levels= c('S','S1','S2','E','M','N','ORF1a','ORF1b','ORF3a',
                                                       'ORF6','ORF7a','ORF7b','ORF8','ORF9b',ordered=TRUE)) 
res_country2 <- res_country2[-which(is.na(res_country2$gene)),]

#devtools::install_github("psyteachr/introdataviz")
res_country3 <- res_country2

res_country3$gene <-factor(res_country3$gene,levels= c('S','S1','S2',ordered=TRUE)) 
res_country3 <- res_country3[-which(is.na(res_country3$gene)),]

pr1 <- ggplot(res_country3, aes(x=gene, y=rho_actual,fill=Type)) + 
  geom_boxplot(notch=F, outlier.shape = 19,outlier.size = 1,color="black",
               outlier.fill = "black",alpha=0.65, width = 0.5)+
  stat_compare_means(test = 'wilcox.test', label = 'p.format',show.legend = T)+
  #stat_compare_means(test = 'wilcox.test', label = 'p.signif',show.legend = T)+
  scale_x_discrete(labels=c("Mortality forces dN/dS ratio"="Mortality forces\n dN/dS ratio",
                            "dN/dS ratio forces mortality"="dN/dS ratio\n forces mortality"))+
  #coord_flip(ylim = c(0, 1))+
  coord_cartesian(ylim = c(0, 1))+
  scale_y_continuous(name=expression(paste("Cross Mapping Skill (",rho,")")), 
                     breaks = seq(0,1,0.2))+
  
  #scale_y_continuous(name="Cross Mapping Skill (p)", breaks = seq(0,1,0.2))+
  labs(x=NULL,color="" ,fill="")+
  mytheme+
  theme(axis.text.x = element_text(size = 20, color="black",vjust = 1,hjust = 0.5, angle = 0))+
  #theme(plot.margin=unit(c(3,0.5,2.5,0.5),'lines'))+
  theme(axis.title.y.left = element_text(size = 20,color="black",vjust=2,face = "bold", margin = margin(0,0.4,0,0,'cm')),
        axis.title.x = element_text(size = 20, color="black",face = "bold",margin = margin(0.4,0,0,0,'cm')))+
  scale_fill_manual(name = '',values = custom_colors) +
  scale_color_manual(name = '',values = custom_colors) +
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+ 
  guides(fill = "none",color="none")
  #guides(color = guide_legend(override.aes = list(shape = NA, size = 3)),
  #       fill=guide_legend(nrow=2,byrow = T))   

pr1

ggsave(pr1,filename = file.path(paste0(outpath,"figure_4a.pdf")),
       width = 7.5,height =5.5)





fisher_p_value <- matrix(NA, ncol= length(type_name)+1,  nrow =s)
fisher_p_value <- as.data.frame(fisher_p_value)
colnames(fisher_p_value) <- c("gene","fisher_dnds_death","fisher_death_dnds",
                              "fisher_dnds_death_rate","fisher_death_rate_dnds")



for(g in 1:s)
  
{
  print(gene_type[g])
  
  res_all <- res_country
  res_all_gene<- res_all[res_all$gene %in% gene_type[g], ]
 
  if (length(which(is.na(res_all_gene$sign))) > 0) {
    res_all_gene6 <- res_all_gene[-which(is.na(res_all_gene$sign)),]
  } else {
    
    res_all_gene6 <- res_all_gene
    
  }
  
  df_dnds_death <- res_all_gene6[which(res_all_gene6$Type=="dnds_xmap_death"),]
  df_death_dnds <- res_all_gene6[which(res_all_gene6$Type=="death_xmap_dnds"),]
  df_dnds_death_rate <- res_all_gene6[which(res_all_gene6$Type=="dnds_xmap_death_rate"),]
  df_death_rate_dnds <- res_all_gene6[which(res_all_gene6$Type=="death_rate_xmap_dnds"),]
  
  
  fisher_p_value[g,"gene"] <- gene_type[g]
  fisher_p_value[g,"fisher_dnds_death"] <- fisher(df_dnds_death$sign)[1]
  fisher_p_value[g,"fisher_death_dnds"] <-fisher(df_death_dnds$sign)[1]
  fisher_p_value[g,"fisher_dnds_death_rate"] <-fisher(df_dnds_death_rate$sign)[1]
  fisher_p_value[g,"fisher_death_rate_dnds"] <-fisher(df_death_rate_dnds$sign)[1]
  
  
  write.csv(fisher_p_value,file=paste0(outpath,"fisher_method_p_mortality_dnds",".csv",sep=""),row.names =F)
  
  
}
