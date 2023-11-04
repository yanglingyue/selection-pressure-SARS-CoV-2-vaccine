# # install.packages("lubridate")
rm(list = ls())
setwd("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine")
dat<-read.csv("01_data/subsample_dnds_vaccine_diversity_low_death.csv")
dat1 <- dat

source("load_package.R")
source("plot_theme.R")
source("load_data_clean.R")
dir.create("03_output/figure_4/Supplementary_fig4_correlation/")
outpath <-  "03_output/figure_4/Supplementary_fig4_correlation/cor_202110_"

##Remove November 2021 to September 2022 as the Omicron variant emerged and circulated during this period.
my1 <- which(dat1$month=="2021-11"|dat1$month=="2021-12"|dat1$month=="2022-01"|dat1$month=="2022-02"|
              dat1$month=="2022-03"|dat1$month=="2022-04"|dat1$month=="2022-05"|dat1$month=="2022-06"|
              dat1$month=="2022-07"|dat1$month=="2022-08"|dat1$month=="2022-09") 
dat1 <-dat1[-my1,]


gn1 <- which(dat1$gene=="Nsp6"|dat1$gene=="RdRp")
dat1 <- dat1[-gn1,]


dat1$gene <- trimws(dat1$gene)
gene_type <- unique(dat1$gene)
s=length(gene_type)
print(gene_type)

colors1 <- c( '#45aaf2','#A3CB38','#1289A7','#D980FA','#F79F1F',
                       '#EE5A24','#009432','#833471','#1B1464','#0652DD',
                       '#006266','#EA2027','#ccae62','#5758BB','#cd6133')
                       

for(g in 1:s)
{
print(gene_type[g])
  
#dir.create(paste0(outpath,gene_type[g]))
#outpath_s1 <-  paste0(outpath,gene_type[g],"/")
  
  
df1= subset(dat1,dat1$gene== gene_type[g])
#dat2= subset(dat1,dat1$gene== "S")
data <- df1
df2 <- df1[,c("country","month","gene","dn_ds_mean", "new_deaths","new_deaths_per_million")]
df2 <- df2 %>%
    group_by(country) %>%
    mutate(
      lag_dn_ds_1 = lag(dn_ds_mean, 1),
      lag_death_1 = lag(new_deaths, 1),
      lag_death_rate_1 = lag(new_deaths_per_million, 1),

      lag_dn_ds_2 = lag(dn_ds_mean, 2),
      lag_death_2 = lag(new_deaths, 2),
      lag_death_rate_2 = lag(new_deaths_per_million, 2),

      lag_dn_ds_3 = lag(dn_ds_mean, 3),
      lag_death_3 = lag(new_deaths, 3),
      lag_death_rate_3 = lag(new_deaths_per_million, 3),
    )
  

  
df2$country <- as.character(df2$country)
country_name <- unique(df2$country)
print(country_name)
n =length(country_name)
  

new_data <- matrix(NA, nrow = n,ncol = 8)
new_data <- as.data.frame(new_data)
  
for(j in 1:n)
  {
    
print(country_name[j]);
Q2 = subset(df2,df2$country== country_name[j])
#Q2 <- slide(Q2, "NPI_NV", NewVar = "NPI_NV", slideBy = -1)
    
R <- matrix(NA, nrow = 1,ncol = 7)
R <- as.data.frame(R)
    
df_Q2 <- Q2[,c("dn_ds_mean","lag_death_rate_2")]
M <- as.data.frame(cor(na.omit(df_Q2),method="spearman"))
#res1 <- as.data.frame(cor.mtest(df_Q2, method = "spearman",conf.level = 0.95,exact=FALSE))
#res1 <- as.data.frame(cor.test(df_Q2$sub_lineage,df_Q2$F_I_BWP, method = "spearman",conf.level = 0.95,exact=FALSE))
res1 <- corr.test(df_Q2$dn_ds_mean, df_Q2$lag_death_rate_2, method = "spearman", ci = TRUE)
    
#R[1,1] <- M[1,2]
#R[1,2] <- res1[2,1]
#R[1,3] <- res1[2,3]
#R[1,4] <- res1[2,5]
    
R[1,1] <- res1$r
R[1,2] <- res1$p
R[1,3] <- res1$ci[1]
R[1,4] <- res1$ci[3]
    
    
fit <- lm(dn_ds_mean~lag_death_rate_2,data=df_Q2)
summary(fit)
R[1,5] <- summary(fit)$coefficients[2,1]
R[1,6] <- summary(fit)$coefficients[2,2]
R[1,7] <- summary(fit)$coefficients[2,4]
colnames(R) <- c("cor","p","low_ci","high_ci","Estimate","SE","p_lm")
    
R$country <- country_name[j]
new_data[j,] <- R
colnames(new_data) <- c("cor","p","low_ci","high_ci","Estimate","SE","p_lm","country")
    
}
  
cor_gene <- new_data
cor_gene$gene <- gene_type[g]

region <- read.csv("02_sup-data/region_code.csv")
colnames(region) <- c("iso_code","continent","country")  
cor_gene <- merge(cor_gene,region,by=c("country"),all.x = T)


dat1_S1_202110= subset(dat1,dat1$gene== "S1"&dat1$month=="2021-10")
dat1_S1_202110 <- dat1_S1_202110[,c("country","F_Vero")]
colnames(dat1_S1_202110)[2] <- "Coverage"
cor_gene <- merge(cor_gene,dat1_S1_202110,by=c("country"),all.x = T)

assign(paste("cor_gene_",gene_type[g], sep=""),cor_gene)

write.csv(cor_gene,paste0(outpath,gene_type[g],"_","lag2.csv"))

}


cor_gene_S_region <- rbind(cor_gene_S,cor_gene_S1,cor_gene_S2)


pr1 <- ggplot(cor_gene_S_region, aes(x=gene, y=cor)) + 
  geom_violin(alpha = 0.35,trim=F,width=0.45,color="black",fill="#E18727B2") +
  geom_jitter(shape = 19, position = position_jitter(0.07),alpha=1,color='#cc8e35')+
  geom_boxplot(notch=F, outlier.shape = NA,outlier.size = 2,fill="#d55e00",color="black",
               outlier.fill = "black",alpha=0.45, width = 0.15)+
  #geom_jitter(shape=19,size=2.5, position=position_jitter(0.2),color="#db3726",alpha=0.55)+
  scale_y_continuous(name="Spearman' correlation coefficient",breaks = seq(-1,1,0.5))+
  labs(x=NULL,color="" ,fill="")+
  mytheme+
  theme(axis.text.x = element_text(size = 20, color="black",vjust = 1,hjust = 0.5, angle = 0))+
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",size=1)+
  theme(axis.title.y.left = element_text(size = 20,color="black",vjust=2,face = "bold", margin = margin(0,0.4,0,0,'cm')),
        axis.title.x = element_text(size = 20, color="black",face = "bold",margin = margin(0.4,0,0,0,'cm')))+
  #scale_fill_manual(name = '',values = custom_colors) +
  #scale_color_manual(name = '',values = custom_colors) +
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+ 
  guides(color = guide_legend(override.aes = list(shape = NA, size = 3)))   

pr1

ggsave(pr1,filename = file.path(paste0("03_output/figure_4/","figure_4b_lag2.pdf")),
       width = 7.63,height =5.5)
