# # install.packages("lubridate")
rm(list = ls())
setwd("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine")
dat<-read.csv("01_data/subsample_adaptation_vaccine_low.csv")
dat1 <- dat

source("load_package.R")
source("plot_theme.R")
source("load_data_clean.R")
dir.create("03_output/figure_4/fig4_cor/")
outpath <-  "03_output/figure_4/fig4_cor/cor_202110_"

##Remove November 2021 to September 2022 as the Omicron variant emerged and circulated during this period.
my1 <- which(dat1$month=="2021-11"|dat1$month=="2021-12"|dat1$month=="2022-01"|dat1$month=="2022-02"|
              dat1$month=="2022-03"|dat1$month=="2022-04"|dat1$month=="2022-05"|dat1$month=="2022-06"|
              dat1$month=="2022-07"|dat1$month=="2022-08"|dat1$month=="2022-09") 
dat1 <-dat1[-my1,]


gn1 <- which(dat1$gene=="S1")
dat1 <- dat1[gn1,]


dat1$gene <- trimws(dat1$gene)
gene_type <- unique(dat1$gene)
s=length(gene_type)
print(gene_type)

colnames(dat1)[4] <- "dnds"
colnames(dat1)[6] <- "muts_nonsyn"
colnames(dat1)[10] <- "Mutational fitness"

var_type <- colnames(dat1)[c(4,6,10)]

dat1$country <- as.character(dat1$country)
country_name <- unique(dat1$country)
print(country_name)
n =length(country_name)

for(g in 1:s)
{
print(gene_type[g])
  
#dir.create(paste0(outpath,gene_type[g]))
#outpath_s1 <-  paste0(outpath,gene_type[g],"/")

df1= subset(dat1,dat1$gene== gene_type[g])
#dat2= subset(dat1,dat1$gene== "S")
data <- df1

final_data <- matrix(NA, nrow = n*length(var_type),ncol = 6)
final_data <- as.data.frame(final_data)


## Calculate the Spearman correlation between three quantified measures of SARS-CoV-2 adaptation for each country 
# and the reported monthly deaths per million people.

for (t in 1:length(var_type))
{

print(var_type[t])
  
df2 <- df1[,c("country","month","gene",var_type[t], "new_deaths_per_million")]

df2 <- df2 %>%
  group_by(country) %>%
  mutate(
    lag_death_rate_1 = lag(new_deaths_per_million, 1),
    
    lag_death_rate_2 = lag(new_deaths_per_million, 2),
    
    lag_death_rate_3 = lag(new_deaths_per_million, 3),
  )


new_data <- matrix(NA, nrow = n,ncol = 6)
new_data <- as.data.frame(new_data)

for(j in c(1:n))
{
  
  print(country_name[j]);
  Q2 = subset(df2,df2$country== country_name[j])
  #Q2 <- slide(Q2, "NPI_NV", NewVar = "NPI_NV", slideBy = -1)
  
  R <- matrix(NA, nrow = 1,ncol = 4)
  R <- as.data.frame(R)
  
  df_Q2 <- Q2[,c(var_type[t],"lag_death_rate_2")]
  M <- as.data.frame(cor(na.omit(df_Q2),method="spearman"))
  res1 <- corr.test(df_Q2[,var_type[t]], df_Q2$lag_death_rate_2, method = "spearman", ci = TRUE)
  
  
  R[1,1] <- res1$r
  R[1,2] <- res1$p
  R[1,3] <- res1$ci[1]
  R[1,4] <- res1$ci[3]
  

  colnames(R) <- c("cor","p","low_ci","high_ci")
  
  R$country <- country_name[j]
  R$vars_type <- var_type[t]
  new_data[j,] <- R
  colnames(new_data) <- colnames(R)
  
}

final_data[(c(n*t-c(n-1)):(n*t)),] <- new_data[,]
colnames(final_data) <- colnames(new_data)
cor_gene <- final_data
cor_gene$gene <- gene_type[g]
region <- read.csv("02_sup-data/region_code.csv")
colnames(region) <- c("iso_code","continent","country")  
cor_gene <- merge(cor_gene,region,by=c("country"),all.x = T)


dat1_S1_202110= subset(dat1,dat1$gene== "S1"&dat1$month=="2021-10")
dat1_S1_202110 <- dat1_S1_202110[,c("country","F_Vero")]
colnames(dat1_S1_202110)[2] <- "Coverage"
cor_gene <- merge(cor_gene,dat1_S1_202110,by=c("country"),all.x = T)

#assign(paste("cor_gene_",gene_type[g], sep=""),cor_gene)
}

write.csv(cor_gene,paste0(outpath,gene_type[g],"_","lag2.csv"),row.names = F)

}


#cor_gene_S_region <- rbind(cor_gene_S,cor_gene_S1,cor_gene_S2)


####Plot Figure 4B#######

cor_gene$vars_type <-factor(cor_gene$vars_type,
                               levels= c("dnds", "muts_nonsyn",
                                         "Mutational fitness",ordered=TRUE)) 
pr1 <- ggplot(cor_gene, aes(x=vars_type, y=cor)) + 
  geom_violin(alpha = 0.35,trim=F,width=0.45,color="black",fill="#E18727B2") +
  geom_jitter(shape = 19, position = position_jitter(0.07),alpha=1,color='#cc8e35')+
  geom_boxplot(notch=F, outlier.shape = NA,outlier.size = 2,fill="#d55e00",color="black",
               outlier.fill = "black",alpha=0.45, width = 0.15)+
  #geom_jitter(shape=19,size=2.5, position=position_jitter(0.2),color="#db3726",alpha=0.55)+
  scale_x_discrete(labels=c("dnds"="Nonsynonymous to \nsynonymous divergence ratio",
                            "muts_nonsyn"="Number of nonsynonymous \nmutations",
                            "Mutational fitness"="Mutational fitness"))+
  scale_y_continuous(name="Spearman' correlation coefficient",breaks = seq(-1,1,0.5))+
  labs(x=NULL,color="" ,fill="")+
  mytheme+
  theme(axis.text.x = element_text(size = 20, color="black",vjust = 1,hjust = 0.5, angle = 0))+
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",size=1)+
  theme(axis.title.y.left = element_text(size = 20,color="black",vjust=2,face = "bold", margin = margin(0,0.4,0,0,'cm')),
        axis.title.x = element_text(size = 20, color="black",face = "bold",margin = margin(0.4,0,0,0,'cm')))+
  #scale_fill_manual(name = '',values = custom_colors) +
  #scale_color_manual(name = '',values = custom_colors) +
  theme(axis.text.x = element_text(angle = 30, hjust = 0.75,vjust = 0.75))+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+ 
  guides(color = guide_legend(override.aes = list(shape = NA, size = 3)))   

pr1

ggsave(pr1,filename = file.path(paste0("03_output/figure_4/","figure_4B.pdf")),
       width = 8.6,height =7.5)




####Plot Figure 4C#######
dat1_SF_202110 <- subset(dat1,dat1$gene== "S1"&dat1$month=="2021-10")

cor_gene<- cor_gene[order(cor_gene$Coverage,cor_gene$country,decreasing=F),]
cor_gene$country <-factor(cor_gene$country,levels=unique(cor_gene$country))

pal <- rev(brewer.pal(11, "RdBu"))
levels <- pretty(c(-1, 1), 40)
col1 <- colorRampPalette(pal[2:6])
col2 <- colorRampPalette(pal[6:10])


cols <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))

cor_gene$country <- factor(cor_gene$country,
                           levels = dat1_SF_202110$country[order(dat1_SF_202110$F_Vero, decreasing = F)],ordered=TRUE)

pc1 <- ggplot(cor_gene, aes(vars_type,country)) +
  geom_tile(aes(fill=cor),color = "white",size=0.3) +
  geom_text(aes(label = ifelse(p < 0.05, "×", "")), color = "black", size = 7) +  # Add this line
  #geom_text(aes(label = values), color = "black",size=3) + 
  scale_fill_gradientn(colors=cols, limits = c(-1,1),na.value="grey88") +
  #facet_wrap(~type1,scales = "free_y",labeller = "label_parsed",nrow=1)+
  scale_x_discrete(labels=c("dnds"="Nonsynonymous to \nsynonymous divergence ratio",
                            "muts_nonsyn"="Number of nonsynonymous \nmutations",
                            "Mutational fitness"="Mutational fitness"),expand = c(0, 0))+
  #coord_flip()+
  xlab("") + 
  ylab("") +
  theme(panel.grid=element_blank(),
        axis.line=element_line(size=0.3,colour="black"))+
  mytheme+
  theme(axis.text.x = element_text(size = 20, color="black",vjust = 0.5,hjust = 1, angle = 90))+
  theme(axis.text.x = element_text(angle = 30, hjust = 0.75,vjust = 0.75))+
  theme(legend.position = "right",    legend.key.size=unit(0.5,'cm'),
        legend.key.width = unit(0.5,"cm"),legend.key.height  = unit(4,"cm"))+
  guides(fill=guide_colorbar("Spearman' correlation coefficient",title.position = "top",
                             title.theme = element_text(size = 20, angle = 90,hjust = 0.5),
                             label.theme= element_text(size = 20),
                             title.vjust = 0.9)) +
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5,size=24, face = "bold"))

pc1


ggsave(pc1,file=paste0("03_output/figure_4/","figure_4C.pdf"),
       width = 10,height =18)