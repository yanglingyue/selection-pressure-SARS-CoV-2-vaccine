setwd("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine")
# Load packages (assumed already installed)
rm(list = ls())
library(data.table);library(corrplot);library(ggcorrplot);library(RColorBrewer)
library(ggpubr);

dat<-read.csv("01_data/subsample_dnds_vaccine_diversity_low_death.csv")
dat1 <- dat

source("load_package.R")
source("plot_theme.R")
source("load_data_clean.R")

dir.create("04_extended_analysis/plot_output/Figure_S2")
outpath <-  "04_extended_analysis/plot_output/Figure_S2/"


colors1 <- c('#45aaf2','#A3CB38','#1289A7','#D980FA','#F79F1F',
             '#EE5A24','#009432','#833471','#1B1464','#0652DD',
             '#006266','#EA2027','#1B1464','#5758BB','#6F1E51')

dat1$gene <- trimws(dat1$gene)
gene_type <- unique(dat1$gene)
s=length(gene_type)
print(gene_type)

gene_dat <- matrix(NA, nrow = s*31,ncol =8)
gene_dat <- as.data.frame(gene_dat)

for(g in 1:s)
  
{
  
  G2= subset(dat1,dat1$gene== gene_type[g])
  
  G3=G2[,c("country","month","gene","FW_BW_P","F_Vero","travel_route","NPI_NV","new_cases_month","dn_ds_mean")]
  
  var <- c("month","country","dn_ds_mean")
  alldata <- G3
  month_name <- unique(alldata$month)
  print(month_name)
  m =length(month_name)
  
  ourdata <-alldata[,var] 
  colnames(ourdata) <- c("month","country","value")
  
  res_ci <- matrix(NA,m,5)
  res_ci <- as.data.frame(res_ci)
  
  for (j in 1:m)
  {
    print(month_name[j]);
    result<- matrix(NA,1,3)
    result <- as.data.frame(result) 
    Q2 = subset(ourdata,trimws(ourdata$month)== month_name[j])
    # Calculate the mean of the sample data
    mean_value <- mean(Q2$value,na.rm = T)
    Q3 <- na.omit(Q2)
    # Compute the size
    n <- length(Q3$value)
    
    # Find the standard deviation
    standard_deviation <- sd(Q3$value)
    
    # Find the standard error
    standard_error <- standard_deviation / sqrt(n)
    alpha = 0.05
    degrees_of_freedom = n - 1
    t_score = qt(p=alpha/2, df=degrees_of_freedom,lower.tail=F)
    margin_error <- t_score * standard_error
    # Calculating lower bound and upper bound
    lower_bound <- mean_value - margin_error
    upper_bound <- mean_value + margin_error
    
    # Print the confidence interval
    print(c(lower_bound,mean_value,upper_bound))
    
    result[,] <- c(lower_bound,mean_value,upper_bound)
    colnames(result) <- c("dnds_lci","dnds_mean","dnds_uci")
    result$month <- month_name[j]
    result$gene <- gene_type[g]
    #result$Method <-method_type[q]
    
    res_ci[j,] <- result
  }
  colnames(res_ci) <- colnames(result)
  
  
  var <- c("month","country","F_Vero")
  alldata <- G3
  
  ourdata <-alldata[,var] 
  colnames(ourdata) <- c("month","country","value")
  vero_ci <- matrix(NA,m,5)
  vero_ci <- as.data.frame(vero_ci)
  
  for (j in 1:m)
  {
    print(month_name[j]);
    result<- matrix(NA,1,3)
    result <- as.data.frame(result) 
    Q2 = subset(ourdata,trimws(ourdata$month)== month_name[j])
    # Calculate the mean of the sample data
    mean_value <- mean(Q2$value,na.rm = T)
    Q3 <- na.omit(Q2)
    # Compute the size
    n <- length(Q3$value)
    
    # Find the standard deviation
    standard_deviation <- sd(Q3$value)
    
    # Find the standard error
    standard_error <- standard_deviation / sqrt(n)
    alpha = 0.05
    degrees_of_freedom = n - 1
    t_score = qt(p=alpha/2, df=degrees_of_freedom,lower.tail=F)
    margin_error <- t_score * standard_error
    # Calculating lower bound and upper bound
    lower_bound <- mean_value - margin_error
    upper_bound <- mean_value + margin_error
    
    # Print the confidence interval
    print(c(lower_bound,mean_value,upper_bound))
    
    result[,] <- c(lower_bound,mean_value,upper_bound)
    colnames(result) <- c("vero_lci","vero_mean","vero_uci")
    result$month <- month_name[j]
    result$gene <- gene_type[g]
    #result$Method <-method_type[q]
    
    vero_ci[j,] <- result
  }
  colnames(vero_ci) <- colnames(result)
  
  
  D1 <- merge(vero_ci,res_ci,by=c("gene","month"),all.y = T)
  
  
  gene_dat[(c(g*31-c(30)):(g*31)),] <- D1[,]
  colnames(gene_dat) <- colnames(D1)
}

dat2 <- gene_dat

dat3= subset(dat2,dat2$gene== "E"|dat2$gene== "M"|dat2$gene== "N"|dat2$gene== "ORF1a"|
               dat2$gene== "ORF1b"|dat2$gene== "ORF3a"|dat2$gene== "ORF6"|dat2$gene== "ORF7a"|
               dat2$gene== "ORF7b"|dat2$gene== "ORF8"|dat2$gene== "ORF9b")

dat3$gene <-factor(dat3$gene,levels= c('E','M','N','ORF1a','ORF1b','ORF3a','ORF6',
                                     'ORF7a','ORF7b','ORF8','ORF9b',ordered=TRUE))
p1 <-ggplot(dat3,aes(x=month,y=dnds_mean))+
  #geom_rect(aes(xmin="2021-01", xmax="2022-09", ymin=-Inf, ymax=Inf),fill='grey80',alpha = 0.01)+
  geom_line(aes(x=month,y=dnds_mean,group=gene,color=gene),size=0.5,stat="identity")+
  geom_ribbon(aes(x=month,ymin = dnds_lci, ymax = dnds_uci,group=gene,fill=gene),alpha=0.7)+
  geom_point(aes(x=month, y=dnds_mean,group=gene,color=gene),size=2,shape=19,alpha=0.9)+
  facet_wrap(~gene,scales = "free",ncol=4)+
  scale_y_continuous(limits=c(0,2),breaks = seq(0,2,0.5))+
  #labels = seq(0,100,25),expand = c(0, 0))+
  labs(x="", y= 'Ratio of nonsynonymous to synonymous divergence',color="")+
  #scale_x_discrete(labels=NULL)+
  scale_color_manual(name = '',values = colors1)+ 
  scale_fill_manual(name = 'Gene',values = colors1)+ 
  mytheme+
  theme(axis.title.y.left = element_text(size = 20,color="black",vjust=2,face = "bold", margin = margin(0,0.4,0,0,'cm')),
        axis.title.x = element_text(size = 20, color="black",face = "bold",margin = margin(0.4,0,0,0,'cm')))+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
    panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+
  #geom_vline(xintercept = "2021-01",lty='dashed',colour ="grey50",size=1)+
  scale_x_discrete(breaks=c("2020-04","2020-11","2021-06","2022-01",
                            "2022-08"),
                   labels=c("Apr\n2020","Nov\n2020","Jun\n2021","Jan\n2022",
                            "Aug\n2022"))+
  theme(legend.position = "none")+
  theme(panel.spacing = unit(1, "cm"))+
  theme(strip.text = element_text(face="bold",size=20))+
  annotate(geom = "text", x = "2021-01", y =1.05,angle=90, label = "Gamma", hjust = 0,size =4.8)+
  annotate(geom = "text", x = "2021-05", y =1.05,angle=90, label = "Delta", hjust =0,size = 4.8)+
  annotate(geom = "text", x = "2020-12", y =1.05,angle=90, label = "Alpha/Beta", hjust = 0,size = 4.8)+
  annotate(geom = "text", x = "2021-11", y =1.05,angle=90, label = "Omicron", hjust =0,size = 4.8)+
  #annotate(geom = "text", x = "2020-04", y =6.8,angle=0, label = "International travel restrictions", hjust ="left",size = 4.5,family = "ArialN")+
  #annotate(geom = "text", x = "2020-04", y =5.8,angle=0, label = "Public health and social measures", hjust ="left",size = 4.5,family = "ArialN")+
  annotate(geom = "text", x = "2021-06", y =1.9,angle=0, label = "Vaccine rollout", hjust =0,size = 5.5)+
  #geom_vline(aes(xintercept="2021-01"),lty='solid',colour ="black",size=0.5)+
  #geom_vline(aes(xintercept="2021-05"),lty='solid',colour ="black",size=0.5)+
  #geom_vline(aes(xintercept="2020-12"),lty='solid',colour ="black",size=0.5)+
  #geom_vline(aes(xintercept="2021-11"),lty='solid',colour ="black",size=0.5)+
  #geom_vline(aes(xintercept="2021-01"),lty='solid',colour ="black",size=0.5)+
  #geom_vline(aes(xintercept="2022-09"),lty='solid',colour ="black",size=0.5)+
  geom_hline(aes(yintercept=1),lty='dashed',colour ="black",size=1)
p1

ggsave(p1,filename =file.path(paste0(outpath,"fig_S2.pdf")),
       width = 20,height =12.5)
