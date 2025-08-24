rm(list = ls())
root_path <- "C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine"
setwd(root_path)
dat1<-read.csv("01_data/subsample_adaptation_vaccine_low.csv")
source("load_package.R")
source("plot_theme.R")
source("load_data_clean.R")

dat2 <- dat1
##Remove November 2021 to September 2022 as the Omicron variant emerged and circulated during this period.
my1 <- which(dat1$month=="2021-11"|dat1$month=="2021-12"|dat1$month=="2022-01"|dat1$month=="2022-02"|
               dat1$month=="2022-03"|dat1$month=="2022-04"|dat1$month=="2022-05"|dat1$month=="2022-06"|
               dat1$month=="2022-07"|dat1$month=="2022-08"|dat1$month=="2022-09") 
dat1 <-dat1[-my1,]


dat1$gene <- trimws(dat1$gene)
gene_type <- unique(dat1$gene)
s=length(gene_type)
print(gene_type)

country1 <- read.csv("02_sup-data/country_mRNA_low_2110.csv")
country1 <- country1[,c("country","Region")]
country2 <- read.csv("02_sup-data/country_inactivated_vaccine_low_2110.csv")
country2 <- country2[,c("country","Region")]
country3 <- read.csv("02_sup-data/country_adenovirus_vector_low_2110.csv")
country3 <- country3[,c("country","Region")]



dat1_vero_mr <- merge(dat1,country1,by=c("country"),all.y = T)
dat1_vero_in <- merge(dat1,country2,by=c("country"),all.y = T)
dat1_vero_ad <- merge(dat1,country3,by=c("country"),all.y = T)

data1 <- dat1_vero_mr 
data2 <- dat1_vero_in 
data3 <- dat1_vero_ad 

data1$type <- "mRNA_vaccine"
data2$type <- "Inactivated_vaccine"
data3$type <- "Viral_vector_vaccine"
alldata <- rbind(data1,data2,data3) 

#Please note that we did not estimate associations between adjusted vaccine coverage for viral vector vaccines and selection pressure. 
#Because a substantial proportion of African countries (7 out of 12) administered viral vector vaccines.
#The vaccine coverage in these African countries was low.
outpath0 <-  "04_extended_analysis/plot_output/Figure_S16"
dir.create(outpath0)
outpath1 <-  paste0(outpath0,"/","dnds_vero_type_202110")
dir.create(outpath1)


outpath2 <-  paste0(outpath0,"/","muts_nonsyn_vero_type_202110")
dir.create(outpath2)


outpath3 <-  paste0(outpath0,"/","fitness_vero_type_202110")
dir.create(outpath3)

##Model the association between the dN/dS ratios and adjusted vaccine coverage with 
##different countries stratified by vaccine types from March 2020 to October 2021


vero_type <- unique(alldata$type)
n=length(vero_type)
print(vero_type)

for (g in 1:s)
  {
  print(gene_type[g])
  alldata1= subset(alldata,alldata$gene== gene_type[g])
 
for(i in 1:2)
{
print(vero_type[i])
dat2= subset(alldata1,alldata1$type== vero_type[i])
##Assign a unique index to each country.
dat2 <- dat2 %>%
  group_by(country) %>%
  mutate(country_index = cur_group_id())
dat2$country_index <- as.factor(dat2$country_index)
data <- dat2

source('load-data-cb-function.R')
  
# set data for models of ratio of nonsynonymous to synonymous divergence
Y1  <- data$dnds_mean_focal # response variable for ratio of nonsynonymous to synonymous divergence
#Y1  <- data$muts_nonsyn_focal # response variable for number of nonsynonymous mutations
#Y3  <- data$fitness_focal # response variable for mutational fitness
P  <- data$NPI_NV # forpublic health and social measure (PHSM)
V <- data$FW_BW_P # for adjusted vaccine coverage
N <- log(data$new_cases_month) # Log transformation of the number of new cases per month 
C <- data$country_index # for country interaction with month random effect
M <- data$month # The month index of the dataset 
NI <- data$NI_IHME_P_omicron # for natural immunity

df <- data.frame(Y1, P, V, N, C, M,NI)

model1 <-  gam(Y1 ~ npi_cb+FBW_cb+N+travel_cb+s(C, bs = "re")+P:V,
              method = "REML",data=df)

summary(model1)
AIC(model1)

A <- (round(max(df$V,na.rm=T),2))*10

cp_FW <- crosspred(FBW_cb, model1, cen=0,by=0.1,at=0:A/10,bylag=0.1,cumul=TRUE) #slices



dat_vero_low <- as.data.frame(cp_FW$cumlow[,"lag3"])
colnames(dat_vero_low) <- "dat_vero_low"
dat_vero_fit <- as.data.frame(cp_FW$cumfit[,"lag3"])
colnames(dat_vero_fit) <- "dat_vero_fit"
dat_vero_high <- as.data.frame(cp_FW$cumhigh[,"lag3"])
colnames(dat_vero_high) <- "dat_vero_high"
result <- cbind(dat_vero_low,dat_vero_fit,dat_vero_high)
result$coverage <- rownames(result)
result <- as.data.frame(result)
result$gene <- gene_type[g]
result$Type <- vero_type[i]

write.csv(result,file=paste0(outpath1,"/fit_vero_type_",gene_type[g],"_",vero_type[i],".csv"),row.names =F)


}

 
  
}





##Model the association between the number of nonsynonymous mutations and adjusted vaccine coverage with 
##different countries stratified by vaccine types from March 2020 to October 2021


for (g in 1:s)
{
  print(gene_type[g])
  alldata1= subset(alldata,alldata$gene== gene_type[g])
  
  for(i in 1:2)
  {
    print(vero_type[i])
    dat2= subset(alldata1,alldata1$type== vero_type[i])
    ##Assign a unique index to each country.
    dat2 <- dat2 %>%
      group_by(country) %>%
      mutate(country_index = cur_group_id())
    dat2$country_index <- as.factor(dat2$country_index)
    data <- dat2
    
    source('load-data-cb-function.R')
    
    # set data for models of ratio of nonsynonymous to synonymous divergence
    #Y1  <- data$dnds_mean_focal # response variable for ratio of nonsynonymous to synonymous divergence
    Y1  <- data$muts_nonsyn_focal # response variable for number of nonsynonymous mutations
    #Y3  <- data$fitness_focal # response variable for mutational fitness
    P  <- data$NPI_NV # forpublic health and social measure (PHSM)
    V <- data$FW_BW_P # for adjusted vaccine coverage
    N <- log(data$new_cases_month) # Log transformation of the number of new cases per month 
    C <- data$country_index # for country interaction with month random effect
    M <- data$month # The month index of the dataset 
    NI <- data$NI_IHME_P_omicron # for natural immunity
    
    df <- data.frame(Y1, P, V, N, C, M,NI)
    
    model1 <-  gam(Y1 ~ npi_cb+FBW_cb+N+travel_cb+s(C, bs = "re")+P:V,
                   method = "REML",data=df)
    
    summary(model1)
    AIC(model1)
    
    A <- (round(max(df$V,na.rm=T),2))*10

    cp_FW <- crosspred(FBW_cb, model1, cen=0,by=0.1,at=0:A/10,bylag=0.1,cumul=TRUE) #slices

    
    
    dat_vero_low <- as.data.frame(cp_FW$cumlow[,"lag3"])
    colnames(dat_vero_low) <- "dat_vero_low"
    dat_vero_fit <- as.data.frame(cp_FW$cumfit[,"lag3"])
    colnames(dat_vero_fit) <- "dat_vero_fit"
    dat_vero_high <- as.data.frame(cp_FW$cumhigh[,"lag3"])
    colnames(dat_vero_high) <- "dat_vero_high"
    result <- cbind(dat_vero_low,dat_vero_fit,dat_vero_high)
    result$coverage <- rownames(result)
    result <- as.data.frame(result)
    result$gene <- gene_type[g]
    result$Type <- vero_type[i]
    
    write.csv(result,file=paste0(outpath2,"/fit_vero_type_",gene_type[g],"_",vero_type[i],".csv"),row.names =F)
    
   
        
  }
  
  
  
}



##Model the association between mutational fitness and adjusted vaccine coverage with 
##different countries stratified by vaccine types from March 2020 to October 2021

for (g in 14)
{
  print(gene_type[g])
  alldata1= subset(alldata,alldata$gene== gene_type[g])
  
  for(i in 1:2)
  {
    print(vero_type[i])
    dat2= subset(alldata1,alldata1$type== vero_type[i])
    ##Assign a unique index to each country.
    dat2 <- dat2 %>%
      group_by(country) %>%
      mutate(country_index = cur_group_id())
    dat2$country_index <- as.factor(dat2$country_index)
    data <- dat2
    
    source('load-data-cb-function.R')
    
    # set data for models of ratio of nonsynonymous to synonymous divergence
    #Y1  <- data$dnds_mean_focal # response variable for ratio of nonsynonymous to synonymous divergence
    #Y1  <- data$muts_nonsyn_focal # response variable for number of nonsynonymous mutations
    Y1  <- data$fitness_focal # response variable for mutational fitness
    P  <- data$NPI_NV # forpublic health and social measure (PHSM)
    V <- data$FW_BW_P # for adjusted vaccine coverage
    N <- log(data$new_cases_month) # Log transformation of the number of new cases per month 
    C <- data$country_index # for country interaction with month random effect
    M <- data$month # The month index of the dataset 
    NI <- data$NI_IHME_P_omicron # for natural immunity
    
    df <- data.frame(Y1, P, V, N, C, M,NI)
    
    model1 <-  gam(Y1 ~ npi_cb+FBW_cb+N+travel_cb+s(C, bs = "re")+P:V,
                   method = "REML",data=df)
    
    summary(model1)
    AIC(model1)
    
    A <- (round(max(df$V,na.rm=T),2))*10

    cp_FW <- crosspred(FBW_cb, model1, cen=0,by=0.1,at=0:A/10,bylag=0.1,cumul=TRUE) #slices

   
    
    dat_vero_low <- as.data.frame(cp_FW$cumlow[,"lag3"])
    colnames(dat_vero_low) <- "dat_vero_low"
    dat_vero_fit <- as.data.frame(cp_FW$cumfit[,"lag3"])
    colnames(dat_vero_fit) <- "dat_vero_fit"
    dat_vero_high <- as.data.frame(cp_FW$cumhigh[,"lag3"])
    colnames(dat_vero_high) <- "dat_vero_high"
    result <- cbind(dat_vero_low,dat_vero_fit,dat_vero_high)
    result$coverage <- rownames(result)
    result <- as.data.frame(result)
    result$gene <- gene_type[g]
    result$Type <- vero_type[i]
    
    write.csv(result,file=paste0(outpath3,"/fit_vero_type_",gene_type[g],"_",vero_type[i],".csv"),row.names =F)
  
        
  }
  
  
  
}




#Plot

S_In <- read.csv(paste0(root_path,"/",outpath1,"/","fit_vero_type_S1_Inactivated_vaccine.CSV"))
colnames(S_In) <- c("In_low","In_fit","In_high","vero_coverage" ,"S","In_type")

S_mRNA <- read.csv(paste0(root_path,"/",outpath1,"/","fit_vero_type_S1_mRNA_vaccine.CSV"))
colnames(S_mRNA) <- c("mRNA_low","mRNA_fit","mRNA_high","vero_coverage" ,"mRNA_type")

#custom_colors <- c("#ED0000FF", "#00468BFF")
custom_colors <- c("#f1a340", "#998ec3")

DF1 <- merge(S_In,S_mRNA,by=c("vero_coverage"),all.y = T)



p1 <-ggplot(DF1,aes(x=vero_coverage,y=In_fit))+
  #geom_rect(aes(xmin=54.6, xmax=76, ymin=-Inf, ymax=Inf),fill='grey90',alpha = 0.01)+
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",linewidth=1)+
  #geom_vline(aes(xintercept=25),lty='dashed',colour ="#92C5DE",linewidth=1)+
  #geom_vline(aes(xintercept=50),lty='dashed',colour ="#4393C3",linewidth=1)+
  #geom_vline(aes(xintercept=75),lty='dashed',colour ="#2166AC",linewidth=1)+
  geom_line(aes(x=vero_coverage,y=In_fit,group=1,color="Inactivated vaccine"),linewidth=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = In_low, ymax = In_high,group=1,fill="Inactivated vaccine"),alpha=0.25)+
  geom_line(aes(x=vero_coverage,y=mRNA_fit,group=1,color="mRNA vaccine"),linewidth=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = mRNA_low, ymax = mRNA_high,group=1,fill="mRNA vaccine"),alpha=0.25)+
  scale_x_continuous(limits=c(0,80),breaks = seq(0,75,25))+
  scale_y_continuous(limits=c(-4.2,4.5),breaks = seq(-4,4,2))+
  #labels = seq(0,100,25),expand = c(0, 0))+
  labs(x="Adjusted vaccine coverage", y= "Effect on nonsynonymous to \nsynonymous divergence ratio" )+  
  ggtitle("")+
  #scale_x_discrete(labels=NULL)+
  scale_color_manual(name = '',values = custom_colors)+ 
  scale_fill_manual(name = '',values = custom_colors)+ 
  mytheme+
  theme(plot.margin=unit(c(0.75,1,0.75,1.5),'lines'))+
  #theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
  #      panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+ 
  theme(panel.grid=element_line(linetype=2,color='grey'))+
  guides(color=guide_legend(nrow=1),fill=guide_legend(nrow=1))+
  theme(axis.title.y.left = element_text( margin = margin(0,0.1,0,0,'cm')),
        axis.title.x = element_text(margin = margin(0.1,0,0,0,'cm')))+
  theme(legend.position = "bottom")

p1

S_In <- read.csv(paste0(root_path,"/",outpath2,"/","fit_vero_type_S1_Inactivated_vaccine.CSV"))
colnames(S_In) <- c("In_low","In_fit","In_high","vero_coverage" ,"S","In_type")

S_mRNA <- read.csv(paste0(root_path,"/",outpath2,"/","fit_vero_type_S1_mRNA_vaccine.CSV"))
colnames(S_mRNA) <- c("mRNA_low","mRNA_fit","mRNA_high","vero_coverage" ,"mRNA_type")


DF2 <- merge(S_In,S_mRNA,by=c("vero_coverage"),all.y = T)


p2 <-ggplot(DF2,aes(x=vero_coverage,y=In_fit))+
  #geom_rect(aes(xmin=54.6, xmax=76, ymin=-Inf, ymax=Inf),fill='grey90',alpha = 0.01)+
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",linewidth=1)+
  #geom_vline(aes(xintercept=25),lty='dashed',colour ="#92C5DE",linewidth=1)+
  #geom_vline(aes(xintercept=50),lty='dashed',colour ="#4393C3",linewidth=1)+
  #geom_vline(aes(xintercept=75),lty='dashed',colour ="#2166AC",linewidth=1)+
  geom_line(aes(x=vero_coverage,y=In_fit,group=1,color="Inactivated vaccine"),linewidth=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = In_low, ymax = In_high,group=1,fill="Inactivated vaccine"),alpha=0.25)+
  geom_line(aes(x=vero_coverage,y=mRNA_fit,group=1,color="mRNA vaccine"),linewidth=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = mRNA_low, ymax = mRNA_high,group=1,fill="mRNA vaccine"),alpha=0.25)+
  scale_x_continuous(limits=c(0,80),breaks = seq(0,75,25))+
  scale_y_continuous(limits=c(-17,34),breaks = seq(-17,34,17))+
  #labels = seq(0,100,25),expand = c(0, 0))+
  labs(x="Adjusted vaccine coverage", y= "Effect on nonsynonymous to \nsynonymous divergence ratio" )+  
  ggtitle("")+
  #scale_x_discrete(labels=NULL)+
  scale_color_manual(name = '',values = custom_colors)+ 
  scale_fill_manual(name = '',values = custom_colors)+ 
  mytheme+
  theme(plot.margin=unit(c(0.75,1,0.75,1.5),'lines'))+
  #theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
  #      panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+ 
  theme(panel.grid=element_line(linetype=2,color='grey'))+
  guides(color=guide_legend(nrow=1),fill=guide_legend(nrow=1))+
  theme(axis.title.y.left = element_text( margin = margin(0,0.1,0,0,'cm')),
        axis.title.x = element_text(margin = margin(0.1,0,0,0,'cm')))+
  theme(legend.position = "bottom")

p2

S_In <- read.csv(paste0(root_path,"/",outpath3,"/","fit_vero_type_S1_Inactivated_vaccine.CSV"))
colnames(S_In) <- c("In_low","In_fit","In_high","vero_coverage" ,"S","In_type")

S_mRNA <- read.csv(paste0(root_path,"/",outpath3,"/","fit_vero_type_S1_mRNA_vaccine.CSV"))
colnames(S_mRNA) <- c("mRNA_low","mRNA_fit","mRNA_high","vero_coverage" ,"mRNA_type")


DF3 <- merge(S_In,S_mRNA,by=c("vero_coverage"),all.y = T)


p3 <-ggplot(DF3,aes(x=vero_coverage,y=In_fit))+
  #geom_rect(aes(xmin=54.6, xmax=76, ymin=-Inf, ymax=Inf),fill='grey90',alpha = 0.01)+
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",linewidth=1)+
  #geom_vline(aes(xintercept=25),lty='dashed',colour ="#92C5DE",linewidth=1)+
  #geom_vline(aes(xintercept=50),lty='dashed',colour ="#4393C3",linewidth=1)+
  #geom_vline(aes(xintercept=75),lty='dashed',colour ="#2166AC",linewidth=1)+
  geom_line(aes(x=vero_coverage,y=In_fit,group=1,color="Inactivated vaccine"),linewidth=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = In_low, ymax = In_high,group=1,fill="Inactivated vaccine"),alpha=0.25)+
  geom_line(aes(x=vero_coverage,y=mRNA_fit,group=1,color="mRNA vaccine"),linewidth=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = mRNA_low, ymax = mRNA_high,group=1,fill="mRNA vaccine"),alpha=0.25)+
  scale_x_continuous(limits=c(0,80),breaks = seq(0,75,25))+
  scale_y_continuous(limits=c(-1,3),breaks = seq(-1,3,1))+
  #labels = seq(0,100,25),expand = c(0, 0))+
  labs(x="Adjusted vaccine coverage", y= "Effect on mutational fitness" )+  
  ggtitle("")+
  #scale_x_discrete(labels=NULL)+
  scale_color_manual(name = '',values = custom_colors)+ 
  scale_fill_manual(name = '',values = custom_colors)+ 
  mytheme+
  theme(plot.margin=unit(c(0.75,1,0.75,1.5),'lines'))+
  #theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
  #      panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+ 
  theme(panel.grid=element_line(linetype=2,color='grey'))+
  guides(color=guide_legend(nrow=1),fill=guide_legend(nrow=1))+
  theme(axis.title.y.left = element_text( margin = margin(0,0.1,0,0,'cm')),
        axis.title.x = element_text(margin = margin(0.1,0,0,0,'cm')))+
  theme(legend.position = "bottom")

p3



p<-ggarrange(p1, p2, p3,ncol =3, nrow = 1,widths = c(1,1,1),heights = c(1,1,1),
             labels = c("A","B","C"), align = "hv",
             legend = "bottom", common.legend = T,
             font.label = list(size = 30, face = "bold"))
p

ggsave(p,filename = file.path(outpath0,"figure_S16.pdf"),
       width = 18.2,height =6.2)


