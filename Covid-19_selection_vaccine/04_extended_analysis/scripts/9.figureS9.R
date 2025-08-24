# # install.packages("lubridate")
rm(list = ls())
root_path <- "C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine"
setwd(root_path)
dat<-read.csv("01_data/subsample_adaptation_vaccine_open.csv")
dat1 <- dat

source("load_package.R")
source("plot_theme.R")
source("load_data_clean.R")

outpath0 <-  "04_extended_analysis/plot_output/Figure_S9"
dir.create(outpath0)

dat1_S1= subset(dat1,dat1$gene== "S1")


var <- c("month","country","travel_route")
alldata <- dat1_S1
alldata$country <- as.character(alldata$country)
country_name <- unique(alldata$country)
print(country_name)
m =length(country_name)

ourdata <-alldata[,var] 
colnames(ourdata) <- c("month","country","value")

res_ci <- matrix(NA,m,4)
res_ci <- as.data.frame(res_ci)

for (j in 1:m)
{
  print(country_name[j]);
  result<- matrix(NA,1,3)
  result <- as.data.frame(result) 
  Q2 = subset(ourdata,trimws(ourdata$country)== country_name[j])
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
  colnames(result) <- c("travel_lci","travel_mean","travel_uci")
  result$country <- country_name[j]
  #result$Method <-method_type[q]
  
  res_ci[j,] <- result
}
colnames(res_ci) <- colnames(result)

res_ci1 <- res_ci
res_ci$type <- NA

res_ci[which(res_ci$travel_mean<=3000),"type"] <- "Monthly travel volume less than 3000"
res_ci[which(res_ci$travel_mean>3000),"type"] <- "Monthly travel volume more than 3000"



res_ci$country <- factor(res_ci$country,
                                 levels = res_ci$country[order(res_ci$travel_mean, decreasing = F)],ordered=TRUE)

unique(res_ci$country)

x3 <- which(res_ci$country=="Mexico")

res_ci <-res_ci[-x3,]

custom_colors <- c("#c25b88","#7570b3","#3cb371","#ee7600","#eec900","#1874cd","#7d26cd")




ps1 <- ggplot(res_ci,aes(x=country,y=travel_mean,group=1))+
  geom_bar(aes(x=country,y=travel_mean,fill=type),stat="identity",alpha=0.75)+
  scale_fill_manual("",values=custom_colors)+
  geom_errorbar(aes(x=country, y=travel_mean,ymin=travel_lci,ymax=travel_uci),
                width=0.2,stat="identity",position = position_dodge(0.5),lwd=0.6,alpha=1,color="black")+
  #scale_y_continuous(limits=c(0,102),breaks = seq(0,100,25))+
  #guides(fill = "none",color="none")+
  guides(fill=guide_legend(nrow=1,byrow = T))+
  xlab("") +
  ylab("International travel") +
  mytheme+
  #geom_hline(aes(yintercept=1000),lty='solid',colour ="black",size=1)+
  geom_hline(aes(yintercept=3000),lty='dashed',colour ="black",size=1)+
  theme(panel.grid=element_line(linetype=3,color='grey'))+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+
  theme(axis.text.x = element_text(size = 18, color="black",vjust = 0.5,hjust = 1, angle = 90))

#ggtitle("Natural immunity")+
ps1

ps2 <- ggplot(res_ci, aes(x=type, y=travel_mean,fill=type)) + 
  geom_boxplot(notch=F, outlier.shape = 19,outlier.size = 2,
               outlier.fill = "black",alpha=0.45, width = 0.3)+
  labs(x="",
       y= "")+
  scale_y_continuous(limits=c(0,12000),breaks = seq(0,12000,4000))+
  scale_fill_manual(name = '',values = custom_colors)+
  scale_color_manual(name = '',values = custom_colors)+
  geom_signif(comparisons = list(c("Monthly travel volume less than 3000", "Monthly travel volume more than 3000")),
              map_signif_level=F, textsize=6, y_position = 11000,test = "wilcox.test",
              color="black")+
  mytheme+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+
  #theme(axis.text.x = element_text(angle = 30, hjust = 0.75,vjust = 0.75))+
  theme(panel.grid=element_line(linetype=3,color='grey'))+
  theme(plot.margin=unit(c(0.5,1,0.5,1),'lines'))+
  guides(fill = "none",color="none")
ps2  

ggsave(ps1,filename = file.path(outpath0,"Figure_S9A1.pdf"),
       width = 17.75,height =7.5)

ggsave(ps2,filename = file.path(outpath0,"Figure_S9A2.pdf"),
       width = 4.2,height =3.7)

it_country <- res_ci[which(res_ci$travel_mean<=3000),]


#write.csv(it_country,paste0(outpath0,"/country_traver_low_2000.csv"),row.names =FALSE)






##Remove November 2021 to September 2022 as the Omicron variant emerged and circulated during this period.
my1 <- which(dat1$month=="2021-11"|dat1$month=="2021-12"|dat1$month=="2022-01"|dat1$month=="2022-02"|
               dat1$month=="2022-03"|dat1$month=="2022-04"|dat1$month=="2022-05"|dat1$month=="2022-06"|
               dat1$month=="2022-07"|dat1$month=="2022-08"|dat1$month=="2022-09") 
dat1 <-dat1[-my1,]


dat1$gene <- trimws(dat1$gene)
gene_type <- unique(dat1$gene)
s=length(gene_type)
print(gene_type)


outpath1 <-  paste0(outpath0,"/","dnds_vaccine_202003_202110")
dir.create(outpath1)


outpath2 <-  paste0(outpath0,"/","muts_nonsyn_vaccine_202003_202110")
dir.create(outpath2)


outpath3 <-  paste0(outpath0,"/","fitness_vaccine_202003_202110")
dir.create(outpath3)


##Remove Africa countries and zero-Covid-19 countries
dat1<- dat1[dat1$country %in% it_country$country, ]


#Model the association between adjusted vaccine coverage and the nonsynonymous to synonymous divergence ratios from March 2020 to October 2021

for(g in 1:s)
{
  print(gene_type[g])
  data= subset(dat1,dat1$gene== gene_type[g])
  #dat2= subset(dat1,dat1$gene== "S")
  
  source('load-data-cb-function.R')
  
  # set data for models of ratio of nonsynonymous to synonymous divergence
  Y1  <- data$dnds_mean_focal # response variable for ratio of nonsynonymous to synonymous divergence
  #Y2  <- data$muts_nonsyn_focal # response variable for number of nonsynonymous mutations
  #Y3  <- data$fitness_focal # response variable for mutational fitness
  P  <- data$NPI_NV # forpublic health and social measure (PHSM)
  V <- data$FW_BW_P # for adjusted vaccine coverage
  N <- log(data$new_cases_month) # Log transformation of the number of new cases per month 
  C <- data$country_index # for country interaction with month random effect
  C1 <- data$country # for country interaction with month random effect
  M <- data$month # The month index of the dataset 
  NI <- data$NI_IHME_P_omicron # for natural immunity
  
  df <- data.frame(Y1, P, V, N, C,C1,NI, M)
  df[which(df$N=="-Inf"),"N"] <- NA
  
  A <- round(max(df$V,na.rm = T),2)*10
  #A <- 75*10
  
  
  model1 <-  gam(Y1 ~ npi_cb+FBW_cb+travel_cb+N+s(C, bs = "re")+P:V,
                 method = "REML",data=df)
  summary(model1)
  AIC(model1)
  cp_FW <- crosspred(FBW_cb, model1, cen=0,by=0.1,at=0:A/10,bylag=0.1,cumul=TRUE) #slices
  cp_npi <- crosspred(npi_cb, model1,cen=0,at=0:1000/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
  cp_travel <- crosspred(travel_cb, model1, cen=0,at=0:150/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
  
  #assign(paste("model_vero_", gene_type[g], sep=""),model1)
  #assign(paste("cp_cases_", gene_type[g], sep=""),cp_cases)
  #assign(paste("cp_FW_vero_", gene_type[g], sep=""),cp_FW)
  
  
  # Extract the all estimates for the association between dN/dS ratio and adjusted vaccine coverage
  #dat_vero_low <- cp_FW$cumlow[,"lag3"]
  dat_vero_low <- as.data.frame(cp_FW$cumlow[,"lag3"])
  colnames(dat_vero_low) <- "dat_vero_low"
  dat_vero_fit <- as.data.frame(cp_FW$cumfit[,"lag3"])
  colnames(dat_vero_fit) <- "dat_vero_fit"
  dat_vero_high <- as.data.frame(cp_FW$cumhigh[,"lag3"])
  colnames(dat_vero_high) <- "dat_vero_high"
  result_vero_all <- cbind(dat_vero_low,dat_vero_fit,dat_vero_high)
  result_vero_all$coverage <- rownames(result_vero_all)
  result_vero_all <- as.data.frame(result_vero_all)
  result_vero_all$gene <- gene_type[g]
  write.csv(result_vero_all,file=paste0(outpath1,"/vero_dnds_slice_",gene_type[g],".csv"),row.names =F)
  
  
}



#Model the association between adjusted vaccine coverage and number of nonsynonymous mutations from March 2020 to October 2021

for(g in 1:s)
{
  
  print(gene_type[g])
  data= subset(dat1,dat1$gene== gene_type[g])
  #dat2= subset(dat1,dat1$gene== "S")
  
  source('load-data-cb-function.R')
  
  # set data for models of ratio of nonsynonymous to synonymous divergence
  #Y1  <- data$dnds_mean_focal # response variable for ratio of nonsynonymous to synonymous divergence
  Y1  <- data$muts_nonsyn_focal # response variable for number of nonsynonymous mutations
  #Y3  <- data$fitness_focal # response variable for mutational fitness
  P  <- data$NPI_NV # forpublic health and social measure (PHSM)
  V <- data$FW_BW_P # for adjusted vaccine coverage
  N <- log(data$new_cases_month) # Log transformation of the number of new cases per month 
  C <- data$country_index # for country interaction with month random effect
  C1 <- data$country # for country interaction with month random effect
  M <- data$month # The month index of the dataset 
  NI <- data$NI_IHME_P_omicron # for natural immunity
  
  df <- data.frame(Y1, P, V, N, C,C1,NI, M)
  df[which(df$N=="-Inf"),"N"] <- NA
  
  A <- round(max(df$V,na.rm = T),2)*10
  #A <- 75*10
  
  
  model1 <-  gam(Y1 ~ npi_cb+FBW_cb+travel_cb+N+s(C, bs = "re")+P:V,
                 method = "REML",data=df)
  summary(model1)
  AIC(model1)
  cp_FW <- crosspred(FBW_cb, model1, cen=0,by=0.1,at=0:A/10,bylag=0.1,cumul=TRUE) #slices
  cp_npi <- crosspred(npi_cb, model1,cen=0,at=0:1000/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
  cp_travel <- crosspred(travel_cb, model1, cen=0,at=0:150/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
  
  #assign(paste("model_vero_", gene_type[g], sep=""),model1)
  #assign(paste("cp_cases_", gene_type[g], sep=""),cp_cases)
  #assign(paste("cp_FW_vero_", gene_type[g], sep=""),cp_FW)
  
  
  # Extract the all estimates for the association between number of nonsynonymous mutations and adjusted vaccine coverage
  #dat_vero_low <- cp_FW$cumlow[,"lag3"]
  dat_vero_low <- as.data.frame(cp_FW$cumlow[,"lag3"])
  colnames(dat_vero_low) <- "dat_vero_low"
  dat_vero_fit <- as.data.frame(cp_FW$cumfit[,"lag3"])
  colnames(dat_vero_fit) <- "dat_vero_fit"
  dat_vero_high <- as.data.frame(cp_FW$cumhigh[,"lag3"])
  colnames(dat_vero_high) <- "dat_vero_high"
  result_vero_all <- cbind(dat_vero_low,dat_vero_fit,dat_vero_high)
  result_vero_all$coverage <- rownames(result_vero_all)
  result_vero_all <- as.data.frame(result_vero_all)
  result_vero_all$gene <- gene_type[g]
  write.csv(result_vero_all,file=paste0(outpath2,"/vero_muts_nonsyn_slice_",gene_type[g],".csv"),row.names =F)
  
  
}




#Model the association between adjusted vaccine coverage and mutational fitness from March 2020 to October 2021

for(g in 14)
{
  
  print(gene_type[g])
  data= subset(dat1,dat1$gene== gene_type[g])
  #dat2= subset(dat1,dat1$gene== "S")
  
  source('load-data-cb-function.R')
  
  # set data for models of mutational fitness
  #Y1  <- data$dnds_mean_focal # response variable for ratio of nonsynonymous to synonymous divergence
  #Y1  <- data$muts_nonsyn_focal # response variable for number of nonsynonymous mutations
  Y1  <- data$fitness_focal # response variable for mutational fitness
  P  <- data$NPI_NV # forpublic health and social measure (PHSM)
  V <- data$FW_BW_P # for adjusted vaccine coverage
  N <- log(data$new_cases_month) # Log transformation of the number of new cases per month 
  C <- data$country_index # for country interaction with month random effect
  C1 <- data$country # for country interaction with month random effect
  M <- data$month # The month index of the dataset 
  NI <- data$NI_IHME_P_omicron # for natural immunity
  
  df <- data.frame(Y1, P, V, N, C,C1,NI, M)
  df[which(df$N=="-Inf"),"N"] <- NA
  
  A <- round(max(df$V,na.rm = T),2)*10
  #A <- 75*10
  
  
  
  model1 <-  gam(Y1 ~ npi_cb+FBW_cb+travel_cb+N+s(C, bs = "re")+P:V,
                 method = "REML",data=df)
  summary(model1)
  AIC(model1)
  cp_FW <- crosspred(FBW_cb, model1, cen=0,by=0.1,at=0:A/10,bylag=0.1,cumul=TRUE) #slices
  #cp_NW <- crosspred(NIW_cb, model1, cen=0,by=0.1,at=0:B/10,bylag=0.1,cumul=TRUE) #slices
  cp_npi <- crosspred(npi_cb, model1,cen=0,at=0:1000/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
  cp_travel <- crosspred(travel_cb, model1, cen=0,at=0:150/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
  
  #assign(paste("model_vero_", gene_type[g], sep=""),model1)
  #assign(paste("cp_cases_", gene_type[g], sep=""),cp_cases)
  #assign(paste("cp_FW_vero_", gene_type[g], sep=""),cp_FW)
  
  
  # Extract the all estimates for the association between mutational fitness and adjusted vaccine coverage
  #dat_vero_low <- cp_FW$cumlow[,"lag3"]
  dat_vero_low <- as.data.frame(cp_FW$cumlow[,"lag3"])
  colnames(dat_vero_low) <- "dat_vero_low"
  dat_vero_fit <- as.data.frame(cp_FW$cumfit[,"lag3"])
  colnames(dat_vero_fit) <- "dat_vero_fit"
  dat_vero_high <- as.data.frame(cp_FW$cumhigh[,"lag3"])
  colnames(dat_vero_high) <- "dat_vero_high"
  result_vero_all <- cbind(dat_vero_low,dat_vero_fit,dat_vero_high)
  result_vero_all$coverage <- rownames(result_vero_all)
  result_vero_all <- as.data.frame(result_vero_all)
  result_vero_all$gene <- gene_type[g]
  write.csv(result_vero_all,file=paste0(outpath3,"/vero_fitness_slice_",gene_type[g],".csv"),row.names =F)
  
  
  
}


#Plot

DF1 <- read.csv(paste0(root_path,"/",outpath1,"/","vero_dnds_slice_S1.CSV"))
colnames(DF1) <- c("S_dat_vero_low","S_dat_vero_fit","S_dat_vero_high","vero_coverage" ,"S")


#custom_colors <- c("#ED0000FF", "#00468BFF")
custom_colors <- c( "#00468BFF","#5f9723")



p1 <-ggplot(DF1,aes(x=vero_coverage,y=S_dat_vero_fit))+
  #geom_rect(aes(xmin=57.2, xmax=75, ymin=-Inf, ymax=Inf),fill='grey90',alpha = 0.01)+
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",linewidth=1)+
  #geom_vline(aes(xintercept=25),lty='dashed',colour ="#92C5DE",linewidth=1)+
  #geom_vline(aes(xintercept=50),lty='dashed',colour ="#4393C3",linewidth=1)+
  #geom_vline(aes(xintercept=75),lty='dashed',colour ="#2166AC",linewidth=1)+
  geom_line(aes(x=vero_coverage,y=S_dat_vero_fit,group=1,color="Pre-Omicron period"),linewidth=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = S_dat_vero_low, ymax = S_dat_vero_high,group=1,fill="Pre-Omicron period"),alpha=0.25)+
  
  scale_x_continuous(limits=c(0,61),breaks = seq(0,60,20))+
  scale_y_continuous(limits=c(-19.8,6),breaks = seq(-18,6,6))+
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
  theme(legend.position = "none")

p1



DF2 <- read.csv(paste0(root_path,"/",outpath2,"/","vero_muts_nonsyn_slice_S1.CSV"))
colnames(DF2) <- c("S1_dat_vero_low","S1_dat_vero_fit","S1_dat_vero_high","vero_coverage" ,"S1")


p2 <-ggplot(DF2,aes(x=vero_coverage,y=S1_dat_vero_fit))+
  #geom_rect(aes(xmin=58.6, xmax=76, ymin=-Inf, ymax=Inf),fill='grey90',alpha = 0.01)+
  #geom_vline(aes(xintercept=25),lty='dashed',colour ="#92C5DE",linewidth=1)+
  #geom_vline(aes(xintercept=50),lty='dashed',colour ="#4393C3",linewidth=1)+
  #geom_vline(aes(xintercept=75),lty='dashed',colour ="#2166AC",linewidth=1)+
  
  geom_line(aes(x=vero_coverage,y=S1_dat_vero_fit,group=1,color="Pre-Omicron period"),size=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = S1_dat_vero_low, ymax = S1_dat_vero_high,group=1,fill="Pre-Omicron"),alpha=0.25)+
  
  scale_x_continuous(limits=c(0,61),breaks = seq(0,60,20))+
  scale_y_continuous(limits=c(-66,20),breaks = seq(-60,20,20))+
  labs(x="Adjusted vaccine coverage", y= "Effect on number \nof nonsynonymous mutations" )+  
  ggtitle("")+
  #scale_x_discrete(labels=NULL)+
  scale_color_manual(name = '',values = custom_colors)+ 
  scale_fill_manual(name = '',values = custom_colors)+ 
  mytheme+
  theme(panel.grid=element_line(linetype=2,color='grey'))+
  guides(color=guide_legend(nrow=1),fill=guide_legend(nrow=1))+
  theme(legend.position = "none")+
  theme(plot.margin=unit(c(0.75,1,0.75,1.5),'lines'))+
  #theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
  #      panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+ 
  theme(axis.title.y.left = element_text( margin = margin(0,0.1,0,0,'cm')),
        axis.title.x = element_text(margin = margin(0.1,0,0,0,'cm')))+
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",size=1)
p2



DF3 <- read.csv(paste0(root_path,"/",outpath3,"/","vero_fitness_slice_S1.CSV"))
colnames(DF3) <- c("S2_dat_vero_low","S2_dat_vero_fit","S2_dat_vero_high","vero_coverage" ,"S2")


p3 <-ggplot(DF3,aes(x=vero_coverage,y=S2_dat_vero_fit))+
  #geom_rect(aes(xmin=51, xmax=75, ymin=-Inf, ymax=Inf),fill='grey90',alpha = 0.01)+
  #geom_vline(aes(xintercept=25),lty='dashed',colour ="#92C5DE",linewidth=1)+
  #geom_vline(aes(xintercept=50),lty='dashed',colour ="#4393C3",linewidth=1)+
  #geom_vline(aes(xintercept=75),lty='dashed',colour ="#2166AC",linewidth=1)+
  
  geom_line(aes(x=vero_coverage,y=S2_dat_vero_fit,group=1,color="Pre-Omicron period"),size=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = S2_dat_vero_low, ymax = S2_dat_vero_high,group=1,fill="Pre-Omicron"),alpha=0.25)+
  scale_x_continuous(limits=c(0,61),breaks = seq(0,60,20))+
  scale_y_continuous(limits=c(-6.6,2),breaks = seq(-6,2,2))+
  #labels = seq(0,100,25),expand = c(0, 0))+
  labs(x="Adjusted vaccine coverage", y= "Effect on mutational fitness" )+  
  ggtitle("")+
  #scale_x_discrete(labels=NULL)+
  scale_color_manual(name = '',values = custom_colors)+ 
  scale_fill_manual(name = '',values = custom_colors)+ 
  mytheme+
  theme(panel.grid=element_line(linetype=2,color='grey'))+
  guides(color=guide_legend(nrow=1),fill=guide_legend(nrow=1))+
  theme(legend.position = "none")+
  #theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
  #      panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+ 
  #theme(legend.position = "none")+
  theme(plot.margin=unit(c(0.75,1,0.75,1.5),'lines'))+
  theme(axis.title.y.left = element_text( margin = margin(0,0.1,0,0,'cm')),
        axis.title.x = element_text(margin = margin(0.1,0,0,0,'cm')))+
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",size=1)
p3



p<-ggarrange(p1, p2, p3,ncol =3, nrow = 1,widths = c(1,1,1),heights = c(1,1,1),
             labels = c("B","C","D"), align = "hv",
             legend = "none", common.legend = T,
             font.label = list(size = 30, face = "bold"))
p

ggsave(p,filename = file.path(outpath0,"figure_S9.pdf"),
       width = 18,height =6)


#p<-ggarrange(ps2,p1, p2, p3,ncol =2, nrow = 2,widths = c(1,1,1,1),heights = c(1,1,1,1),
#             labels = c("A","B","C","D"), align = "hv",
#             legend = "none", common.legend = T,
#             font.label = list(size = 30, face = "bold"))
#p

#ggsave(p,filename = file.path(outpath0,"figure_S9.pdf"),
#       width = 18,height =18)




