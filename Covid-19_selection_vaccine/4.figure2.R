rm(list = ls())
root_path <- "C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine"
setwd(root_path)
dat1<-read.csv("01_data/subsample_dnds_vaccine_diversity_low_death.csv")
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

#outpath <-  "02_model-output/all_vaccine_202003_202110"
outpath0 <- "03_output/figure_2_3"
dir.create(outpath0)
outpath <-  paste0(outpath0,"/","all_vaccine_202003_202110")
dir.create(outpath)
outpath1 <-  paste0(outpath,"/","vero_result",sep="")
dir.create(outpath1)
outpath2 <-  paste0(outpath0,"/","natural_immunity_202003_202110")
dir.create(outpath2)
outpath3 <-   paste0(outpath2,"/","NI_result",sep="")
dir.create(outpath3)

outpath4 <-  paste0(outpath0,"/","all_vaccine_202201_202209")
dir.create(outpath4)

outpath5 <-  paste0(outpath4,"/","Omicron_result",sep="")
dir.create(outpath5)


#Model the association between adjusted vaccine coverage and the nonsynonymous to synonymous divergence ratios

for(g in 1:s)
{
  
print(gene_type[g])
data= subset(dat1,dat1$gene== gene_type[g])
#dat2= subset(dat1,dat1$gene== "S")

source('load-data-cb-function.R')

# set data for models of growth rate of cases and lineage diversity
Y1  <- data$dn_ds_mean # response variable (growth rate of cases)
P  <- data$NPI_NV # forpublic health and social measure (PHSM)
V <- data$FW_BW_P # for adjusted vaccine coverage
N <- log(data$new_cases_month) # Log transformation of the number of new cases per month 
C <- data$country_index # for country interaction with month random effect
C1 <- data$country # for country interaction with month random effect
M <- data$month # The month index of the dataset 
df <- data.frame(Y1, P, V, N, C,C1, M)
df[which(df$N=="-Inf"),"N"] <- NA

A <- round(max(df$V,na.rm = T),2)*10

model1 <-  gam(Y1 ~ npi_cb+FBW_cb+travel_cb+N+s(C, bs = "re")+P:V,
              method = "REML",data=df)
summary(model1)
AIC(model1)
cp_FW <- crosspred(FBW_cb, model1, cen=0,by=0.1,at=0:A/10,bylag=0.1,cumul=TRUE) #slices

assign(paste("model_vero_", gene_type[g], sep=""),model1)
#assign(paste("cp_cases_", gene_type[g], sep=""),cp_cases)
assign(paste("cp_FW_vero_", gene_type[g], sep=""),cp_FW)


vero_low <- cp_FW$cumlow["25","lag3"]
vero_fit <- cp_FW$cumfit["25","lag3"]
vero_high <- cp_FW$cumhigh["25","lag3"]
result <- cbind(vero_low,vero_fit,vero_high)
result <- as.data.frame(result)
#result$Lag <- rownames(result)
result$coverage1 <- "25"
result$coverage <- "Low coverage (25)"
result$gene <- gene_type[g]
result_vero_25 <- result



vero_low <- cp_FW$cumlow["50","lag3"]
vero_fit <- cp_FW$cumfit["50","lag3"]
vero_high <- cp_FW$cumhigh["50","lag3"]
result <- cbind(vero_low,vero_fit,vero_high)
result <- as.data.frame(result)
#result$Lag <- rownames(result)
result$coverage1 <- "50"
result$coverage <- "Medium coverage (50)"
result$gene <- gene_type[g]
result_vero_50 <- result


vero_low <- cp_FW$cumlow["75","lag3"]
vero_fit <- cp_FW$cumfit["75","lag3"]
vero_high <- cp_FW$cumhigh["75","lag3"]
result <- cbind(vero_low,vero_fit,vero_high)
result <- as.data.frame(result)
#result$Lag <- rownames(result)
result$coverage1 <- "75"
result$coverage <- "High coverage (75)"
result$gene <- gene_type[g]
result_vero_75 <- result

result_vero <- rbind(result_vero_25,result_vero_50,result_vero_75)
write.csv(result_vero,file=paste0(outpath1,"/vero_pressure_selection_error_bar_",gene_type[g],".csv"),row.names =F)


#dat_vero_low <- cp_FW$cumlow[,"lag3"]
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
write.csv(result,file=paste0(outpath,"/vero_pressure_selection_slice_",gene_type[g],".csv"),row.names =F)




}




#Model the association between natural immunity coverage and the nonsynonymous to synonymous divergence ratios
full_vero_2110 <- dat1[which(dat1$gene=="S1"&dat1$month=="2021-10"),]

#Filter out countries with a full-dose COVID-19 vaccine coverage rate below 25% as of September 2022.
country_low <- full_vero_2110[which(full_vero_2110$F_Vero<=25),"country"]
country_low$country <- as.character(country_low$country)
dat1_ni<- dat1[dat1$country %in% country_low$country, ]

#country_low_vero <- full_vero_2110[which(full_vero_2110$F_Vero<=25),c("country","F_Vero")]
#write.csv(country_low_vero,"02_sup-data/country_less_25_low_ni.csv",row.names = F)


for(g in 1:s)
{
print(gene_type[g])
data= subset(dat1_ni,dat1_ni$gene== gene_type[g])
#dat2= subset(dat1,dat1$gene== "S")
source('load-data-cb-function.R')
  
# set data for models of growth rate of cases and lineage diversity
  Y1  <- data$dn_ds_mean # response variable (growth rate of cases)
  P  <- data$NPI_NV # forpublic health and social measure (PHSM)
  C <- data$country_index # for country interaction with month random effect
  M <- data$month # The month index of the dataset 
  NI <- data$NI_IHME_P_omicron # for natural immunity
  
  df <- data.frame(Y1, P, C, M, NI)

  
B <- (round(max(df$NI,na.rm=T),2))*10
#B <- 750
  
model2 <-  gam(Y1 ~ npi_cb+NIW2_cb+travel_cb+s(C, bs = "re"),
               method = "REML",data=df)
summary(model2)
AIC(model2)
cp_NW <- crosspred(NIW2_cb, model2, cen=0,by=0.1,at=0:B/10,bylag=0.1,cumul=TRUE) #slices


assign(paste("model_NI_",gene_type[g], sep=""),model1)
assign(paste("cp_NW_", gene_type[g], sep=""),cp_NW)


NI_low <- cp_NW$cumlow["25","lag3"]
NI_fit <- cp_NW$cumfit["25","lag3"]
NI_high <- cp_NW$cumhigh["25","lag3"]
result <- cbind(NI_low,NI_fit,NI_high)
result <- as.data.frame(result)
result$Lag <- rownames(result)
result$coverage1 <- "25"
result$coverage <- "Low coverage (25)"
result$gene <- gene_type[g]
result_NI_25 <- result

NI_low <- cp_NW$cumlow["50","lag3"]
NI_fit <- cp_NW$cumfit["50","lag3"]
NI_high <- cp_NW$cumhigh["50","lag3"]
result <- cbind(NI_low,NI_fit,NI_high)
result <- as.data.frame(result)
result$Lag <- rownames(result)
result$coverage1 <- "50"
result$coverage <- "Medium coverage (50)"
result$gene <- gene_type[g]
result_NI_50 <- result

result <- matrix(NA, nrow = 1,ncol = 3)
result <- as.data.frame(result)
colnames(result) <- c("NI_low","NI_fit","NI_high")
result <- as.data.frame(result)
result$Lag <- rownames(result)
result$coverage1 <- "75"
result$coverage <- "High coverage (75)"
result$gene <- gene_type[g]
result_NI_75 <- result


result_vero <- rbind(result_NI_25,result_NI_50,result_NI_75)
result_vero$Lag <- "lag3"

write.csv(result_vero,file=paste0(outpath3,"/natural_immunity_pressure_selection_error_bar_",gene_type[g],".csv"),row.names =F)


dat_NI_low <- as.data.frame(cp_NW$cumlow[,"lag3"])
colnames(dat_NI_low) <- "dat_NI_low"
dat_NI_fit <- as.data.frame(cp_NW$cumfit[,"lag3"])
colnames(dat_NI_fit) <- "dat_NI_fit"
dat_NI_high <- as.data.frame(cp_NW$cumhigh[,"lag3"])
colnames(dat_NI_high) <- "dat_NI_high"
result <- cbind(dat_NI_low,dat_NI_fit,dat_NI_high)
result$coverage <- rownames(result)
result <- as.data.frame(result)
result$gene <- gene_type[g]

write.csv(result,file=paste0(outpath2,"/NI_pressure_selection_slice_",gene_type[g],".csv"),row.names =F)




}


#Model the association between adjusted vaccine coverage and the nonsynonymous to synonymous divergence ratios for Omicron
my2 <- which(dat2$month=="2022-01"|dat2$month=="2022-02"|dat2$month=="2022-03"|dat2$month=="2022-04"|
              dat2$month=="2022-05"|dat2$month=="2022-06"|dat2$month=="2022-07"|dat2$month=="2022-08"|
              dat2$month=="2022-09") 
dat_Omicron <-dat2[my2,]


for(g in 1:s)
{
  
  print(gene_type[g])
  data= subset(dat_Omicron,dat_Omicron$gene== gene_type[g])
  #dat2= subset(dat1,dat1$gene== "S")

  source('load-data-cb-function.R')
  
  # set data for models of growth rate of cases and lineage diversity
  Y1  <- data$dn_ds_mean # response variable (growth rate of cases)
  
  P  <- data$NPI_NV # forpublic health and social measure (PHSM)
  V <- data$FW_BW_P # for adjusted vaccine coverage
  N <- log(data$new_cases_month) # Log transformation of the number of new cases per month 
  C <- data$country # for country interaction with month random effect
  M <- data$month # The month index of the dataset 
  
  df <- data.frame(Y1, P, V, N, C, M)
  df[which(df$N=="-Inf"),"N"] <- NA
  
  
  A <- round(max(df$V,na.rm = T),2)*10
  
  model1 <-  gam(Y1 ~ npi_cb+FBW_cb+N+travel_cb+s(C, bs = "re")+P:V,method = "REML",data=df)
  summary(model1)
  AIC(model1)
  cp_npi <- crosspred(npi_cb, model1, cen=0,by=0.1,at=0:1000/10,bylag=0.1,cumul=TRUE) #slices
  cp_travel <- crosspred(travel_cb, model1, cen=0,by=0.1,at=0:150/10,bylag=0.1,cumul=TRUE) #slices
  #cp_cases <- crosspred(cases_cb, model1,cen=0,by=0.1,bylag=0.1,at=0:180/10,cumul=TRUE) #slices
  cp_FW <- crosspred(FBW_cb, model1, cen=0,by=0.1,at=0:A/10,bylag=0.1,cumul=TRUE) #slices
  
  
  
  assign(paste("model_", gene_type[g], sep=""),model1)
  assign(paste("cp_npi_", gene_type[g], sep=""),cp_npi)
  assign(paste("cp_travel_", gene_type[g], sep=""),cp_travel)
  #assign(paste("cp_cases_", gene_type[g], sep=""),cp_cases)
  assign(paste("cp_FW_", gene_type[g], sep=""),cp_FW)
  
  #########fig.3A################
  
  
  
  vero_low <- cp_FW$cumlow["25","lag3"]
  vero_fit <- cp_FW$cumfit["25","lag3"]
  vero_high <- cp_FW$cumhigh["25","lag3"]
  result <- cbind(vero_low,vero_fit,vero_high)
  result <- as.data.frame(result)
  result$Lag <- rownames(result)
  result$coverage1 <- "25"
  result$coverage <- "Low coverage (25)"
  result$gene <- gene_type[g]
  result_vero_25 <- result
  
  
  
  vero_low <- cp_FW$cumlow["50","lag3"]
  vero_fit <- cp_FW$cumfit["50","lag3"]
  vero_high <- cp_FW$cumhigh["50","lag3"]
  result <- cbind(vero_low,vero_fit,vero_high)
  result <- as.data.frame(result)
  result$Lag <- rownames(result)
  result$coverage1 <- "50"
  result$coverage <- "Medium coverage (50)"
  result$gene <- gene_type[g]
  result_vero_50 <- result
  
  
  vero_low <- cp_FW$cumlow["75","lag3"]
  vero_fit <- cp_FW$cumfit["75","lag3"]
  vero_high <- cp_FW$cumhigh["75","lag3"]
  result <- cbind(vero_low,vero_fit,vero_high)
  result <- as.data.frame(result)
  result$Lag <- rownames(result)
  result$coverage1 <- "75"
  result$coverage <- "High coverage (75)"
  result$gene <- gene_type[g]
  result_vero_75 <- result
  
  result_vero <- rbind(result_vero_25,result_vero_50,result_vero_75)
  
  
  write.csv(result_vero,file=paste0(outpath5,"/","vero_pressure_selection_error_bar_",gene_type[g],".csv",sep=""),row.names =F)
  
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
  write.csv(result,file=paste0(outpath4,"/","vero_pressure_selection_slice_2022_",gene_type[g],".csv",sep=""),row.names =F)
  
  
  

}

#Plot

S <- read.csv(paste0(root_path,"/",outpath,"/","vero_pressure_selection_slice_S.CSV"))
colnames(S) <- c("S_dat_vero_low","S_dat_vero_fit","S_dat_vero_high","vero_coverage" ,"S")

S_O <- read.csv(paste0(root_path,"/",outpath4,"/","vero_pressure_selection_slice_2022_S.CSV"))
colnames(S_O) <- c("SO_dat_vero_low","SO_dat_vero_fit","SO_dat_vero_high","vero_coverage" ,"SO")

#custom_colors <- c("#ED0000FF", "#00468BFF")
custom_colors <- c("#5f9723", "#00468BFF")

DF1 <- merge(S,S_O,by=c("vero_coverage"),all.y = T)



p1 <-ggplot(DF1,aes(x=vero_coverage,y=S_dat_vero_fit))+
  geom_rect(aes(xmin=55, xmax=68, ymin=-Inf, ymax=Inf),fill='grey90',alpha = 0.01)+
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",linewidth=1)+
  geom_vline(aes(xintercept=25),lty='dashed',colour ="#92C5DE",linewidth=1)+
  geom_vline(aes(xintercept=50),lty='dashed',colour ="#4393C3",linewidth=1)+
  geom_vline(aes(xintercept=75),lty='dashed',colour ="#2166AC",linewidth=1)+
  geom_line(aes(x=vero_coverage,y=S_dat_vero_fit,group=1,color="Pre-Omicron period"),linewidth=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = S_dat_vero_low, ymax = S_dat_vero_high,group=1,fill="Pre-Omicron period"),alpha=0.25)+
  geom_line(aes(x=vero_coverage,y=SO_dat_vero_fit,group=1,color="Omicron period"),linewidth=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = SO_dat_vero_low, ymax = SO_dat_vero_high,group=1,fill="Omicron period"),alpha=0.25)+
  scale_x_continuous(limits=c(0,80),breaks = seq(0,75,25))+
  scale_y_continuous(limits=c(-6,3),breaks = seq(-6,3,3))+
  #labels = seq(0,100,25),expand = c(0, 0))+
  labs(x="Adjusted vaccine coverage", y= "Effect on ratio of nonsynonymous\nto synonymous divergence" )+  
  ggtitle("S")+
  #scale_x_discrete(labels=NULL)+
  scale_color_manual(name = '',values = custom_colors)+ 
  scale_fill_manual(name = '',values = custom_colors)+ 
  mytheme+
  theme(plot.margin=unit(c(0.75,1.5,0.75,1.5),'lines'))+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+ 
  theme(panel.grid=element_line(linetype=2,color='grey'))+
  guides(color=guide_legend(nrow=1),fill=guide_legend(nrow=1))+
  theme(axis.title.y.left = element_text( margin = margin(0,0.4,0,0,'cm')),
        axis.title.x = element_text(margin = margin(0.4,0,0,0,'cm')))+
  theme(legend.position = "bottom")

p1



S1 <- read.csv(paste0(root_path,"/",outpath,"/","vero_pressure_selection_slice_S1.CSV"))
colnames(S1) <- c("S1_dat_vero_low","S1_dat_vero_fit","S1_dat_vero_high","vero_coverage" ,"S1")

S1_O <- read.csv(paste0(root_path,"/",outpath4,"/","vero_pressure_selection_slice_2022_S1.CSV"))
colnames(S1_O) <- c("S1O_dat_vero_low","S1O_dat_vero_fit","S1O_dat_vero_high","vero_coverage" ,"S1O")



DF2 <- merge(S1,S1_O,by=c("vero_coverage"),all.y = T)

p2 <-ggplot(DF2,aes(x=vero_coverage,y=S1_dat_vero_fit))+
  geom_rect(aes(xmin=57, xmax=76, ymin=-Inf, ymax=Inf),fill='grey90',alpha = 0.01)+
  geom_vline(aes(xintercept=25),lty='dashed',colour ="#92C5DE",linewidth=1)+
  geom_vline(aes(xintercept=50),lty='dashed',colour ="#4393C3",linewidth=1)+
  geom_vline(aes(xintercept=75),lty='dashed',colour ="#2166AC",linewidth=1)+
  
  geom_line(aes(x=vero_coverage,y=S1_dat_vero_fit,group=1,color="Pre-Omicron period"),size=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = S1_dat_vero_low, ymax = S1_dat_vero_high,group=1,fill="Pre-Omicron"),alpha=0.25)+
  geom_line(aes(x=vero_coverage,y=S1O_dat_vero_fit,group=1,color="Omicron period"),size=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = S1O_dat_vero_low, ymax = S1O_dat_vero_high,group=1,fill="Omicron"),alpha=0.25)+
  scale_x_continuous(limits=c(0,80),breaks = seq(0,75,25))+
  scale_y_continuous(limits=c(-6,3),breaks = seq(-6,3,3))+
  #labels = seq(0,100,25),expand = c(0, 0))+
  labs(x="Adjusted vaccine coverage", y= "Effect on ratio of nonsynonymous\nto synonymous divergence" )+  
  ggtitle("S1")+
  #scale_x_discrete(labels=NULL)+
  scale_color_manual(name = '',values = custom_colors)+ 
  scale_fill_manual(name = '',values = custom_colors)+ 
  mytheme+
  theme(panel.grid=element_line(linetype=2,color='grey'))+
  guides(color=guide_legend(nrow=1),fill=guide_legend(nrow=1))+
  theme(legend.position = "bottom")+
  theme(plot.margin=unit(c(0.75,1.5,0.75,1.5),'lines'))+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+ 
  theme(axis.title.y.left = element_text( margin = margin(0,0.4,0,0,'cm')),
        axis.title.x = element_text(margin = margin(0.4,0,0,0,'cm')))+
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",size=1)
p2



S2 <- read.csv(paste0(root_path,"/",outpath,"/","vero_pressure_selection_slice_S2.CSV"))
colnames(S2) <- c("S2_dat_vero_low","S2_dat_vero_fit","S2_dat_vero_high","vero_coverage" ,"S2")

S2_O <- read.csv(paste0(root_path,"/",outpath4,"/","vero_pressure_selection_slice_2022_S2.CSV"))
colnames(S2_O) <- c("S2O_dat_vero_low","S2O_dat_vero_fit","S2O_dat_vero_high","vero_coverage" ,"S2O")


DF3 <- merge(S2,S2_O,by=c("vero_coverage"),all.y = T)

p3 <-ggplot(DF3,aes(x=vero_coverage,y=S2_dat_vero_fit))+
  geom_rect(aes(xmin=50, xmax=68, ymin=-Inf, ymax=Inf),fill='grey90',alpha = 0.01)+
  geom_vline(aes(xintercept=25),lty='dashed',colour ="#92C5DE",linewidth=1)+
  geom_vline(aes(xintercept=50),lty='dashed',colour ="#4393C3",linewidth=1)+
  geom_vline(aes(xintercept=75),lty='dashed',colour ="#2166AC",linewidth=1)+
  
  geom_line(aes(x=vero_coverage,y=S2_dat_vero_fit,group=1,color="Pre-Omicron period"),size=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = S2_dat_vero_low, ymax = S2_dat_vero_high,group=1,fill="Pre-Omicron"),alpha=0.25)+
  scale_x_continuous(limits=c(0,80),breaks = seq(0,75,25))+
  geom_line(aes(x=vero_coverage,y=S2O_dat_vero_fit,group=1,color="Omicron period"),size=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = S2O_dat_vero_low, ymax = S2O_dat_vero_high,group=1,fill="Omicron"),alpha=0.25)+
  scale_y_continuous(limits=c(-2,1),breaks = seq(-2,1,1))+
  #labels = seq(0,100,25),expand = c(0, 0))+
  labs(x="Adjusted vaccine coverage", y= "Effect on ratio of nonsynonymous\nto synonymous divergence" )+  
  ggtitle("S2")+
  #scale_x_discrete(labels=NULL)+
  scale_color_manual(name = '',values = custom_colors)+ 
  scale_fill_manual(name = '',values = custom_colors)+ 
  mytheme+
  theme(panel.grid=element_line(linetype=2,color='grey'))+
  guides(color=guide_legend(nrow=1),fill=guide_legend(nrow=1))+
  theme(legend.position = "bottom")+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+ 
  #theme(legend.position = "none")+
  theme(plot.margin=unit(c(0.75,1.5,0.75,1.5),'lines'))+
  theme(axis.title.y.left = element_text( margin = margin(0,0.4,0,0,'cm')),
        axis.title.x = element_text(margin = margin(0.4,0,0,0,'cm')))+
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",size=1)
p3




p<-ggarrange(p1, p2, p3,ncol =3, nrow = 1,widths = c(1,1,1),heights = c(1,1,1),
             labels = c("A","B","C"), align = "hv",
             legend = "bottom", common.legend = T,
             font.label = list(size = 24, face = "bold"))
p

ggsave(p,filename = file.path(outpath0,"figure_2.pdf"),
       width = 18,height =6.2)
