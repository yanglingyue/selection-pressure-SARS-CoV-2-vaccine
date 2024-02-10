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

#outpath <-  "02_model-output/all_vaccine_202003_202110"
outpath0 <- "03_output/figure_2_3"
dir.create(outpath0)
outpath <-  paste0(outpath0,"/","muts_nonsyn_vaccine_202003_202110")
dir.create(outpath)
outpath1 <-  paste0(outpath,"/","vero_result",sep="")
dir.create(outpath1)
outpath2 <-  paste0(outpath0,"/","muts_nonsyn_NI_202003_202110")
dir.create(outpath2)
outpath3 <-   paste0(outpath2,"/","NI_result",sep="")
dir.create(outpath3)
outpath4 <-  paste0(outpath0,"/","muts_nonsyn_vaccine_202201_202209")
dir.create(outpath4)



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


model1 <-  gam(Y1 ~ npi_cb+FBW_cb+travel_cb+N+s(C, bs = "re")+P:V,
              method = "REML",data=df)
summary(model1)
AIC(model1)
cp_FW <- crosspred(FBW_cb, model1, cen=0,by=0.1,at=0:A/10,bylag=0.1,cumul=TRUE) #slices

#assign(paste("model_vero_", gene_type[g], sep=""),model1)
#assign(paste("cp_cases_", gene_type[g], sep=""),cp_cases)
#assign(paste("cp_FW_vero_", gene_type[g], sep=""),cp_FW)

# Extract estimates for the number of nonsynonymous mutation model when the adjusted vaccine coverage is 25%
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


# Extract estimates for the number of nonsynonymous mutation model when the adjusted vaccine coverage is 50%
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

# Extract estimates for the number of nonsynonymous mutation model when the adjusted vaccine coverage is 75%
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
write.csv(result_vero,file=paste0(outpath1,"/vero_muts_nonsyn_error_bar_",gene_type[g],".csv"),row.names =F)


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
write.csv(result_vero_all,file=paste0(outpath,"/vero_muts_nonsyn_slice_",gene_type[g],".csv"),row.names =F)




}





#Model the association between natural immunity coverage and number of nonsynonymous mutations
full_vero_2110 <- dat1[which(dat1$gene=="S1"&dat1$month=="2021-10"),]

#Filter out countries with a full-dose COVID-19 vaccine coverage rate below 25% as of September 2022.
country_low <- full_vero_2110[which(full_vero_2110$F_Vero<=25),"country"]
country_low$country <- as.character(country_low$country)
dat1_ni<- dat1[dat1$country %in% country_low$country, ]


for(g in 1:s)
{
  print(gene_type[g])
  data= subset(dat1_ni,dat1_ni$gene== gene_type[g])
  #dat2= subset(dat1,dat1$gene== "S")
  source('load-data-cb-function.R')
  
  # set data for models of growth rate of cases and lineage diversity
  #Y1  <- data$dnds_mean_focal # response variable for ratio of nonsynonymous to synonymous divergence
  Y1  <- data$muts_nonsyn_focal # response variable for number of nonsynonymous mutations
  #Y1  <- data$fitness_focal # response variable for mutational fitness
  P  <- data$NPI_NV # forpublic health and social measure (PHSM)
  C <- data$country_index # for country interaction with month random effect
  M <- data$month # The month index of the dataset 
  NI <- data$NI_IHME_P_omicron # for natural immunity
  
  df <- data.frame(Y1, P, C, M, NI)
  
  
  B <- (round(max(df$NI,na.rm=T),2))*10
  #B <- 750
  
  model2 <-  gam(Y1 ~ npi_cb+NIW_cb+travel_cb+s(C, bs = "re"),
                 method = "REML",data=df)
  summary(model2)
  AIC(model2)
  cp_NW <- crosspred(NIW_cb, model2, cen=0,by=0.1,at=0:B/10,bylag=0.1,cumul=TRUE) #slices

  
  assign(paste("model_NI_",gene_type[g], sep=""),model2)
  assign(paste("cp_NW_", gene_type[g], sep=""),cp_NW)
  
  
  
  
  # Extract estimates for the number of nonsynonymous mutation model when the natural immunity is 25%
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
  
  # Extract estimates for the number of nonsynonymous mutation model when the natural immunity is 50%
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
  
  # Extract estimates for the number of nonsynonymous mutation model when the natural immunity is 75%
  result <- matrix(NA, nrow = 1,ncol = 3)
  result <- as.data.frame(result)
  colnames(result) <- c("NI_low","NI_fit","NI_high")
  result <- as.data.frame(result)
  result$Lag <- rownames(result)
  result$coverage1 <- "75"
  result$coverage <- "High coverage (75)"
  result$gene <- gene_type[g]
  result_NI_75 <- result
  
  
  result_NI <- rbind(result_NI_25,result_NI_50,result_NI_75)
  result_NI$Lag <- "lag3"
  
  write.csv(result_NI,file=paste0(outpath3,"/NI_muts_nonsyn_error_bar_",gene_type[g],".csv"),row.names =F)
  
  # Extract the all estimates for the association between number of nonsynonymous mutations and natural immunity
  
  dat_NI_low <- as.data.frame(cp_NW$cumlow[,"lag3"])
  colnames(dat_NI_low) <- "dat_NI_low"
  dat_NI_fit <- as.data.frame(cp_NW$cumfit[,"lag3"])
  colnames(dat_NI_fit) <- "dat_NI_fit"
  dat_NI_high <- as.data.frame(cp_NW$cumhigh[,"lag3"])
  colnames(dat_NI_high) <- "dat_NI_high"
  result_NI_all <- cbind(dat_NI_low,dat_NI_fit,dat_NI_high)
  result_NI_all$coverage <- rownames(result_NI_all)
  result_NI_all <- as.data.frame(result_NI_all)
  result_NI_all$gene <- gene_type[g]
  
  write.csv(result_NI_all,file=paste0(outpath2,"/NI_muts_nonsyn_slice_",gene_type[g],".csv"),row.names =F)
  
  

  
  
        
}



#Model the association between adjusted vaccine coverage and the number of nonsynonymous mutations for Omicron
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
  
# set data for models of ratio of nonsynonymous to synonymous divergence
#Y1 <- data$dnds_mean_focal # response variable for ratio of nonsynonymous to synonymous divergence
Y1 <- data$muts_nonsyn_focal # response variable for number of nonsynonymous mutations
#Y3 <- data$fitness_focal # response variable for mutational fitness
P <- data$NPI_NV # forpublic health and social measure (PHSM)
V <- data$FW_BW_P # for adjusted vaccine coverage
N <- log(data$new_cases_month) # Log transformation of the number of new cases per month 
C <- data$country_index # for country interaction with month random effect
C1 <- data$country # for country interaction with month random effect
M <- data$month # The month index of the dataset 
NI <- data$NI_IHME_P_omicron # for natural immunity
  
df <- data.frame(Y1, P, V, N, C,C1,NI, M)
df[which(df$N=="-Inf"),"N"] <- NA
  
A <- round(max(df$V,na.rm = T),2)*10

  
model3 <-  gam(Y1 ~ npi_cb+FBW_cb+N+travel_cb+s(C, bs = "re")+P:V,method = "REML",data=df)
summary(model3)
AIC(model3)
cp_FW <- crosspred(FBW_cb, model3, cen=0,by=0.1,at=0:A/10,bylag=0.1,cumul=TRUE) #slices
#cp_NW <- crosspred(NIW_cb, model1, cen=0,by=0.1,at=0:B/10,bylag=0.1,cumul=TRUE) #slices


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
write.csv(result,file=paste0(outpath4,"/","vero_muts_nonsyn_slice_2022_",gene_type[g],".csv",sep=""),row.names =F)
  
  
  
  


  
}
