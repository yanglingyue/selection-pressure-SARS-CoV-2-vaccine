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
outpath0 <-  "04_extended_analysis/plot_output/Figure_S11"
dir.create(outpath0)
outpath <-  paste0(outpath0,"/","all_vaccine_202003_202110_IHME")
dir.create(outpath)
outpath1 <-  paste0(outpath,"/","vero_result",sep="")
dir.create(outpath1)


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
N <- log(data$cases_IHME_month) # Log transformation of the number of monthly new estimated infection from IHME
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



#Plot
setwd(paste0(root_path,"/",outpath1,sep=""))
#rm(list = ls())
datapath <- paste0(root_path,"/",outpath1,sep="")
datafile <- list.files(datapath, pattern = "*.csv$", full.names = TRUE)
files=list()
read.txt <- function(x) {
  des <- readLines(x)                   
  return(paste(des, collapse = ""))     
}
review <- lapply(datafile, read.csv)
docname <- list.files(datapath, pattern = "*.csv$")
myfiles = do.call(rbind, lapply(docname, function(x) read.csv(x, stringsAsFactors = FALSE)))
D1_vero <- myfiles




D1_vero$gene <-factor(D1_vero$gene,levels= c('S','S1','S2','E','M','N','ORF1a','ORF1b','ORF3a','ORF6','ORF7a','ORF7b','ORF8','ORF9b',
                                             ordered=TRUE))

D1_vero$coverage <-factor(D1_vero$coverage,levels= c('Low coverage (25)','Medium coverage (50)','High coverage (75)',
                                                     ordered=TRUE))
D1_vero <- na.omit(D1_vero)


p1  <- ggplot(D1_vero,aes(x=gene, y=vero_fit,color=coverage,group=coverage))+ 
  geom_point(aes(x=gene,y=vero_fit,fill=coverage),position=position_dodge(0.75),size=7,alpha=1,shape=21,color="black")+
  geom_errorbar(aes(x=gene, y=vero_fit,ymin=vero_low,ymax=vero_high),
                width=0.2,stat="identity",position = position_dodge(0.75),lwd=0.6,color="black")+
  scale_y_continuous(limits = c(-8,4),breaks = seq(-8,4,4))+
  labs(x="", y= "Effect on ratio of nonsynonymous\nto synonymous divergence" )+  
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",linewidth=0.5)+
  scale_color_manual(name = 'Adjusted vaccine coverage',values = c("#92C5DE","#4393C3","#2166AC"))+
  scale_fill_manual(name = 'Adjusted vaccine coverage',values = c("#92C5DE","#4393C3","#2166AC"))+ 
  mytheme+
  ggtitle("Adjusted vaccine coverage")+
  theme(strip.text = element_text(face="bold",size=18))+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+ 
  #theme(legend.position = "none")+
  guides(color=guide_legend(nrow=1))
#ggtitle("Correlation between Ka/Ks and international travel")

p1

ggsave(p1,filename = file.path(root_path,"/",outpath0,"/","figure_S11.pdf"),
       width = 15.5,height =6.5)



##Perform a significance test 
##Determining if there is a significant difference in associations between adjusted vaccine coverage and selection pressure on S and non-S proteins.
setwd(paste0(root_path,"/",outpath,sep=""))
datapath <- paste0(root_path,"/",outpath,sep="")
datafile <- list.files(datapath, pattern = "*.csv$", full.names = TRUE)
files=list()
read.txt <- function(x) {
  des <- readLines(x)                   
  return(paste(des, collapse = ""))     
}
review <- lapply(datafile, read.csv)
docname <- list.files(datapath, pattern = "*.csv$")
myfiles = do.call(rbind, lapply(docname, function(x) read.csv(x, stringsAsFactors = FALSE)))
D1_vero_all <- myfiles

D1_vero_all$gene <-factor(D1_vero_all$gene,levels= c('S','S1','S2','E','M','N','ORF1a','ORF1b','ORF3a',
                                                     'ORF6','ORF7a','ORF7b','ORF8','ORF9b',ordered=TRUE))
D1_vero_all <- na.omit(D1_vero_all)

S_ind1 <- which(D1_vero_all$gene=="S"|D1_vero_all$gene=="S1"|D1_vero_all$gene=="S2")

S_vero <- D1_vero_all[S_ind1,]
non_S_vero <- D1_vero_all[-S_ind1,]

t.test(S_vero$dat_vero_fit , non_S_vero$dat_vero_fit, var.equal = T)

