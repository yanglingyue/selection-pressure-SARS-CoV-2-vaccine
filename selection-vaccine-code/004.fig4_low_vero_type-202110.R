rm(list = ls())
library(ggplot2);library(ggsci);library(tidyr);library(car)
library(corrplot);library(ggcorrplot);library(RColorBrewer);library(tsModel)
library(naniar);library(mgcv);library(MASS);library(Cairo)
library(lubridate);library(dlnm);library(splines);library(showtext)
library(ggpubr);library(dplyr)

showtext_auto()
font_add("ArialN","arial.ttf")
font_families()

setwd("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/selection-vaccine-code-v8")
dat1 <-read.csv("01_data/input_dnds_subsample_low_200_BV_Nextstrain.csv")# Load data

##Remove January and February 2020, as the submitted sequences of most countries are unavailable.
y1 <- which(dat1$month=="2020-01"|dat1$month=="2020-02") 
dat1 <-dat1[-y1,]

##Remove November 2021 to September 2022 as the Omicron variant emerged and circulated during this period.
y2 <- which(dat1$month=="2021-11"|dat1$month=="2021-12"|dat1$month=="2022-01"|dat1$month=="2022-02"|
              dat1$month=="2022-03"|dat1$month=="2022-04"|dat1$month=="2022-05"|dat1$month=="2022-06"|
              dat1$month=="2022-07"|dat1$month=="2022-08"|dat1$month=="2022-09") 
dat1 <-dat1[-y2,]

##Remove countries with missing data for other indexes such as vaccines or PHSM (Public Health and Social Measures).
x <- which(dat1$country=="Puerto Rico"|dat1$country=="Canary Islands"|dat1$country=="French Guiana"|
             dat1$country=="Reunion"|dat1$country=="Curacao"|dat1$country=="North Macedonia"|
             dat1$country=="Slovakia")
dat1 <-dat1[-x,]


dat1$country <- as.factor(dat1$country)
dat1<- dat1[order(dat1$country,dat1$month,decreasing=F),]


dat1$gene <- trimws(dat1$gene)
gene_type <- unique(dat1$gene)
s=length(gene_type)
print(gene_type)


country1 <- read.csv("04-sup-data/country_mRNA_low.csv")
country1 <- country1[,c("country","Region")]
country2 <- read.csv("04-sup-data/country_inactivated_vaccine_low.csv")
country2 <- country2[,c("country","Region")]
country3 <- read.csv("04-sup-data/country_adenovirus_vector_low.csv")
country3 <- country3[,c("country","Region")]



dat1$gene <- trimws(dat1$gene)
gene_type <- unique(dat1$gene)
s=length(gene_type)
print(gene_type)

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


outpath <-  "02_model-output/vero_type_202003_202110"
dir.create(outpath)

outpath1 <-  paste0(outpath,"/","vero_type_result",sep="")
dir.create(outpath1)



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

  source('001.load-data-cb-function.R')
  
# set data for models of growth rate of cases and lineage diversity
Y1  <- data$dn_ds_mean # response variable (growth rate of cases)
#Y2  <- data$mf # response variable (growth rate of cases)
P  <- data$NPI_NV # forpublic health and social measure (PHSM)
V <- data$FW_BW_P # for adjusted vaccine coverage
N <- log(data$new_cases_month) # Log transformation of the number of new cases per month 
C <- data$country_index # for country interaction with month random effect
M <- data$month # The month index of the dataset 

df <- data.frame(Y1, P, V, N, C, M)

model1 <-  gam(Y1 ~ npi_cb+FBW_cb+N+travel_cb+s(C, bs = "re")+P:V,
              method = "REML",data=df)

summary(model1)
AIC(model1)

A <- (round(max(df$V,na.rm=T),2))*10

cp_FW <- crosspred(FBW_cb, model1, cen=0,by=0.1,at=0:A/10,bylag=0.1,cumul=TRUE) #slices


assign(paste("model_",vero_type[i],"_",gene_type[g], sep=""),model1)
assign(paste("cp_FW_", vero_type[i], sep=""),cp_FW)

#########fig.3A################

if (vero_type[i]=="Inactivated_vaccine"){
  


vero_low <- cp_FW$cumlow["50","lag3"]
vero_fit <- cp_FW$cumfit["50","lag3"]
vero_high <- cp_FW$cumhigh["50","lag3"]
result <- cbind(vero_low,vero_fit,vero_high)
result <- as.data.frame(result)
#result$Lag <- rownames(result)
result$Type <- vero_type[i]
result$coverage1 <- "50"
result$coverage <- "Medium coverage (50)"
result$gene <- gene_type[g]
result_vero_50 <- result


vero_low <- cp_FW$cumlow["25","lag3"]
vero_fit <- cp_FW$cumfit["25","lag3"]
vero_high <- cp_FW$cumhigh["25","lag3"]
result <- cbind(vero_low,vero_fit,vero_high)
result <- as.data.frame(result)
#result$Lag <- rownames(result)
result$Type <- vero_type[i]
result$coverage1 <- "25"
result$coverage <- "Low coverage (25)"
result$gene <- gene_type[g]
result_vero_25 <- result

result <- matrix(NA, nrow = 1,ncol = 3)
result <- as.data.frame(result)
colnames(result) <- c("vero_low","vero_fit","vero_high")
#result$Lag <- rownames(result)
result$Type <- vero_type[i]
result$coverage1 <- "75"
result$coverage <- "High coverage (75)"
result$gene <- gene_type[g]
result_vero_75 <- result

}

else{

if (vero_type[i]=="mRNA_vaccine")
{
  

  vero_low <- cp_FW$cumlow["75","lag3"]
  vero_fit <- cp_FW$cumfit["75","lag3"]
  vero_high <- cp_FW$cumhigh["75","lag3"]
  result <- cbind(vero_low,vero_fit,vero_high)
  result <- as.data.frame(result)
  #result$Lag <- rownames(result)
  result$Type <- vero_type[i]
  result$coverage1 <- "75"
  result$coverage <- "High coverage (75)"
  result$gene <- gene_type[g]
  result_vero_75 <- result 
  
 
  
  vero_low <- cp_FW$cumlow["50","lag3"]
  vero_fit <- cp_FW$cumfit["50","lag3"]
  vero_high <- cp_FW$cumhigh["50","lag3"]
  result <- cbind(vero_low,vero_fit,vero_high)
  result <- as.data.frame(result)
  #result$Lag <- rownames(result)
  result$Type <- vero_type[i]
  result$coverage1 <- "50"
  result$coverage <- "Medium coverage (50)"
  result$gene <- gene_type[g]
  result_vero_50 <- result
  
  
  vero_low <- cp_FW$cumlow["25","lag3"]
  vero_fit <- cp_FW$cumfit["25","lag3"]
  vero_high <- cp_FW$cumhigh["25","lag3"]
  result <- cbind(vero_low,vero_fit,vero_high)
  result <- as.data.frame(result)
  #result$Lag <- rownames(result)
  result$Type <- vero_type[i]
  result$coverage1 <- "25"
  result$coverage <- "Low coverage (25)"
  result$gene <- gene_type[g]
  result_vero_25 <- result
}
  
  else{
    
    result <- matrix(NA, nrow = 1,ncol = 3)
    result <- as.data.frame(result)
    colnames(result) <- c("vero_low","vero_fit","vero_high")
    #result$Lag <- rownames(result)
    result$Type <- vero_type[i]
    result$coverage1 <- "75"
    result$coverage <- "High coverage (75)"
    result$gene <- gene_type[g]
    result_vero_75 <- result
    
   
    
    result <- matrix(NA, nrow = 1,ncol = 3)
    result <- as.data.frame(result)
    colnames(result) <- c("vero_low","vero_fit","vero_high")
    result <- as.data.frame(result)
    #result$Lag <- rownames(result)
    result$Type <- vero_type[i]
    result$coverage1 <- "50"
    result$coverage <- "Medium coverage (50)"
    result$gene <- gene_type[g]
    result_vero_50 <- result
    
    
    result <- matrix(NA, nrow = 1,ncol = 3)
    result <- as.data.frame(result)
    colnames(result) <- c("vero_low","vero_fit","vero_high")
    result <- as.data.frame(result)
    #result$Lag <- rownames(result)
    result$Type <- vero_type[i]
    result$coverage1 <- "25"
    result$coverage <- "Low coverage (25)"
    result$gene <- gene_type[g]
    result_vero_25 <- result
  }
}


result_vero <- rbind(result_vero_25,result_vero_50,result_vero_75)

write.csv(result_vero,file=paste0(outpath1,"/model_selection_",gene_type[g],"_",vero_type[i],".csv"),row.names =F)



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

write.csv(result,file=paste0(outpath,"/fit_vero_type_",gene_type[g],"_",vero_type[i],".csv"),row.names =F)



}

 
  
}



#Plot
setwd(paste0("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/selection-vaccine-code-v8","/",outpath1,sep=""))
#rm(list = ls())
datapath <- paste0("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/selection-vaccine-code-v8","/",outpath1,sep="")

datafile <- list.files(datapath, pattern = "*.csv$", full.names = TRUE)
files=list()
read.txt <- function(x) {
  des <- readLines(x)                   
  return(paste(des, collapse = ""))     
}
review <- lapply(datafile, read.csv)
docname <- list.files(datapath, pattern = "*.csv$")
myfiles = do.call(rbind, lapply(docname, function(x) read.csv(x, stringsAsFactors = FALSE)))

mytheme <- theme_minimal()+
  theme(
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    # plot.title=element_text(size=12, vjust = 0.5, hjust = 0.5),
    legend.position = c("bottom"),
    # legend.key.size = unit(0.25,"cm"),
    #legend.key.width = unit(0.2,"cm"),
    #legend.key.height  = unit(0.5,"cm"),
    legend.background = element_rect(fill=NA, linewidth=0,color=NA),
    legend.text=element_text(size=18),
    legend.title=element_text(face="bold",size=18),
    axis.line.y=element_line(linetype=1,color='black',linewidth=0.5),
    axis.line.x=element_line(linetype=1,color='black',linewidth=0.5),
    axis.ticks = element_line(linetype=2,color='black'),
    panel.grid=element_line(linetype=2,color='grey'),
    plot.title = element_text(hjust = 0.5,size=20, face = "bold"),
    axis.title.y.left = element_text(size = 18,color="black",vjust=2.5,face = "bold"),
    axis.title.y.right =element_text(size = 18,color="black",vjust=0,angle=90,face = "bold"),
    axis.title.x = element_text(size = 18, color="black",vjust = 0),
    axis.text.y.left = element_text(size = 18,color="black",vjust = 0.5,angle = 0),
    axis.text.y.right = element_text(size = 18,color="black",vjust = 0.5,angle = 0),
    axis.text.x = element_text(size = 18, color="black",vjust = 0.5, angle = 0),
    axis.ticks.length=unit(0.15,"cm"),
    axis.ticks.y=element_line(color="black",linewidth=.5),
    axis.ticks.x=element_line(color="black",linewidth=.5),
    plot.margin=unit(c(0.5,0.5,0.5,0.5),'lines'),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
  )


DF_vero <- myfiles
DF_vero[which(DF_vero$Type=="Inactivated_vaccine"),"Type"] <- "Inactivated vaccine"
DF_vero[which(DF_vero$Type=="mRNA_vaccine"),"Type"] <- "mRNA vaccine"
DF_vero[which(DF_vero$Type=="Viral_vector_vaccine"),"Type"] <- "Viral vector vaccine"


#data1$gene <-factor(data1$gene,levels= c('S','RBD',
#                                         ordered=TRUE))

DF_vero$gene <-factor(DF_vero$gene,levels= c('S','S1','S2','E','M','N','ORF1a','ORF1b','ORF3a','ORF6','ORF7a','ORF7b','ORF8','ORF9b',
                                         ordered=TRUE))

DF_vero$coverage <-factor(DF_vero$coverage,levels= c('Low coverage (25)','Medium coverage (50)','High coverage (75)',
                                                 ordered=TRUE))

x <- which(is.na(DF_vero$gene)|is.na(DF_vero$coverage))
D1 <- DF_vero[-x,]
D1_mrna= subset(D1,D1[,"Type"]== "mRNA vaccine")
D1_ad= subset(D1,D1[,"Type"]== "Viral vector vaccine")
D1_in= subset(D1,D1[,"Type"]== "Inactivated vaccine")


p1  <- ggplot(D1_mrna,aes(x=gene, y=vero_fit,color=coverage,group=coverage))+ 
  geom_point(aes(x=gene,y=vero_fit,color=coverage),position=position_dodge(0.5),size=4.5,alpha=0.8)+
  geom_errorbar(aes(x=gene, y=vero_fit,ymin=vero_low,ymax=vero_high),
                width=0.2,stat="identity",position = position_dodge(0.5),lwd=0.6,color="black")+
  scale_y_continuous(limits = c(-6,6),breaks = seq(-6,6,3))+
  #coord_cartesian(ylim = c(-20, 20))+
  #facet_wrap(~Type,scales = "free",nrow=1)+
  labs(x="", y= "Effect on ratio of nonsynonymous\nto synonymous divergence" )+  
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",linewidth=0.5)+
  scale_color_manual(name = 'Adjusted vaccine coverage',values = c("#EA2027","#45aaf2","#504190"))+ 
  mytheme+
  theme(strip.background = element_blank(), strip.text = element_blank())+
  theme(strip.text = element_text(face="bold",size=18))+
  guides(color=guide_legend(nrow=1))+
  ggtitle("mRNA vaccine")

p1


p2  <- ggplot(D1_in,aes(x=gene, y=vero_fit,color=coverage,group=coverage))+ 
  geom_point(aes(x=gene,y=vero_fit,color=coverage),position=position_dodge(0.5),size=4.5,alpha=0.8)+
  geom_errorbar(aes(x=gene, y=vero_fit,ymin=vero_low,ymax=vero_high),
                width=0.2,stat="identity",position = position_dodge(0.5),lwd=0.6,color="black")+
  scale_y_continuous(limits = c(-3,9),breaks = c(-3,0,3,6,9))+
  #coord_cartesian(ylim = c(-20, 20))+
  #facet_wrap(~Type,scales = "free",nrow=1)+
  labs(x="", y= "Effect on ratio of nonsynonymous\nto synonymous divergence" )+  
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",linewidth=0.5)+
  scale_color_manual(name = 'Adjusted vaccine coverage',values = c("#EA2027","#45aaf2","#504190"))+ 
  mytheme+
  theme(strip.background = element_blank(), strip.text = element_blank())+
  theme(strip.text = element_text(face="bold",size=18))+
  guides(color=guide_legend(nrow=1))+
  ggtitle("Inactivated vaccine")
p2

p<-ggarrange(p1,p2,ncol =1, nrow = 2,
             labels = c("a","b"),  vjust=1,
             common.legend=TRUE,legend = "bottom",align = "hv",
             font.label = list(size = 24, face = "bold"))
p


ggsave(p,filename = file.path("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/selection-vaccine-code-v8/03_fig-main"
                              ,"fig4-vaccine-type-202003-202110.pdf"),
       width = 13.5,height =11)

## Load results
setwd(paste0("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/selection-vaccine-code-v8","/",outpath,sep=""))
datapath <- paste0("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/selection-vaccine-code-v8","/",outpath,sep="")
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


##Perform a significance test
##Determining if there is a significant difference in associations between adjusted vaccine coverage and selection pressure on S and non-S proteins.
##For adjusted vaccine coverage for mRNA vaccines

D1_vero_mRNA <- subset(D1_vero_all,D1_vero_all$Type=="mRNA_vaccine")
S_ind1 <- which(D1_vero_mRNA$gene=="S"|D1_vero_mRNA$gene=="S1"|D1_vero_mRNA$gene=="S2")

S_vero1 <- D1_vero_mRNA[S_ind1,]
non_S_vero1 <- D1_vero_mRNA[-S_ind1,]

t.test(S_vero1$dat_vero_fit, non_S_vero1$dat_vero_fit, var.equal = T)



##Perform a significance test
##Determining if there is a significant difference in associations between adjusted vaccine coverage and selection pressure on S and non-S proteins.
##For adjusted vaccine coverage for inactivated vaccines

D1_vero_in <- subset(D1_vero_all,D1_vero_all$Type=="Inactivated_vaccine")
S_ind2 <- which(D1_vero_in$gene=="S"|D1_vero_in$gene=="S1"|D1_vero_in$gene=="S2")

S_vero2 <- D1_vero_in[S_ind2,]
non_S_vero2 <- D1_vero_in[-S_ind2,]

t.test(S_vero2$dat_vero_fit, non_S_vero2$dat_vero_fit, var.equal = T)
