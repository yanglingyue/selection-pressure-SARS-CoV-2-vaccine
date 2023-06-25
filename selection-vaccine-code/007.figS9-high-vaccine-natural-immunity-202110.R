rm(list = ls())
library(ggplot2);library(ggsci);library(tidyr);library(car)
library(corrplot);library(ggcorrplot);library(RColorBrewer);library(tsModel)
library(naniar);library(mgcv);library(MASS);library(Cairo)
library(lubridate);library(dlnm);library(splines);library(showtext)
library(ggpubr);library(dplyr)

showtext_auto()
font_add("ArialN","arial.ttf")
font_families()

setwd("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/selection-vaccine-code")
dat1 <-read.csv("01_data/input_dnds_subsample_low_200_BV_Nextstrain_high.csv")# Load data

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

##Assign a unique index to each country.
dat1 <- dat1 %>%
  group_by(country) %>%
  mutate(country_index = cur_group_id())
dat1$country_index <- as.factor(dat1$country_index)



dat1$gene <- trimws(dat1$gene)
gene_type <- unique(dat1$gene)
s=length(gene_type)
print(gene_type)
outpath <-  "02_model-output/all_vaccine_high_202003_202110"
dir.create(outpath)
outpath1 <-  paste0(outpath,"/","vero_result",sep="")
dir.create(outpath1)
outpath2 <-  "02_model-output/natural_immunity_high_202003_202110"
dir.create(outpath2)
outpath3 <-   paste0(outpath2,"/","NI_result",sep="")
dir.create(outpath3)


#Model the association between adjusted vaccine coverage and the nonsynonymous to synonymous divergence ratios

for(g in 1:s)
{
  
print(gene_type[g])
dat2= subset(dat1,dat1$gene== gene_type[g])
#dat2= subset(dat1,dat1$gene== "S")
data <- dat2

source('001.load-data-cb-function.R')

# set data for models of growth rate of cases and lineage diversity
Y1  <- data$dn_ds_mean # response variable (growth rate of cases)
P  <- data$NPI_NV # forpublic health and social measure (PHSM)
V <- data$FW_BW_P # for adjusted vaccine coverage
N <- log(data$new_cases_month) # Log transformation of the number of new cases per month 
C <- data$country_index # for country interaction with month random effect
M <- data$month # The month index of the dataset 
NI <- data$NI_IHME_P # for natural immunity
df <- data.frame(Y1, P, V, N, C, M, NI)
df[which(df$N=="-Inf"),"N"] <- NA

#A <- round(max(df$V,na.rm = T),2)*10
A <- 750
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

country1 <- read.csv("04-sup-data/country_less_40_low_ni.csv")
country1 <- country1[,c("country","Vero_AD")]
colnames(country1) <- c("country","vero_prop")
country1$vero_prop <- NA
dat1_ni <- merge(dat1,country1,by=c("country"),all.y = T)


for(g in 1:s)
{
print(gene_type[g])
dat2= subset(dat1_ni,dat1_ni$gene== gene_type[g])
#dat2= subset(dat1,dat1$gene== "S")
data <- dat2  
source('001.load-data-cb-function.R')
  
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
#result$Lag <- rownames(result)
result$coverage1 <- "25"
result$coverage <- "Low coverage (25)"
result$gene <- gene_type[g]
result_NI_25 <- result

NI_low <- cp_NW$cumlow["50","lag3"]
NI_fit <- cp_NW$cumfit["50","lag3"]
NI_high <- cp_NW$cumhigh["50","lag3"]
result <- cbind(NI_low,NI_fit,NI_high)
result <- as.data.frame(result)
#result$Lag <- rownames(result)
result$coverage1 <- "50"
result$coverage <- "Medium coverage (50)"
result$gene <- gene_type[g]
result_NI_50 <- result

result <- matrix(NA, nrow = 1,ncol = 3)
result <- as.data.frame(result)
colnames(result) <- c("NI_low","NI_fit","NI_high")
result <- as.data.frame(result)
#result$Lag <- rownames(result)
result$coverage1 <- "75"
result$coverage <- "High coverage (75)"
result$gene <- gene_type[g]
result_NI_75 <- result


result_NI <- rbind(result_NI_25,result_NI_50,result_NI_75)

write.csv(result_NI,file=paste0(outpath3,"/natural_immunity_pressure_selection_error_bar_",gene_type[g],".csv"),row.names =F)


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



col <- "#42B540B2"
tcol <- do.call(rgb, c(as.list(col2rgb(col)), alpha = 255/3, max = 255))
        
png(file=paste0(outpath2,"/fig_natural_immunity_slice_",gene_type[g],".png"),width=750,height=750,res=110)
par(mfrow = c(1, 1), oma = c(0.5, 0.5, 3, 0.5))
        
plot(cp_NW,"slices",lag=3,cumul=T,
    col=col, ci.arg=list(col = tcol),mgp=c(3,1,0),
    cex.axis=1.6,cex.lab=1.6,cex.main=1.6,lwd=3,family = "ArialN",
    ci.level=0.95,col=2,xlab="Natural immunity",
    ylab="Effect on selection pressure")
        
mtext(paste0(gene_type[g]), side = 3, font=2,family = "ArialN",cex=1.4 ,line = -1, outer = T)
dev.off()

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
D1_vero <- myfiles


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
    axis.title.x = element_text(size = 18, color="black",vjust = 0,face = "bold"),
    axis.text.y.left = element_text(size = 18,color="black",vjust = 0.5,angle = 0),
    axis.text.y.right = element_text(size = 18,color="black",vjust = 0.5,angle = 0),
    axis.text.x = element_text(size = 18, color="black",vjust = 0.5, angle = 0),
    axis.ticks.length=unit(0.15,"cm"),
    axis.ticks.y=element_line(color="black",linewidth=.5),
    axis.ticks.x=element_line(color="black",linewidth=.5),
    plot.margin=unit(c(1,0.75,1,0.75),'lines'),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
  )


D1_vero$gene <-factor(D1_vero$gene,levels= c('S','S1','S2','E','M','N','ORF1a','ORF1b','ORF3a','ORF6','ORF7a','ORF7b','ORF8','ORF9b',
                                         ordered=TRUE))

D1_vero$coverage <-factor(D1_vero$coverage,levels= c('Low coverage (25)','Medium coverage (50)','High coverage (75)',
                                                 ordered=TRUE))
D1_vero <- na.omit(D1_vero)


p1  <- ggplot(D1_vero,aes(x=gene, y=vero_fit,color=coverage,group=coverage))+ 
  geom_point(aes(x=gene,y=vero_fit,color=coverage),position=position_dodge(0.5),size=4.5,alpha=1)+
  geom_errorbar(aes(x=gene, y=vero_fit,ymin=vero_low,ymax=vero_high),
                width=0.2,stat="identity",position = position_dodge(0.5),lwd=0.6,color="black")+
  scale_y_continuous(limits = c(-12,6),breaks = seq(-12,6,6))+
  labs(x="", y= "Effect on ratio of nonsynonymous\nto synonymous divergence" )+  
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",linewidth=0.5)+
  scale_color_manual(name = 'Adjusted vaccine coverage',values = c("#92C5DE","#4393C3","#2166AC"))+ 
  mytheme+
  ggtitle("Adjusted vaccine coverage")+
  theme(strip.text = element_text(face="bold",size=18))+
  #theme(legend.position = "none")+
  guides(color=guide_legend(nrow=1))
#ggtitle("Correlation between Ka/Ks and international travel")

p1


setwd(paste0("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/selection-vaccine-code-v8","/",outpath3,sep=""))
datapath <- paste0("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/selection-vaccine-code-v8","/",outpath3,sep="")
datafile <- list.files(datapath, pattern = "*.csv$", full.names = TRUE)
files=list()
read.txt <- function(x) {
  des <- readLines(x)                   
  return(paste(des, collapse = ""))     
}
review <- lapply(datafile, read.csv)
docname <- list.files(datapath, pattern = "*.csv$")
myfiles = do.call(rbind, lapply(docname, function(x) read.csv(x, stringsAsFactors = FALSE)))
D1_Ni <- myfiles

D1_Ni$gene <-factor(D1_Ni$gene,levels= c('S','S1','S2','E','M','N','ORF1a','ORF1b','ORF3a','ORF6','ORF7a','ORF7b','ORF8','ORF9b',
                                         ordered=TRUE))

D1_Ni$coverage <-factor(D1_Ni$coverage,levels= c('Low coverage (25)','Medium coverage (50)','High coverage (75)',
                                                 ordered=TRUE))

D1_Ni <- na.omit(D1_Ni)


p2  <- ggplot(D1_Ni,aes(x=gene, y=NI_fit,color=coverage,group=coverage))+ 
  geom_point(aes(x=gene,y=NI_fit,color=coverage),position=position_dodge(0.5),size=4.5,alpha=1)+
  geom_errorbar(aes(x=gene, y=NI_fit,ymin=NI_low,ymax=NI_high),
                width=0.2,stat="identity",position = position_dodge(0.5),lwd=0.6,color="black")+
  scale_y_continuous(limits = c(-2,6),breaks = seq(-2,6,2))+
  labs(x="", y= "Effect on ratio of nonsynonymous\nto synonymous divergence" )+  
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",linewidth=0.5)+
  scale_color_manual(name = 'Natural immunnity coverage',values = c("#F4A582","#D6604D","#B2182B" ))+ 
  mytheme+
  ggtitle("Natural immunity")+
  theme(strip.text = element_text(face="bold",size=18))+
  guides(color=guide_legend(nrow=1))

#ggtitle("Correlation between Ka/Ks and international travel")

p2



p<-ggarrange(p1, p2,ncol =1, nrow = 2,widths = c(1,1),heights = c(1,1),
             labels = c("a","b"),  vjust=1.1,align = "hv",
             #legend = "bottom",
             font.label = list(size = 24, face = "bold"))
p

ggsave(p,filename = file.path("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/selection-vaccine-code-v8/03_fig-main"
                              ,"figS9-high-all-vaccine-202003-202110.pdf"),
       width = 13.5,height =12)


##Perform a significance test 
##Determining if there is a significant difference in associations between adjusted vaccine coverage and selection pressure on S and non-S proteins.
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

S_ind1 <- which(D1_vero_all$gene=="S"|D1_vero_all$gene=="S1"|D1_vero_all$gene=="S2")

S_vero <- D1_vero_all[S_ind1,]
non_S_vero <- D1_vero_all[-S_ind1,]

t.test(S_vero$dat_vero_fit , non_S_vero$dat_vero_fit, var.equal = T)


##Perform a significance test 
##Determining if there is a significant difference in associations between natural immunity and selection pressure on S and non-S proteins.
setwd(paste0("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/selection-vaccine-code-v8","/",outpath2,sep=""))
datapath <- paste0("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/selection-vaccine-code-v8","/",outpath2,sep="")
datafile <- list.files(datapath, pattern = "*.csv$", full.names = TRUE)
files=list()
read.txt <- function(x) {
  des <- readLines(x)                   
  return(paste(des, collapse = ""))     
}
review <- lapply(datafile, read.csv)
docname <- list.files(datapath, pattern = "*.csv$")
myfiles = do.call(rbind, lapply(docname, function(x) read.csv(x, stringsAsFactors = FALSE)))
D1_Ni_all <- myfiles
D1_Ni_all$gene <-factor(D1_Ni_all$gene,levels= c('S','S1','S2','E','M','N','ORF1a','ORF1b','ORF3a',
                                                     'ORF6','ORF7a','ORF7b','ORF8','ORF9b',ordered=TRUE))
D1_Ni_all <- na.omit(D1_Ni_all)

S_ind2 <- which(D1_Ni_all$gene=="S"|D1_Ni_all$gene=="S1"|D1_Ni_all$gene=="S2")
S_Ni <- D1_Ni_all[S_ind2,]
non_S_Ni <- D1_Ni_all[-S_ind2,]

t.test(S_Ni$dat_NI_fit, non_S_Ni$dat_NI_fit, var.equal = T)


