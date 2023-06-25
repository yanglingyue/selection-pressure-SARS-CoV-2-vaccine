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
dat1 <-read.csv("01_data/input_dnds_subsample_low_200_BV_Nextstrain.csv")# Load data

##Select the time period completely dominated by the Omicron variants.
y1 <- which(dat1$month=="2022-01"|dat1$month=="2022-02"|dat1$month=="2022-03"|dat1$month=="2022-04"|
              dat1$month=="2022-05"|dat1$month=="2022-06"|dat1$month=="2022-07"|dat1$month=="2022-08"|
              dat1$month=="2022-09") 
dat1 <-dat1[y1,]

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

outpath <-  "02_model-output/all_vaccine_202201_202209"
dir.create(outpath)

outpath1 <-  paste0(outpath,"/","result",sep="")
dir.create(outpath1)



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


write.csv(result_vero,file=paste0(outpath1,"/","vero_pressure_selection_error_bar_",gene_type[g],".csv",sep=""),row.names =F)

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
write.csv(result,file=paste0(outpath,"/","vero_pressure_selection_slice_2022_",gene_type[g],".csv",sep=""),row.names =F)



}


setwd("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/selection-vaccine-code-v8")
rm(list = ls())

library(data.table);library(corrplot);library(ggcorrplot);library(RColorBrewer)
library(ggpubr)


mytheme <- theme_minimal()+
  theme(
    #panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    #panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    # plot.title=element_text(size=12, vjust = 0.5, hjust = 0.5),
    legend.position = c("bottom"),
    # legend.key.size = unit(0.25,"cm"),
    #legend.key.width = unit(0.2,"cm"),
    #legend.key.height  = unit(0.5,"cm"),
    legend.background = element_rect(fill=NA, size=0,color=NA),
    legend.text=element_text(size=18),
    legend.title=element_text(face="bold",size=18),
    axis.line.y=element_line(linetype=1,color='black',linewidth=0.5),
    axis.line.x=element_line(linetype=1,color='black',linewidth=0.5),
    axis.ticks = element_line(linetype=2,color='black'),
    panel.grid=element_line(linetype=2,color='grey'),
    plot.title = element_text(hjust = 0.5,size=20, face = "bold"),
    axis.title.y.left = element_text(size = 18,color="black",face = "bold",vjust=0),
    axis.title.y.right =element_text(size = 18,color="black",vjust=0,face = "bold",angle=90),
    axis.title.x = element_text(size = 18, color="black",vjust = 0,face = "bold"),
    axis.text.y.left = element_text(size = 18,color="black",vjust = 0.5,angle = 0),
    axis.text.y.right = element_text(size = 18,color="black",vjust = 0.5,angle = 0),
    axis.text.x = element_text(size = 18, color="black",vjust = 0.5, angle = 0),
    axis.ticks.length=unit(0.15,"cm"),
    axis.ticks.y=element_line(color="black",linewidth=.5),
    axis.ticks.x=element_line(color="black",linewidth=.5),
    plot.margin=unit(c(2,0.25,1,0.25),'lines'),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
  )


S <- read.csv("02_model-output/all_vaccine_202003_202110/vero_pressure_selection_slice_S.CSV")
colnames(S) <- c("S_dat_vero_low","S_dat_vero_fit","S_dat_vero_high","vero_coverage" ,"S")

S_O <- read.csv("02_model-output/all_vaccine_202201_202209/vero_pressure_selection_slice_2022_S.CSV")
colnames(S_O) <- c("SO_dat_vero_low","SO_dat_vero_fit","SO_dat_vero_high","vero_coverage" ,"SO")

custom_colors <- c( "#00468BFF","#ED0000FF")

DF1 <- merge(S,S_O,by=c("vero_coverage"),all.y = T)

p1 <-ggplot(DF1,aes(x=vero_coverage,y=S_dat_vero_fit))+
  geom_rect(aes(xmin=55, xmax=68, ymin=-Inf, ymax=Inf),fill='grey90',alpha = 0.01)+
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",linewidth=1)+
  geom_line(aes(x=vero_coverage,y=S_dat_vero_fit,group=1,color="Pre-Omicron period"),linewidth=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = S_dat_vero_low, ymax = S_dat_vero_high,group=1,fill="Pre-Omicron period"),alpha=0.25)+
  geom_line(aes(x=vero_coverage,y=SO_dat_vero_fit,group=1,color="Omicron period"),linewidth=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = SO_dat_vero_low, ymax = SO_dat_vero_high,group=1,fill="Omicron period"),alpha=0.25)+
  scale_x_continuous(limits=c(0,80),breaks = seq(0,75,25))+
  scale_y_continuous(limits=c(-6,3),breaks = seq(-6,3,3))+
  #labels = seq(0,100,25),expand = c(0, 0))+
  labs(x="", y= "Effect on ratio of nonsynonymous\nto synonymous divergence" )+  
  ggtitle("S")+
  #scale_x_discrete(labels=NULL)+
  scale_color_manual(name = '',values = custom_colors)+ 
  scale_fill_manual(name = '',values = custom_colors)+ 
  mytheme+
  theme(panel.grid=element_line(linetype=2,color='grey'))+
  guides(color=guide_legend(nrow=1),fill=guide_legend(nrow=1))+
  theme(legend.position = "bottom")
#geom_vline(aes(xintercept=55),lty='solid',colour ="grey75",size=0.5)+
#geom_vline(aes(xintercept=68),lty='solid',colour ="grey75",size=0.5)+
#geom_vline(aes(xintercept=75),lty='solid',colour ="grey75",size=0.5)+
p1



S1 <- read.csv("02_model-output/all_vaccine_202003_202110/vero_pressure_selection_slice_S1.CSV")
colnames(S1) <- c("S1_dat_vero_low","S1_dat_vero_fit","S1_dat_vero_high","vero_coverage" ,"S1")

S1_O <- read.csv("02_model-output/all_vaccine_202201_202209/vero_pressure_selection_slice_2022_S1.CSV")
colnames(S1_O) <- c("S1O_dat_vero_low","S1O_dat_vero_fit","S1O_dat_vero_high","vero_coverage" ,"S1O")



DF2 <- merge(S1,S1_O,by=c("vero_coverage"),all.y = T)

p2 <-ggplot(DF2,aes(x=vero_coverage,y=S1_dat_vero_fit))+
  geom_rect(aes(xmin=57, xmax=76, ymin=-Inf, ymax=Inf),fill='grey90',alpha = 0.01)+
  
  geom_line(aes(x=vero_coverage,y=S1_dat_vero_fit,group=1,color="Pre-Omicron period"),size=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = S1_dat_vero_low, ymax = S1_dat_vero_high,group=1,fill="Pre-Omicron period"),alpha=0.25)+
  geom_line(aes(x=vero_coverage,y=S1O_dat_vero_fit,group=1,color="Omicron period"),size=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = S1O_dat_vero_low, ymax = S1O_dat_vero_high,group=1,fill="Omicron period"),alpha=0.25)+
  scale_x_continuous(limits=c(0,80),breaks = seq(0,75,25))+
  scale_y_continuous(limits=c(-6,3),breaks = seq(-6,3,3))+
  #labels = seq(0,100,25),expand = c(0, 0))+
  labs(x="Adjusted vaccine coverage", y= "" )+  
  ggtitle("S1")+
  #scale_x_discrete(labels=NULL)+
  scale_color_manual(name = '',values = custom_colors)+ 
  scale_fill_manual(name = '',values = custom_colors)+ 
  mytheme+
  theme(panel.grid=element_line(linetype=2,color='grey'))+
  guides(color=guide_legend(nrow=1),fill=guide_legend(nrow=1))+
  theme(legend.position = "bottom")+
  #theme(legend.position = "none")+
  #geom_vline(aes(xintercept=25),lty='solid',colour ="grey75",size=0.5)+
  #geom_vline(aes(xintercept=50),lty='solid',colour ="grey75",size=0.5)+
  #geom_vline(aes(xintercept=75),lty='solid',colour ="grey75",size=0.5)+
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",size=1)
p2



S2 <- read.csv("02_model-output/all_vaccine_202003_202110/vero_pressure_selection_slice_S2.CSV")
colnames(S2) <- c("S2_dat_vero_low","S2_dat_vero_fit","S2_dat_vero_high","vero_coverage" ,"S2")

S2_O <- read.csv("02_model-output/all_vaccine_202201_202209/vero_pressure_selection_slice_2022_S2.CSV")
colnames(S2_O) <- c("S2O_dat_vero_low","S2O_dat_vero_fit","S2O_dat_vero_high","vero_coverage" ,"S2O")


DF3 <- merge(S2,S2_O,by=c("vero_coverage"),all.y = T)

p3 <-ggplot(DF3,aes(x=vero_coverage,y=S2_dat_vero_fit))+
  geom_rect(aes(xmin=50, xmax=68, ymin=-Inf, ymax=Inf),fill='grey90',alpha = 0.01)+
  
  geom_line(aes(x=vero_coverage,y=S2_dat_vero_fit,group=1,color="Pre-Omicron period"),size=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = S2_dat_vero_low, ymax = S2_dat_vero_high,group=1,fill="Pre-Omicron period"),alpha=0.25)+
  scale_x_continuous(limits=c(0,80),breaks = seq(0,75,25))+
  geom_line(aes(x=vero_coverage,y=S2O_dat_vero_fit,group=1,color="Omicron period"),size=1,stat="identity")+
  geom_ribbon(aes(x=vero_coverage,ymin = S2O_dat_vero_low, ymax = S2O_dat_vero_high,group=1,fill="Omicron period"),alpha=0.25)+
  scale_y_continuous(limits=c(-1.4,0.7),breaks = seq(-1.4,0.7,0.7))+
  #labels = seq(0,100,25),expand = c(0, 0))+
  labs(x="", y= "" )+  
  ggtitle("S2")+
  #scale_x_discrete(labels=NULL)+
  scale_color_manual(name = '',values = custom_colors)+ 
  scale_fill_manual(name = '',values = custom_colors)+ 
  mytheme+
  theme(panel.grid=element_line(linetype=2,color='grey'))+
  guides(color=guide_legend(nrow=1),fill=guide_legend(nrow=1))+
  theme(legend.position = "bottom")+
  #theme(legend.position = "none")+
  #geom_vline(aes(xintercept=25),lty='solid',colour ="grey75",size=0.5)+
  #geom_vline(aes(xintercept=50),lty='solid',colour ="grey75",size=0.5)+
  #geom_vline(aes(xintercept=75),lty='solid',colour ="grey75",size=0.5)+
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",size=1)
p3




p<-ggarrange(p1, p2, p3,ncol =3, nrow = 1,widths = c(1,1,1),heights = c(1,1,1),
             labels = c("a","b","c"),  vjust=1.75,hjust= -1.5,align = "hv",
             legend = "bottom", common.legend = T,
             font.label = list(size = 24, face = "bold"))
p

ggsave(p,filename = file.path("03_fig-main","fig2-202003-202209.pdf"),
       width = 16,height =6)


