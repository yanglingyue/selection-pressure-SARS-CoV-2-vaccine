setwd("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/selection-vaccine-code-v8")
rm(list = ls())
#Sys.setenv(LANGUAGE = "en")
#Sys.setlocale(category = "LC_ALL", locale = "English_United States.1252")
library(ggplot2);library(ggsci);library(RColorBrewer);library(tsModel);
library(tidyverse);library(tidyr);library(lubridate);library(mgcv);
library(MASS);library(splines);library(dlnm);library(showtext);
library(car);library(R330);library(boot);library(dplyr)
showtext_auto()
font_add("ArialN","arial.ttf")
font_families()
dat1 <-read.csv("01_data/input_dnds_subsample_low_200_BV_Nextstrain.csv")# Load data

#y <- which(dat1$month=="2020-01"|dat1$month=="2020-02")
y1 <- which(dat1$month=="2020-01"|dat1$month=="2020-02"|
              dat1$month=="2021-11"|dat1$month=="2021-12"|dat1$month=="2022-01"|dat1$month=="2022-02"|
              dat1$month=="2022-03"|dat1$month=="2022-04"|dat1$month=="2022-05"|dat1$month=="2022-06"|
              dat1$month=="2022-07"|dat1$month=="2022-08"|dat1$month=="2022-09") 
dat1 <-dat1[-y1,]
#dat1 <-dat1[-c(y),]
x <- which(dat1$country=="Puerto Rico"|dat1$country=="Canary Islands"|dat1$country=="French Guiana"|
             dat1$country=="Reunion"|dat1$country=="Curacao"|dat1$country=="North Macedonia"|
             dat1$country=="Slovakia")
dat1 <-dat1[-x,]

dat1 <- dat1 %>%
  group_by(country) %>%
  mutate(country_index = cur_group_id())
dat1$country_index <- as.factor(dat1$country_index)

dat1$travel2 <- log(dat1$travel_route)
dat1$case2 <- log(dat1$new_cases_month)

dat1$country <- as.factor(dat1$country)
dat1<- dat1[order(dat1$country,dat1$month,decreasing=F),]
dat2= subset(dat1,dat1$gene== "S")
data <- dat2


# Set maximum lag
nlag = 3

# Creating lagged variables
# Define the monthly mean covariate lag term matrix
lag_npi <- tsModel::Lag(data$NPI_NV, group = data$country_index, k = 0:nlag)
lag_travel <- tsModel::Lag(data$travel2, group = data$country_index, k = 0:nlag)
lag_FW_vero <- tsModel::Lag(data$FW_BW_P, group = data$country_index, k = 0:nlag)
#lag_mf <- tsModel::Lag(data$mf, group = data$country, k = 0:nlag)


# Set the lag section

var <- lag_npi
npi_cb <- crossbasis(var,
                      argvar = list(fun = "lin"),
                      arglag = list(fun = "ns",df=2))
summary(npi_cb)



var <- lag_travel
travel_cb <- crossbasis(var,
                        argvar = list(fun = "lin"),
                        arglag = list(fun = "ns", df=2))
summary(travel_cb)



############adjusted vaccine coverage consider immunity waning and booster#######################################

var <- lag_FW_vero
FW_cb1 <- crossbasis(var,
                     argvar = list(fun = "lin"),
                     arglag = list(fun = "ns", df=2))
summary(FW_cb1)


var <- lag_FW_vero
FW_cb2 <- crossbasis(var,
                     argvar = list(fun = "poly",degree=2),
                     arglag = list(fun = "ns", df=2))
summary(FW_cb2)


var <- lag_FW_vero
FW_cb3 <- crossbasis(var,
                     argvar = list(fun = "poly",degree=3),
                     arglag = list(fun = "ns", df=2))
summary(FW_cb3)


var <- lag_FW_vero
FW_cb4 <- crossbasis(var,
                     argvar = list(fun = "ns",
                                   knots = equalknots(data$FW_BW_P, 1)),
                     arglag = list(fun = "ns", df=2))
summary(FW_cb4)


var <- lag_FW_vero
FW_cb5 <- crossbasis(var,
                     argvar = list(fun = "ns",
                                   knots = equalknots(data$FW_BW_P, 2)),
                     arglag = list(fun = "ns", df=2))
summary(FW_cb5)




# assign unique column names to cross-basis matrix for gam() model
# note: not necessary for glm(), gam() or glm.nb() models
#colnames(npi_cb) = paste0("npi_cb.", colnames(npi_cb))

#colnames(travel_cb) = paste0("travel_cb.", colnames(travel_cb))


#colnames(FW_cb1) = paste0("FW_cb1.", colnames(FW_cb1))
#colnames(FW_cb2) = paste0("FW_cb2.", colnames(FW_cb2))
#colnames(FW_cb3) = paste0("FW_cb3.", colnames(FW_cb3))
#colnames(FW_cb4) = paste0("FW_cb4.", colnames(FW_cb4))
#colnames(FW_cb5) = paste0("FW_cb5.", colnames(FW_cb5))


# set data for models of growth rate of cases and lineage diversity
Y1  <- data$dn_ds_mean # response variable (growth rate of cases)
#Y2  <- data$mf # response variable (Shannon's index of lineage diversity based on the subsampling strategy 1)
P  <- data$NPI_NV # forpublic health and social measure (PHSM)
V <- data$FW_BW_P # for adjusted vaccine coverage
N <- log(data$new_cases_month) # Log transformation of the number of new cases per month 
C <- data$country_index # for country interaction with month random effect
M <- data$month # The month index of the dataset 

df <- data.frame(Y1, P, V, N, C, M)


# Stepwise methods were used to screen growth rate model variables using AIC criterion---lineage diversity and adjusted vaccine coverage
mylist <- list( FW_cb1,FW_cb2,FW_cb3,FW_cb4,FW_cb5)


               
names(mylist)<- c( "FW_cb1","FW_cb2","FW_cb3","FW_cb4","FW_cb5")

                 
                 
gam_list <- lapply(mylist, function(x){                                                              
  mgcv::gam(Y1 ~ + x +travel_cb+npi_cb+s(C, bs = "re")+N+P:V, data = df,   method = "REML") 
})   

result <- sapply(X = gam_list , FUN = AIC)
result <- as.data.frame(result) 
result
min(result[,1])
rownames(result)[which.min(result$result)]





model <-  gam(Y1 ~ npi_cb+FW_cb3+travel_cb+N+s(C, bs = "re")+P:V
                 ,method = "REML",data=df)

#pdf("report/plot/report_growth_rate_fig5.pdf",height=6,width=10)
par(oma = c(0, 0, 1, 0),cex.axis=1,cex.lab=1,family = "ArialN")
plot(model,family = "ArialN",cex.lab = 1,cex.main=1,cex.axis=1)
#dev.off()

#pdf("report/plot/report_growth_rate_fig6.pdf",height=5,width=10)

model.res <- residuals(model,type="deviance")
par(mfrow = c(1,2), mar = c(5,4,1,2)+0.2)
hist(model.res,
     xlab = "Residuals", main = "",
     nclass = 20, prob = T,xlim=c(-2.2,4.4))
model.res =
  model.res[which(model.res < abs(min(model.res)))]
normFitDens = dnorm(-22:44/10, mean =  8.497202e-15, sd = 0.8218212)
lines(-22:44/10, normFitDens)
qqnorm(model.res, main ="", col = "#00000077", pch = 16,ylim=c(-2,2))
qqline(model.res, col="red")
#dev.off()


