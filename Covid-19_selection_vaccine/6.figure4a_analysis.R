# # install.packages("lubridate")
rm(list = ls())
setwd("D:/YLY/ccm_dnds_death")
getwd()
# install.packages("devtools")
# devtools::install_github("ha0ye/rEDM")
library(lubridate);library(rEDM);library(tidyr)
library(ggpubr);library(reshape2);library(ggplot2);library(gplots)
library(RColorBrewer);library(readxl);library(lubridate);library(gtable)
library(reshape);library(scales);library(showtext);library(dplyr);library(Kendall) 

sessionInfo()

mytheme <- theme_bw()+
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
    legend.background = element_rect(fill=NA, size=0,color=NA),
    legend.text=element_text(size=16),
    legend.title=element_text(face="bold",size=16),
    axis.line.y=element_line(linetype=1,color='black',size=0.5),
    axis.line.x=element_line(linetype=1,color='black',size=0.5),
    axis.ticks = element_line(linetype=2,color='black'),
    panel.grid=element_line(linetype=2,color='grey'),
    plot.title = element_text(hjust = 0.5,size=18, face = "bold"),
    axis.title.y.left = element_text(size = 16,color="black",vjust=2),
    axis.title.y.right =element_text(size = 16,color="black",vjust=0,angle=90),
    axis.title.x = element_text(size = 16, color="black",vjust = 0),
    axis.text.y.left = element_text(size = 16,color="black",vjust = 0.5,angle = 0),
    axis.text.y.right = element_text(size = 16,color="black",vjust = 0.5,angle = 0),
    axis.text.x = element_text(size = 16, color="black",vjust = 0.5, angle = 0),
    axis.ticks.length=unit(0.15,"cm"),
    axis.ticks.y=element_line(color="black",size=.5),
    axis.ticks.x=element_line(color="black",size=.5),
    plot.margin=unit(c(1,1,1,1),'lines'),
    panel.border = element_rect(colour = "white", fill=NA, size=1)
  )
dir.create("02_output/ccm_2179_death_E2_4")
outpath <-  "02_output/ccm_2179_death_E2_4/"

dir.create(paste0(outpath,"All_result_p_2179"))
outpath1 <-  paste0(outpath,"All_result_p_2179","/")


dir.create(paste0(outpath,"All_parameter_selection"))
outpath4 <-  paste0(outpath,"All_parameter_selection","/")

seed_value <- 2179

dat<-read.csv("01_data/subsample_dnds_vaccine_diversity_low_death.csv")
dat1 <- dat
#view(data)
dat1$date <- parse_date_time(dat1$month, "ym")
dat1$date <-as.Date(as.POSIXct(dat1$date,tz="Hongkong"))
dat1$month1 <- month(dat1$date)
dat1$year <- year(dat1$date)
dat1$NPI <- dat1$NPI_NV

colnames(dat1)[19]<-"death"
colnames(dat1)[20]<-"death_rate"
colnames(dat1)[4] <- "dnds"
colnames(dat1)[6] <- "natural_immunity"
colnames(dat1)[9] <- "vaccine"

##Remove January and February 2020, as the submitted sequences of most countries are unavailable.
y1 <- which(dat1$month=="2020-01"|dat1$month=="2020-02") 
dat1 <-dat1[-y1,]

##Remove November 2021 to September 2022 as the Omicron variant emerged and circulated during this period.
y2 <- which(dat1$month=="2021-11"|dat1$month=="2021-12"|dat1$month=="2022-01"|dat1$month=="2022-02"|
              dat1$month=="2022-03"|dat1$month=="2022-04"|dat1$month=="2022-05"|dat1$month=="2022-06"|
              dat1$month=="2022-07"|dat1$month=="2022-08"|dat1$month=="2022-09") 
dat1 <-dat1[-y2,]

#y1 <- which(dat1$month=="2019-12")
#y2 <- which(dat1$month=="2020-01")
#y3 <- which(dat1$month=="2020-02")
#y4 <- which(dat1$month=="2022-09")
#dat1 <-dat1[-c(y1,y2,y3,y4),]

##Remove countries with missing data for other indexes such as vaccines or PHSM (Public Health and Social Measures).
x <- which(dat1$country=="Puerto Rico"|dat1$country=="Canary Islands"|dat1$country=="French Guiana"|
             dat1$country=="Reunion"|dat1$country=="Curacao"|dat1$country=="North Macedonia"|
             dat1$country=="Slovakia")
dat1 <-dat1[-x,]

dat1$travel2 <- log(dat1$travel_route)
dat1$travel2<-as.numeric(dat1$travel2)

dat1$gene <- trimws(dat1$gene)
gene_type <- unique(dat1$gene)
s=length(gene_type)
print(gene_type)


country_name <- trimws(unique(dat1[,"country"]))
n =length(country_name)
print(country_name)


#for(co1 in c(1:40,42:81,83:84,86:length(country_name)))

for(co1 in c(1:length(country_name)))
#for (co1 in 84)
  {
#dat2= subset(dat1,dat1$gene== "S")
print(country_name[co1])
dat2= subset(dat1,dat1$country== country_name[co1])

ccm_p_matrix_gene <- matrix(NA, nrow = s, ncol = 6)
colnames(ccm_p_matrix_gene) <- c("dnds_xmap_death","death_xmap_dnds","dnds_xmap_death_rate",
                                 "death_rate_xmap_dnds",
                                 "gene","country")



E_parameter_matrix <- matrix(NA, nrow = s, ncol = 17)
E_parameter_matrix <- as.data.frame(E_parameter_matrix)
colnames(E_parameter_matrix) <- c("E_dnds","E_death","E_death_rate", 
                                  "theta_dnds","theta_death","theta_death_rate",
                                  "theta_best_rho_dnds","theta_best_rho_death","theta_best_rho_death_rate",
                                  "theta_0_rho_dnds","theta_0_rho_death","theta_0_rho_death_rate",
                                  "theta_diff_dnds","theta_diff_death","theta_diff_death_rate",
                                 "gene","country")

for(g in 1:s)
#for(g in 8)
  {
  print(gene_type[g])
  dat3= subset(dat2,dat2$gene== gene_type[g])
  
 
  #dat2= subset(dat1,dat1$gene== "S") 
  
  dir.create(paste0(outpath,country_name[co1]))
  outpath_c1 <-  paste0(outpath,country_name[co1],"/")
  
  dir.create(paste0(outpath,country_name[co1],"/","data_visualization"))
  outpath_c2 <-  paste0(outpath,country_name[co1],"/","data_visualization","/")
  
  dir.create(paste0(outpath,country_name[co1],"/","parameter_selection"))
  outpath_c3 <-  paste0(outpath,country_name[co1],"/","parameter_selection","/")
  
  dir.create(paste0(outpath,country_name[co1],"/","ccm_result"))
  outpath_c4 <-  paste0(outpath,country_name[co1],"/","ccm_result","/")
  
  dir.create(paste0(outpath,country_name[co1],"/","res_p_value"))
  outpath_c5 <-  paste0(outpath,country_name[co1],"/","res_p_value","/")
  
  
data <- dat3
#data <- subset(dat1,dat1$gene== "S")

data <- data %>%
  group_by(month) %>%
  mutate(month1 = cur_group_id())

c1='#34ace0'
c2='#e74c3c'
c3='#f1c40f'
c4='#706fd3'
c5='#009432'


#c5="#ED0000B2"
dat_plot <- data 
A=max(dat_plot$dnds,na.rm = T)/max(dat_plot$death,na.rm = T)

pd1<-ggplot(dat_plot,aes(month,group=1)) +
  geom_line(aes(x=month,y=dnds,color="dN/dS ratio"),size=1.5)+
  geom_line(aes(x=month,y=dat_plot$death*max(dat_plot$dnds,na.rm = T)*A,color="Death"),size=1.5)+
  #facet_wrap(~country,scales = "free",ncol=8)+
  scale_y_continuous(sec.axis = sec_axis(~./(max(dat_plot$dnds,na.rm = T)*A),name="Deaths"))+
  labs(x="", y= "dN/dS ratio",color="" )+
  geom_hline(yintercept = 1,lty='dashed',colour ="black",size=1)+
  scale_x_discrete(breaks=c("2020-03","2020-09","2021-03","2021-09"),
                   labels=c("Mar'20","Sep'20","Mar'21","Sep'21"))+
  scale_colour_manual(values=c(c1,c2,c3,c4))+
  ggtitle(paste0(country_name[co1]," (",gene_type[g],")"))+
  mytheme+
  theme(strip.text = element_text(face="bold", family = "serif",size=16),
        strip.background = element_rect(fill="white", colour="white"))+
  theme(panel.background=element_rect(fill = "white"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  guides(color=guide_legend(nrow=1,byrow = T))

#ggsave(pd1,file=paste0(outpath_c2,"diversity_dnds_",country_name[co1],"_",gene_type[g],".pdf"),
#       width = 8,height =5) 


B=max(dat_plot$dnds,na.rm = T)/max(dat_plot$death_rate,na.rm = T)

pd2<-ggplot(dat_plot,aes(month,group=1)) +
  geom_line(aes(x=month,y=dnds,color="dN/dS ratio"),size=1.5)+
  geom_line(aes(x=month,y=dat_plot$death_rate*max(dat_plot$dnds,na.rm = T)*B,color="Deaths per million"),size=1.5)+
  #facet_wrap(~country,scales = "free",ncol=8)+
  scale_y_continuous(sec.axis = sec_axis(~./(max(dat_plot$dnds,na.rm = T)*B),name="Deaths per million"))+
  labs(x="", y= "dN/dS ratio",color="" )+
  geom_hline(yintercept = 1,lty='dashed',colour ="black",size=1)+
  scale_x_discrete(breaks=c("2020-03","2020-09","2021-03","2021-09"),
                   labels=c("Mar'20","Sep'20","Mar'21","Sep'21"))+
  scale_colour_manual(values=c(c4,c2))+
  ggtitle(paste0(country_name[co1]," (",gene_type[g],")"))+
  mytheme+
  theme(strip.text = element_text(face="bold", family = "serif",size=16),
        strip.background = element_rect(fill="white", colour="white"))+
  theme(panel.background=element_rect(fill = "white"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  guides(color=guide_legend(nrow=1,byrow = T))

#ggsave(pd1,file=paste0(outpath_c2,"diversity_dnds_",country_name[co1],"_",gene_type[g],".pdf"),
#       width = 8,height =5) 




pd3<-ggplot(dat_plot,aes(month,group=1)) +
  geom_line(aes(x=month,y=vaccine,color="Adjusted vaccine coveage"),size=1.5)+
  geom_line(aes(x=month,y=natural_immunity,color="Natural immunity"),size=1.5)+
  #facet_wrap(~country,scales = "free",ncol=8)+
  labs(x="", y= "Coverage",color="" )+
  geom_hline(yintercept = 1,lty='dashed',colour ="black",size=1)+
  scale_x_discrete(breaks=c("2020-09","2021-03","2021-09"),
                   labels=c("Sep'20","Mar'21","Sep'21"))+
  scale_colour_manual(values=c(c3,c5))+
  ggtitle(paste0(country_name[co1]))+
  mytheme+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  guides(color=guide_legend(nrow=1,byrow = T))

#ggsave(pd2,file=paste0(outpath_c2,"vaccine_NI_",country_name[co1],".pdf"),
#       width = 8,height =5) 



normalize <- function(x,na.rm=TRUE, ...) {
  (x - mean(x,na.rm=TRUE, ...))/sd(x,na.rm=TRUE, ...)
}
# separate time column from data
vars <- c("dnds","death","death_rate")
composite_ts <- data[, vars]
#colnames(composite_ts) <- (c("Shannon", "Growth_Rate", "NPI", "vaccine","travel2"))
# normalize each time series within a country
data_by_country <- split(composite_ts, data$country)
normalized_data <- lapply(data_by_country, function(df) sapply(df, normalize))
#w=(1:30)

#data$month1=rep(w,63)
composite_ts <- cbind(month1 = data$month1, data.frame(do.call(rbind, normalized_data)))


segments_end <- cumsum(sapply(data_by_country, NROW))
segments_begin <- c(1, segments_end[-length(segments_end)] + 1)
segments <- cbind(segments_begin, segments_end)


# Choose random segments for prediction
set.seed(seed_value)

#rndlib <- sample(1:NROW(segments), floor(NROW(segments) * 0.75))
#composite_lib <- segments[rndlib, ]
#composite_pred <- segments[-rndlib, ]
composite_lib <- c(1, NROW(composite_ts))
composite_pred <- c(1, NROW(composite_ts))
#Quantifying predictability 

simplex_out <- lapply(vars, function(var) {
  simplex(composite_ts[, c("month1", var)], silent = TRUE,E = 2:4, lib = composite_lib, pred = composite_pred)
})

names(simplex_out) <- vars

par(mfrow = c(2, 2))
for (var in names(simplex_out)) {
  plot(simplex_out[[var]]$E, simplex_out[[var]]$rho, type = "l", xlab = "Embedding Dimension (E)", 
       ylab = "Forecast Skill (rho)", main = var)
}


best_E <- sapply(simplex_out, function(df) {
  df$E[which.max(df$rho-df$mae)]
})
best_E



####################
pdf(file=paste0(outpath_c3,"fig_ccm_E_selection","_",gene_type[g],"_",country_name[co1],".pdf"), width =12, height = 8)
par(mfrow = c(2, 3))
par(mar = c(3.3,5,5,4.5), mgp = c(2,0.5,0), oma = c(0.5,0,0,0))
for (var in names(simplex_out)) {
  plot(simplex_out[[var]]$E, simplex_out[[var]]$rho, plot.main=NULL, type="l",  lwd=3, col="#235389",
       axes = FALSE,
       xlab = "", ylab = "")
  box(col="black")
  axis(1, col="black", col.ticks="black", col.axis="black", cex.axis=1, tck = 0.02)
  mtext("Embedding dimension (E)", side=1, line=2, col="black", cex=1.2)
  
  axis(2, col="#235389", col.ticks="#235389", col.axis="#235389", cex.axis=1, tck = 0.02)
  mtext(expression(paste("Correlation coefficient (",rho,")")), side=2, line=2, col="#235389", cex=1.2)
  par(new = T)
  plot(simplex_out[[var]]$E, simplex_out[[var]]$mae, type="l", lwd=3, col="#D59730", axes = F, xlab = NA, ylab = NA)
  axis(side = 4, col="#D59730", col.ticks="#D59730", col.axis="#D59730", cex.axis=1, tck = 0.02)
  mtext(side = 4, line = 2, "Mean Absolute Error", col = "#D59730", cex=1.2)
  
  abline(v = best_E[var], col="#922326", lwd=3, lty=2)
  title(var,cex.main=2)
}



#quantifying nonlinearity
smap_out <- lapply(vars, function(var) {
  s_map(composite_ts[, c("month1", var)], E = best_E[var], lib = composite_lib, 
        pred = composite_pred)
})
names(smap_out) <- names(simplex_out)


#par(mfrow = c(2, 2))
for (var in names(smap_out)) {
  plot(smap_out[[var]]$theta, smap_out[[var]]$rho, plot.main=NULL, type = "l",lwd=2,
       axes = FALSE,
       xlab = "", ylab = "")
  box(col="black")
  axis(1, col="black", col.ticks="black", col.axis="black", cex.axis=1, tck = 0.02)
  mtext("Nonlinearity (theta)", side=1, line=2, col="black", cex=1.2)
  
  axis(2, col="black", col.ticks="black", col.axis="black", cex.axis=1, tck = 0.02)
  mtext("Forecast Skill (rho)", side=2, line=2, col="black", cex=1.2)
  title(var,cex.main=2)
}

dev.off()

#####################
E_parameter_matrix[g,"E_dnds"] <- best_E["dnds"]
E_parameter_matrix[g,"E_death"] <- best_E["death"]
E_parameter_matrix[g,"E_death_rate"] <- best_E["death_rate"]

E_parameter_matrix[g,"theta_0_rho_dnds"] <- smap_out[["dnds"]]$rho[1]
E_parameter_matrix[g,"theta_0_rho_death"] <- smap_out[["death"]]$rho[1]
E_parameter_matrix[g,"theta_0_rho_death_rate"] <- smap_out[["death_rate"]]$rho[1]

E_parameter_matrix[g,"theta_best_rho_dnds"] <- max(smap_out[["dnds"]]$rho)
E_parameter_matrix[g,"theta_best_rho_death"] <- max(smap_out[["death"]]$rho)
E_parameter_matrix[g,"theta_best_rho_death_rate"] <- max(smap_out[["death_rate"]]$rho)

E_parameter_matrix[g,"theta_dnds"] <- smap_out[["dnds"]]$theta[which.max(smap_out[["dnds"]]$rho)]
E_parameter_matrix[g,"theta_death"] <-smap_out[["death"]]$theta[which.max(smap_out[["death"]]$rho)]
E_parameter_matrix[g,"theta_death_rate"] <- smap_out[["death_rate"]]$theta[which.max(smap_out[["death_rate"]]$rho)]

E_parameter_matrix[g,"theta_diff_dnds"] <- max(smap_out[["dnds"]]$rho) - smap_out[["dnds"]]$rho[1]
E_parameter_matrix[g,"theta_diff_death"] <- max(smap_out[["death"]]$rho) - smap_out[["death"]]$rho[1]
E_parameter_matrix[g,"theta_diff_death_rate"] <- max(smap_out[["death_rate"]]$rho) - smap_out[["death_rate"]]$rho[1]

E_parameter_matrix[g,"gene"] <-  gene_type[g]
E_parameter_matrix[g,"country"] <-  country_name[co1]

write.table(E_parameter_matrix,file =paste0(outpath4,"parameter_E_theta","_",country_name[co1],".csv"),row.names=FALSE, sep=",")

#####################
#Convergent Cross Mapping 

lib_sizes <- c(seq(2, 20, by = 2))


png(file=paste0(outpath_c4,"/fig_ccm_death_dnds","_",gene_type[g],"_",country_name[co1],".png"), width = 950, height = 500,res=150)

par(mfrow=c(1,2),mar = c(4, 4, 4, 1),oma = c(0, 2, 0, 0))

####
death_xmap_dnds <- ccm(composite_ts, lib = segments, pred = segments, lib_column = "death", 
                   target_column = "dnds", replace = FALSE,E = best_E["death"], lib_sizes = lib_sizes, 
                   silent = TRUE)
dnds_xmap_death <- ccm(composite_ts, lib = segments, pred = segments, lib_column = "dnds", 
                   target_column = "death", replace = FALSE, E = best_E["dnds"], lib_sizes = lib_sizes, 
                   silent = TRUE)

death_xmap_dnds1 <- death_xmap_dnds
death_xmap_dnds1$gene <- gene_type[g]
death_xmap_dnds1$country <- country_name[co1]
death_xmap_dnds1$seed <- seed_value

dnds_xmap_death1 <- dnds_xmap_death
dnds_xmap_death1$gene <- gene_type[g]
dnds_xmap_death1$country <- country_name[co1]
dnds_xmap_death1$seed <- seed_value

write.table(death_xmap_dnds1,file =paste0(outpath_c1,"death_xmap_dnds","_",gene_type[g],"_",country_name[co1],".csv"),row.names=FALSE, sep=",")
write.table(dnds_xmap_death1,file =paste0(outpath_c1,"dnds_xmap_death","_",gene_type[g],"_",country_name[co1],".csv"),row.names=FALSE, sep=",")
death_xmap_dnds_means <- ccm_means(death_xmap_dnds)
dnds_xmap_death_means <- ccm_means(dnds_xmap_death)
#par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0)) #
plot(death_xmap_dnds_means$lib_size, pmax(0, death_xmap_dnds_means$rho), type = "l", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", col = "red", ylim = c(0, 
                                                                                  1), lwd = 2)
lines(dnds_xmap_death_means$lib_size, pmax(0, dnds_xmap_death_means$rho), col = "blue", 
      lwd = 2)
legend(x = "topleft", col = c("red", "blue"), lwd = 2, legend = c("Deaths xmap dNdS", 
                                                                  "dNdS xmap deaths"), inset = 0.02, bty = "n", cex = 0.8)

abline(h = 0, lty = 3)



####
dea_rate_xmap_dnds <- ccm(composite_ts, lib = segments, pred = segments, lib_column = "death_rate", 
                      target_column = "dnds",replace = FALSE, E = best_E["death_rate"], lib_sizes = lib_sizes, 
                      silent = TRUE)
dnds_xmap_dea_rate <- ccm(composite_ts, lib = segments, pred = segments, lib_column = "dnds", 
                          target_column = "death_rate",replace = FALSE, E = best_E["dnds"], lib_sizes = lib_sizes, 
                          silent = TRUE)

dea_rate_xmap_dnds1 <- dea_rate_xmap_dnds
dea_rate_xmap_dnds1$gene <- gene_type[g]
dea_rate_xmap_dnds1$country <- country_name[co1]
dea_rate_xmap_dnds1$seed <- seed_value


dnds_xmap_dea_rate1 <- dnds_xmap_dea_rate
dnds_xmap_dea_rate1$gene <- gene_type[g]
dnds_xmap_dea_rate1$country <- country_name[co1]
dnds_xmap_dea_rate1$seed <- seed_value


write.table(dea_rate_xmap_dnds1,file = paste0(outpath_c1,"dea_rate_xmap_dnds","_",gene_type[g],"_",country_name[co1],".csv"),row.names=FALSE, sep=",")
write.table(dnds_xmap_dea_rate1,file = paste0(outpath_c1,"dnds_xmap_dea_rate","_",gene_type[g],"_",country_name[co1],".csv"),row.names=FALSE, sep=",")

dea_rate_xmap_dnds_means <- ccm_means(dea_rate_xmap_dnds)
dnds_xmap_dea_rate_means <- ccm_means(dnds_xmap_dea_rate)

plot(dea_rate_xmap_dnds_means$lib_size, pmax(0, dea_rate_xmap_dnds_means$rho), type = "l", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", col = "red", ylim = c(0, 
                                                                                  1), lwd = 2)
lines(dnds_xmap_dea_rate_means$lib_size, pmax(0, dnds_xmap_dea_rate_means$rho), col = "blue", 
      lwd = 2)
legend(x = "topleft", col = c("red", "blue"), lwd = 2, legend = c("Death rate xmap dNdS", 
                                                                  "dNdS xmap death rate"), inset = 0.02, bty = "n", cex = 0.8)
abline(h = 0, lty = 3)


dev.off()


######Null model######
n <- NROW(composite_ts)
ccm_rho_matrix <- matrix(NA, nrow = length(vars), ncol = length(vars), dimnames = list(vars,vars))

ccm_all_matrix <- matrix(NA,  nrow = length(vars)*(length(vars)-1), ncol=11)
ccm_all_matrix <- as.data.frame(ccm_all_matrix)
rownames(ccm_all_matrix) <- c("death:dnds", "dnds:death","death_rate:dnds","dnds:death_rate",
                              "death:death_rate", "death_rate:death")


for (ccm_from in vars) {
  for (ccm_to in vars[vars != ccm_from]) {
    out_temp <- ccm(composite_ts, E = best_E[ccm_from], lib = segments, pred = segments,
                    lib_column = ccm_from, target_column = ccm_to, 
                    lib_sizes = n, replace = FALSE, silent = TRUE)
    ccm_rho_matrix[ccm_from, ccm_to] <- out_temp$rho
    ccm_all_matrix[paste0(ccm_from,":",ccm_to),] <- out_temp
  }
}


colnames(ccm_all_matrix) <- colnames(out_temp)

num_surr <- 1000


#Generate surrogate time series with the random shuffled based on composite time series data
surr_dnds <- make_surrogate_data(composite_ts$dnds, method = "random_shuffle", 
                                 num_surr = num_surr)

surr_death <- make_surrogate_data(composite_ts$death, method = "random_shuffle", 
                                     num_surr = num_surr)

surr_death_rate <- make_surrogate_data(composite_ts$death_rate, method = "random_shuffle", 
                                     num_surr = num_surr)

ccm_rho_surr <- data.frame(dnds_xmap_death = numeric(num_surr), death_xmap_dnds = numeric(num_surr),
                           dnds_xmap_death_rate= numeric(num_surr),death_rate_xmap_dnds= numeric(num_surr)
)


#plot surrogate time series


c3='#7f8c8d'
c4='#ff793f'


dfR <- cbind(composite_ts$month1,composite_ts$death,surr_death[,1])
dfR <- as.data.frame(dfR)
colnames(dfR) <- c("month","Observed","Random")
dfR$Observed <- as.numeric(dfR$Observed)
dfR$Random <- as.numeric(dfR$Random)
p1_d1 <- ggplot(dfR,aes(month,group=1)) +
  geom_line(aes(x=month,y=Observed,color="Observed"),size=1.5)+
  geom_line(aes(x=month,y=Random,color="random_shuffle"),size=1.5)+
  geom_point(aes(x=month,y=Observed,color="Observed"),size=5)+
  geom_point(aes(x=month,y=Random,color="random_shuffle"),size=5)+
  scale_x_continuous(breaks=c(1,7,13,19),
                   labels=c("Mar'20","Sep'20","Mar'21","Sep'21"))+
  scale_colour_manual(values=c(c3,c4))+
  labs(color="" )+
  ggtitle("Deaths")+
  mytheme+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  guides(color=guide_legend(nrow=1,byrow = T))


dfR <- cbind(composite_ts$month,composite_ts$death_rate,surr_death_rate[,1])
dfR <- as.data.frame(dfR)
colnames(dfR) <- c("month","Observed","Random")
dfR$Observed <- as.numeric(dfR$Observed)
dfR$Random <- as.numeric(dfR$Random)
p1_d2 <- ggplot(dfR,aes(month,group=1)) +
  geom_line(aes(x=month,y=Observed,color="Observed"),size=1.5)+
  geom_line(aes(x=month,y=Random,color="random_shuffle"),size=1.5)+
  geom_point(aes(x=month,y=Observed,color="Observed"),size=5)+
  geom_point(aes(x=month,y=Random,color="random_shuffle"),size=5)+
  scale_x_continuous(breaks=c(1,7,13,19),
                     labels=c("Mar'20","Sep'20","Mar'21","Sep'21"))+
  scale_colour_manual(values=c(c3,c4))+
  labs(color="" )+
  ggtitle("Richness")+
  mytheme+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  guides(color=guide_legend(nrow=1,byrow = T))



dfR <- cbind(composite_ts$month,composite_ts$dnds,surr_dnds[,1])
dfR <- as.data.frame(dfR)
colnames(dfR) <- c("month","Observed","Random")
dfR$Observed <- as.numeric(dfR$Observed)
dfR$Random <- as.numeric(dfR$Random)
p1_d3 <- ggplot(dfR,aes(month,group=1)) +
  geom_line(aes(x=month,y=Observed,color="Observed"),size=1.5)+
  geom_line(aes(x=month,y=Random,color="random_shuffle"),size=1.5)+
  geom_point(aes(x=month,y=Observed,color="Observed"),size=5)+
  geom_point(aes(x=month,y=Random,color="random_shuffle"),size=5)+
  scale_x_continuous(breaks=c(1,7,13,19),
                     labels=c("Mar'20","Sep'20","Mar'21","Sep'21"))+
  scale_colour_manual(values=c(c3,c4))+
  labs(color="" )+
  ggtitle("dN/dS ratio")+
  mytheme+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  guides(color=guide_legend(nrow=1,byrow = T))


pd3<-ggarrange(pd1,pd2,pd3,p1_d1,p1_d2,p1_d3,ncol =3, nrow = 2,
               labels = c("","","","",""),align = "hv",
               font.label = list(size = 16, face = "bold"))

ggsave(pd3,file=paste0(outpath_c2,"data_surr_normalization_",country_name[co1],"_",gene_type[g],".pdf"),
       width = 24,height =10) 

##############################################################
for (i in 1:num_surr) {
  ccm_rho_surr$dnds_xmap_death[i] <- ccm(cbind(composite_ts$dnds, surr_death[, i]),lib_column = 1,
                                      E = best_E["dnds"],  target_column = 2, lib_sizes = NROW(composite_ts), 
                                      replace = FALSE, silent = TRUE)$rho
  
  ccm_rho_surr$death_xmap_dnds[i] <- ccm(cbind(composite_ts$death, surr_dnds[,i]), lib_column = 1, 
                                      E = best_E["death"], target_column = 2, lib_sizes = NROW(composite_ts), 
                                      replace = FALSE, silent = TRUE)$rho
  
  ccm_rho_surr$dnds_xmap_death_rate[i] <- ccm(cbind(composite_ts$dnds, surr_death_rate[, i]),lib_column = 1,
                                        E = best_E["dnds"],  target_column = 2, lib_sizes = NROW(composite_ts), 
                                        replace = FALSE, silent = TRUE)$rho
  
  ccm_rho_surr$death_rate_xmap_dnds[i] <- ccm(cbind(composite_ts$death_rate, surr_death_rate[,i]), lib_column = 1, 
                                        E = best_E["death_rate"], target_column = 2, lib_sizes = NROW(composite_ts), 
                                        replace = FALSE, silent = TRUE)$rho
  
}



p_dnds_xmap_death <- (sum(ccm_rho_matrix['dnds', 'death'] < ccm_rho_surr$dnds_xmap_death) + 1) /
  (length(ccm_rho_surr$dnds_xmap_death) + 1)


p_death_xmap_dnds <- (sum(ccm_rho_matrix['death', 'dnds'] < ccm_rho_surr$death_xmap_dnds) + 1) /
  (length(ccm_rho_surr$death_xmap_dnds) + 1)


p_dnds_xmap_death_rate <- (sum(ccm_rho_matrix['dnds', 'death_rate'] < ccm_rho_surr$dnds_xmap_death_rate) + 1) /
  (length(ccm_rho_surr$dnds_xmap_death_rate) + 1)


p_death_rate_xmap_dnds <- (sum(ccm_rho_matrix['death_rate', 'dnds'] < ccm_rho_surr$death_rate_xmap_dnds) + 1) /
  (length(ccm_rho_surr$death_rate_xmap_dnds) + 1)


#ccm_p_matrix_gene <- matrix(NA, nrow = 1, ncol = 7)
#colnames(ccm_p_matrix_gene) <- c("dnds_xmap_SI","SI_xmap_dnds","dnds_xmap_rich",
#                                 "rich_xmap_dnds","dnds_xmap_even","even_xmap_dnds",
#                                 "gene")

ccm_p_matrix_gene[g,"dnds_xmap_death"] <- p_dnds_xmap_death
ccm_p_matrix_gene[g,"death_xmap_dnds"] <- p_death_xmap_dnds
ccm_p_matrix_gene[g,"dnds_xmap_death_rate"] <- p_dnds_xmap_death_rate
ccm_p_matrix_gene[g,"death_rate_xmap_dnds"] <- p_death_rate_xmap_dnds

ccm_p_matrix_gene[g,"gene"] <- gene_type[g]
ccm_p_matrix_gene[g,"country"] <- country_name[co1]


write.csv(ccm_p_matrix_gene,file=paste0(outpath_c4,"dnds_death_p_value","_",country_name[co1],".csv",sep=""),row.names =F)


df_ccm_rho_surr <- gather(ccm_rho_surr,Type,rho_surr,1:4)
df_ccm_rho_surr$rho_actual <- NULL
df_ccm_rho_surr$sign <- NULL
df_ccm_rho_surr$sign1 <- NULL
df_ccm_rho_surr$kend_tau <- NULL
df_ccm_rho_surr$kend_p <- NULL

df_ccm_rho_surr[which(df_ccm_rho_surr$Type=="dnds_xmap_death"),"rho_actual"] <- ccm_rho_matrix['dnds', 'death']
df_ccm_rho_surr[which(df_ccm_rho_surr$Type=="death_xmap_dnds"),"rho_actual"] <- ccm_rho_matrix['death', 'dnds']
df_ccm_rho_surr[which(df_ccm_rho_surr$Type=="dnds_xmap_death_rate"),"rho_actual"] <- ccm_rho_matrix['dnds', 'death_rate']
df_ccm_rho_surr[which(df_ccm_rho_surr$Type=="death_rate_xmap_dnds"),"rho_actual"] <- ccm_rho_matrix['death_rate', 'dnds']

df_ccm_rho_surr[which(df_ccm_rho_surr$Type=="dnds_xmap_death"),"kend_tau"] <- MannKendall(dnds_xmap_death$rho)[1]
df_ccm_rho_surr[which(df_ccm_rho_surr$Type=="death_xmap_dnds"),"kend_tau"] <- MannKendall(death_xmap_dnds$rho)[1]
df_ccm_rho_surr[which(df_ccm_rho_surr$Type=="dnds_xmap_death_rate"),"kend_tau"] <- MannKendall(dnds_xmap_dea_rate$rho)[1]
df_ccm_rho_surr[which(df_ccm_rho_surr$Type=="death_rate_xmap_dnds"),"kend_tau"] <- MannKendall(dea_rate_xmap_dnds$rho)[1]

df_ccm_rho_surr[which(df_ccm_rho_surr$Type=="dnds_xmap_death"),"kend_p"] <- MannKendall(dnds_xmap_death$rho)[2]
df_ccm_rho_surr[which(df_ccm_rho_surr$Type=="death_xmap_dnds"),"kend_p"] <- MannKendall(death_xmap_dnds$rho)[2]
df_ccm_rho_surr[which(df_ccm_rho_surr$Type=="dnds_xmap_death_rate"),"kend_p"] <- MannKendall(dnds_xmap_dea_rate$rho)[2]
df_ccm_rho_surr[which(df_ccm_rho_surr$Type=="death_xmap_dnds"),"kend_p"] <- MannKendall(dea_rate_xmap_dnds$rho)[2]



df_ccm_rho_surr[which(df_ccm_rho_surr$Type=="dnds_xmap_death"),"sign"] <- p_dnds_xmap_death
df_ccm_rho_surr[which(df_ccm_rho_surr$Type=="death_xmap_dnds"),"sign"] <- p_death_xmap_dnds
df_ccm_rho_surr[which(df_ccm_rho_surr$Type=="dnds_xmap_death_rate"),"sign"] <- p_dnds_xmap_death_rate
df_ccm_rho_surr[which(df_ccm_rho_surr$Type=="death_rate_xmap_dnds"),"sign"] <- p_death_rate_xmap_dnds


df_ccm_rho_surr[which(df_ccm_rho_surr$sign<=0.05),"sign1"] <- 1
df_ccm_rho_surr[which(df_ccm_rho_surr$sign>0.05),"sign1"] <- 0
df_ccm_rho_surr$sign1 <- as.factor(df_ccm_rho_surr$sign1)

df_ccm_rho_surr$gene <- gene_type[g]
df_ccm_rho_surr$country <- country_name[co1]
df_ccm_rho_surr$seed <- seed_value

df_ccm_rho_surr$Type <-factor(df_ccm_rho_surr$Type,levels= c("dnds_xmap_death","death_xmap_dnds",
                                                                 "dnds_xmap_death_rate","death_rate_xmap_dnds",
                                                                 ordered=TRUE))  

write.csv(df_ccm_rho_surr,file=paste0(outpath_c5,"result_death_dnds_surr","_",country_name[co1],"_",gene_type[g],".csv",sep=""),row.names =F)
write.csv(df_ccm_rho_surr,file=paste0(outpath1,"result_p_death_dnds_surr","_",country_name[co1],"_",gene_type[g],".csv",sep=""),row.names =F)

custom_colors <- c("#e41a1c","#f4a582","#377eb8","#92c5de","#4daf4a","#99d594")


pr1 <- ggplot(df_ccm_rho_surr, aes(x=Type, y=rho_surr,fill=Type)) + 
geom_boxplot(notch=F, outlier.shape = NA,outlier.size = 1,
             outlier.fill = "black")+  
geom_point(aes(x=Type,y=rho_actual,shape=factor(sign1)),color="black",size=3,stroke = 2)+
geom_point(position=position_jitterdodge(jitter.width=2, dodge.width = 0), 
           pch=21, aes(fill=factor(Type),color=factor(Type)), show.legend = F,size=0.5)+
coord_cartesian(ylim = c(0, 1))+
scale_y_continuous(name="Cross Mapping Skill (p)", breaks = seq(0,1,0.2))+
labs(x=NULL,color="" ,fill="")+
mytheme+
scale_fill_manual(name = '',values = custom_colors) +
scale_color_manual(name = '',values = custom_colors) +
scale_shape_manual(name = '',values = c(1, 19),label=c("p > 0.05","p < 0.05")) +
ggtitle(paste0(country_name[co1]," (",gene_type[g],")"))+
#guides(shape = FALSE,color = guide_legend(override.aes = list(shape = NA, size = 3)))   
guides(color = guide_legend(override.aes = list(shape = NA, size = 3)))   

pr1 

ggsave(file=paste0(outpath_c5,"fig_ccm_surr_",gene_type[g],"_",country_name[co1],".pdf"), plot=pr1, width=15, height=5)





}


}



