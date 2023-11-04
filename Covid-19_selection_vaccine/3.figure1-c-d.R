# # install.packages("lubridate")
rm(list = ls())
setwd("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine")
dat<-read.csv("01_data/subsample_dnds_vaccine_diversity_low_death.csv")
dat1 <- dat

source("load_package.R")
source("plot_theme.R")
source("load_data_clean.R")

dir.create("03_output/figure_1")
outpath <-  "03_output/figure_1/"

#Filter out the countries with a full-dose vaccine coverage rate exceeding 75% as of September 2022.
full_vero_2110 <- dat1[which(dat1$gene=="S1"&dat1$month=="2021-10"),]
country_high <- full_vero_2110[which(full_vero_2110$F_Vero>=75),"country"]

#Filter out countries with a full-dose COVID-19 vaccine coverage rate below 50% as of September 2022.
country_low <- full_vero_2110[which(full_vero_2110$F_Vero<=25),"country"]

dat2= subset(dat1,dat1$gene== "S1")
dat2$Income_group <- NA
dat2$country <- as.character(dat2$country)
dat2[which(dat2$country %in% country_high$country),"Income_group"] <- "High coverage region"
dat2[which(dat2$country %in% country_low$country),"Income_group"] <- "Low coverage region"



dat3 <- dat2[-which(is.na(dat2$Income_group)),]
dat3 <- dat3[,c("country","month","gene","F_Vero","dn_ds_mean","Income_group")]


month_name <- unique(dat3$month)
print(month_name)
m =length(month_name)


income_name <- unique(dat3$Income_group)
print(income_name)

dat_final <- matrix(NA,m*length(income_name),5)
dat_final <- as.data.frame(dat_final)


for (ig in 1:length(income_name))
{
print(income_name[ig])
D4= subset(dat3,dat3$Income_group== income_name[ig])
  
var <- c("month","country","dn_ds_mean")
ourdata <-D4[,var] 
colnames(ourdata) <- c("month","country","value")
  
df_ci <- matrix(NA,m,4)
df_ci <- as.data.frame(df_ci)
  
  for (j in 1:m)
  {
    print(month_name[j]);
    result<- matrix(NA,1,3)
    result <- as.data.frame(result) 
    Q2 = subset(ourdata,trimws(ourdata$month)== month_name[j])
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
    colnames(result) <- c("lci","mean","uci")
    result$month <- month_name[j]
    df_ci[j,] <- result
  }
  colnames(df_ci) <- colnames(result)
  df_ci$income_type <- income_name[ig]
  
  dat_final[(c(31*ig-30):(31*ig)),] <- df_ci[,]
  colnames(dat_final) <- colnames(df_ci)
  
}
  

dat_final1 <- dat_final

dat_final1 <- dat_final1[-which(dat_final1$month=="2022-01"|dat_final1$month=="2022-02"|
                                dat_final1$month=="2022-03"|dat_final1$month=="2022-04"|dat_final1$month=="2022-05"|dat_final1$month=="2022-06"|
                                dat_final1$month=="2022-07"|dat_final1$month=="2022-08"|dat_final1$month=="2022-09"), ]

dat_final1$income_type <-factor(dat_final1$income_type,
                                  levels= c('Low coverage region','High coverage region',ordered=TRUE))
p1<-ggplot(dat_final1,aes(x=month,y=mean,group=income_type)) +
  geom_line(aes(x=month,y=mean,color=income_type),size=0.8,alpha=1)+
  geom_ribbon(aes(x=month,ymin = lci, ymax = uci,fill=income_type),alpha=0.55) +
  scale_y_continuous(limits=c(0,8),breaks = seq(0,8,2))+
  labs(x="Month",
       y= "Ratio of nonsynonymous\n to synonymous divergence",color="",fill="")+
  #geom_hline(yintercept = 1,lty='dashed',colour ="black",size=1)+
  scale_x_discrete(breaks=c("2020-03","2020-07","2020-11","2021-03","2021-07","2021-11"),
                   labels=c("Mar\n2020","Jul\n2020","Nov\n2020","Mar\n2021","Jul\n2021","Nov\n2021"))+
  #scale_colour_manual(values= c("#fb9a99","#e31a1c"))+
  scale_colour_manual(values= c("#d4b9da","#5c4591"))+
  scale_fill_manual(values= c( "#d4b9da","#5c4591"))+
  annotate(geom = "text", x = "2021-01", y =2.5,angle=90, label = "Gamma", hjust = 0,size =5)+
  annotate(geom = "text", x = "2021-05", y =2.5,angle=90, label = "Delta", hjust =0,size = 5)+
  annotate(geom = "text", x = "2020-12", y =2.5,angle=90, label = "Alpha/Beta", hjust = 0,size = 5)+
  annotate(geom = "text", x = "2021-11", y =2.5,angle=90, label = "Omicron", hjust =0,size = 5)+
  #annotate(geom = "text", x = "2020-04", y =6.8,angle=0, label = "International travel restrictions", hjust ="left",size = 4.5,family = "ArialN")+
  #annotate(geom = "text", x = "2020-04", y =5.8,angle=0, label = "Public health and social measures", hjust ="left",size = 4.5,family = "ArialN")+
  annotate(geom = "text", x = "2021-03", y =7.5,angle=0, label = "Vaccine rollout", hjust =0,size = 6.5)+
  geom_vline(aes(xintercept="2021-01"),lty='solid',colour ="black",size=0.5)+
  geom_vline(aes(xintercept="2021-05"),lty='solid',colour ="black",size=0.5)+
  geom_vline(aes(xintercept="2020-12"),lty='solid',colour ="black",size=0.5)+
  geom_vline(aes(xintercept="2021-11"),lty='solid',colour ="black",size=0.5)+
  geom_hline(aes(yintercept=1),lty='dashed',colour ="black",size=1)+
  mytheme+
  #theme(axis.title.y.left = element_text(size = 18,color="black",vjust=2,face = "bold", margin = margin(0,0.4,0,0,'cm')),
  #      axis.title.x = element_text(size = 18, color="black",face = "bold",margin = margin(0.4,0,0,0,'cm')))+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+
  theme(legend.position = c(0.2,0.88))+
  guides(color=guide_legend(nrow=2,byrow = T))
  #guides(fill = "none",color="none")
p1 

#ggsave(file=paste0(outpath,"high_75_low_50_visualization_",gene_type[g],".pdf"), plot=p, width=10, height=8)
dat1_S1_202110= subset(dat1,dat1$gene== "S1"&dat1$month=="2021-10")


region <- read.csv("02_sup-data/region_code.csv")
colnames(region) <- c("iso_code","continent","country")
dat1_S1_202110 <- merge(dat1_S1_202110,region,by=c("country"),all.x = T)

selectname <- dat1_S1_202110 %>% filter(country %in% c("India","United Kingdom","Kenya","South Africa",
                                             "Israel","United States","China"))


library(mgcv)


dat4 <- dat1

dat4$total_death_per_million <- dat4$new_deaths_per_million
dat4[which(is.na(dat4$total_death_per_million)),"total_death_per_million"] <- 0
dat4 <- dat4 %>% 
  group_by(country) %>% 
  arrange(country,month) %>% 
  mutate(total_death_per_million = cumsum(total_death_per_million)) %>% 
  ungroup()

a1 <- which(is.na(dat4$total_death_per_million))
dat4[a1,"total_death_per_million"] <- NA

##dat4$total_IHME_case_rate <- dat4$cases_IHME_total/dat4$population*1000000
#dat4$new_deaths <- dat4$cases_IHME_total/dat4$population*1000000

dat4_S1_202110 <- subset(dat4,dat4$gene== "S1"&dat4$month=="2021-10")

dat4_S1_202110 <- merge(dat4_S1_202110,region,by=c("country"),all.x = T)

selectname <- dat4_S1_202110 %>% filter(country %in% c("India","United Kingdom","Senegal","South Africa",
                                                       "Japan","United States","China"))

dat4_S1_202110$type <- NA
#dat4_S1_202209[which(dat4_S1_202209$continent=="Africa"),"type"] <- "Countries in Africa"
#dat4_S1_202209[-which(dat4_S1_202209$continent=="Africa"),"type"] <- "Countries outside Africa"
#country_high <- full_vero_2110[which(full_vero_2110$F_Vero>=75),"country"]
dat4_S1_202110[which(dat4_S1_202110$F_Vero>=75),"type"] <- "High coverage region"
dat4_S1_202110[-which(dat4_S1_202110$F_Vero<=25|dat4_S1_202110$F_Vero>=75),"type"] <- "Medium coverage region"
dat4_S1_202110[which(dat4_S1_202110$F_Vero<=25),"type"] <- "Low coverage region"

#write.csv(dat4_S1_202110,"vaccine_coverage_202110.csv",row.names = F)

dat4_S1_202110$continent <-factor(dat4_S1_202110$continent,
                                  levels= c('Africa','Asia','Europe','North America',
                                            'South America','Oceania',ordered=TRUE))
#dat4_S1_202209 <- dat4_S1_202209[-which(dat4_S1_202209$continent=="Oceania"),]
dat4_S1_202110$type <-factor(dat4_S1_202110$type,
                                levels= c('Low coverage region',"Medium coverage region",'High coverage region',
                                          ordered=TRUE))
p2 <- ggplot(dat4_S1_202110) +
  geom_point(aes(x = F_Vero, y = total_death_per_million,color=type), size=6,shape=19,alpha=0.75) +
  xlab("Full-dose vaccine coverage (%)") +
  ylab("Total reported mortality (× 10,000)") +
  scale_y_continuous(labels = c(expression(italic(0)),expression(''*'5'),expression(''*'10')),
                     expand = c(0,0),breaks = seq(0,100000,50000),limits = c(0,100000))+
  stat_cor(aes(x = F_Vero, y = total_death_per_million,color=type),method = "spearman",size=5)+
  stat_smooth(aes(x = F_Vero, y =total_death_per_million,color=type),method = "lm",se=F,size=1.5,
              linetype = "solid",formula = y ~ x)+
  geom_label_repel(data = selectname,
                   aes(x =  F_Vero, y = total_death_per_million, label = country), min.segment.length = 0,
                   seed = 43, box.padding = 0.4, segment.size = 0.5,
                   size = 5.5, label.size = NA, fill = "transparent") +
  #scale_x_continuous(limits=c(0,65),breaks = c(0,30,60))+
  mytheme+
  #theme(axis.title.y.left = element_text(size = 18,color="black",vjust=2,face = "bold", margin = margin(0,0.4,0,0,'cm')),
  #      axis.title.x = element_text(size = 18, color="black",face = "bold",margin = margin(0.4,0,0,0,'cm')))+
  #guides(fill = "none",color="none")+
  geom_vline(aes(xintercept=25),lty='solid',colour ="black",size=0.5)+
  geom_vline(aes(xintercept=50),lty='solid',colour ="black",size=0.5)+
  geom_vline(aes(xintercept=75),lty='solid',colour ="black",size=0.5)+
  scale_colour_manual("",values=c(c("#d4b9da","#c25b88","#5c4591")))+
  scale_fill_manual("",values=c(c( "#d4b9da","#c25b88","#5c4591")))+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+
  geom_hline(aes(yintercept=1),lty='dashed',colour ="black",size=1)+
  theme(legend.position = c(0.75,0.78))+
  guides(colour =guide_legend(nrow=6))

p2

fig1<-ggarrange(p1,p2,ncol =2, nrow = 1,widths = c(0.68,0.35),
                labels = c("C","D"),  vjust=1,align = "hv",
                font.label = list(size = 24, face = "bold"))
fig1
ggsave(fig1,filename = file.path(paste0(outpath,"fig1-c-d.pdf")),
       width = 17.5,height =6)

