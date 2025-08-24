# # install.packages("lubridate")
rm(list = ls())
setwd("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine")
dat<-read.csv("01_data/subsample_adaptation_vaccine_low.csv")
dat1 <- dat
library(rstatix)

source("load_package.R")
source("plot_theme.R")
source("load_data_clean.R")

dir.create("03_output/figure_1")
outpath <-  "03_output/figure_1/"
#Filter out the countries with a full-dose vaccine coverage rate exceeding 75% as of September 2022.
full_vero_2110 <- dat1[which(dat1$gene=="S1"&dat1$month=="2021-10"),]
country_high <- full_vero_2110[which(full_vero_2110$F_Vero>=75),"country"]

#Filter out countries with a full-dose COVID-19 vaccine coverage rate below 25% as of September 2022.
country_low <- full_vero_2110[which(full_vero_2110$F_Vero<=25),"country"]

#Filter out countries with a full-dose COVID-19 vaccine coverage rate between 25% and 75% as of September 2022.
country_mid <- full_vero_2110[which(full_vero_2110$F_Vero>25&full_vero_2110$F_Vero<75),"country"]




dat2= subset(dat1,dat1$gene== "S1")
dat2$Income_group <- NA
dat2$country <- as.character(dat2$country)
dat2[which(dat2$country %in% country_high$country),"Income_group"] <- "High coverage region"
dat2[which(dat2$country %in% country_mid$country),"Income_group"] <- "Medium coverage region"
dat2[which(dat2$country %in% country_low$country),"Income_group"] <- "Low coverage region"

dat2$Income_group <-factor(dat2$Income_group,levels= c('Low coverage region',
                                                       'Medium coverage region','High coverage region',ordered=TRUE))
#dat2 <- dat2[-which(is.na(dat2$Income_group)),]


var_type <- c("dnds_mean_focal","muts_nonsyn_focal","fitness_focal")


month_name <- unique(dat1$month)
print(month_name)
m =length(month_name)


dat_final <- matrix(NA,m*length(var_type),5)
dat_final <- as.data.frame(dat_final)


for (s in 1:length(var_type))
{
  print(var_type[s])
 
  var <- c("month","country",var_type[s])
  ourdata <-dat2[,var] 
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
  
  assign(paste("df_", var_type[s], sep=""),df_ci)
  
  
  df_ci$var_type <- var_type[s]
  #df_ci$var_type <- as.character(df_ci$var_type)
  
  dat_final[(c(31*s-30):(31*s)),] <- df_ci
  colnames(dat_final) <- colnames(df_ci)
  
}
colnames(df_dnds_mean_focal)[1:3] <- c("dnds_lci","dnds_mean","dnds_uci")
colnames(df_fitness_focal)[1:3] <- c("fitness_lci","fitness_mean","fitness_uci")
colnames(df_muts_nonsyn_focal)[1:3] <- c("nonsyn_lci","nonsyn_mean","nonsyn_uci")


dat_final1 <- merge(df_dnds_mean_focal,df_fitness_focal,by=c("month"),all.x = T)
dat_final2 <- merge(dat_final1,df_muts_nonsyn_focal,by=c("month"),all.x = T)

#cols <- c('#1b9e77',"#c25b88",'#7570b3')
cols <- c("#c25b88","#849f84","#cfb848")

p1<-ggplot(dat_final2,aes(x=month,y=dnds_mean,group=1)) +
  #geom_rect(aes(xmin="2021-11", xmax="2022-09", ymin=-Inf, ymax=Inf),fill='grey80',alpha = 0.01)+
  geom_line(aes(x=month,y=dnds_mean,color="Nonsynonymous to synonymous divergence ratio"),size=0.8,alpha=0.9)+
  geom_ribbon(aes(x=month,ymin = dnds_lci, ymax = dnds_uci,fill="Nonsynonymous to synonymous divergence ratio"),alpha=0.55) +
  geom_point(aes(x=month,y=dnds_mean,color="Nonsynonymous to synonymous divergence ratio"), size=3,shape=19,alpha=0.75) +
  geom_line(aes(x=month,y=fitness_mean,color="Mutational fitness"),size=0.8,alpha=0.9)+
  geom_ribbon(aes(x=month,ymin = fitness_lci, ymax = fitness_uci,fill="Mutational fitness"),alpha=0.55) +
  geom_point(aes(x=month,y=fitness_mean,color="Mutational fitness"), size=3,shape=19,alpha=0.75) +
  
  geom_line(aes(x=month,y=nonsyn_mean*max(dnds_mean,na.rm = T)*0.06,color="Number of nonsynonymous mutations"),size=0.8,alpha=0.9)+
  geom_ribbon(aes(x=month,ymin = nonsyn_lci*max(dnds_mean,na.rm = T)*0.06, 
                  ymax = nonsyn_uci*max(dnds_mean,na.rm = T)*0.06,fill="Number of nonsynonymous mutations"),alpha=0.55) +
  geom_point(aes(x=month,y=nonsyn_mean*max(dnds_mean,na.rm = T)*0.06,color="Number of nonsynonymous mutations"), size=3,shape=19,alpha=0.75) +
  annotate(geom = "text", x = "2021-01", y =4.5,angle=90, label = "Gamma", hjust = 0,size =5)+
  annotate(geom = "text", x = "2021-05", y =4.5,angle=90, label = "Delta", hjust =0,size = 5)+
  annotate(geom = "text", x = "2020-12", y =4.5,angle=90, label = "Alpha/Beta", hjust = 0,size = 5)+
  annotate(geom = "text", x = "2021-11", y =4.5,angle=90, label = "Omicron", hjust =0,size = 5)+
  annotate(geom = "text", x = "2021-03", y =7.5,angle=0, label = "Vaccine rollout", hjust =0,size = 6.5)+
  geom_vline(aes(xintercept="2021-01"),lty='solid',colour ="black",size=0.5)+
  geom_vline(aes(xintercept="2021-05"),lty='solid',colour ="black",size=0.5)+
  geom_vline(aes(xintercept="2020-12"),lty='solid',colour ="black",size=0.5)+
  geom_vline(aes(xintercept="2021-11"),lty='solid',colour ="black",size=0.5)+
  geom_hline(aes(yintercept=1),lty='dashed',colour ="black",size=1)+
  scale_y_continuous(sec.axis = sec_axis(~./(max(dat_final2$dnds_mean,na.rm = T)*0.06),
                                         name="Number of nonsynonymous \nmutations", breaks = seq(0, 28, 7)),
                     breaks = seq(0, 8, 2))+
  
  #facet_wrap(~var_type,scales = "free",ncol=2)+
  labs(x="",
       y= "Nonsynonymous to synonymous\n divergence ratio / Mutational fitness")+
  scale_x_discrete(breaks=c("2020-04","2020-08","2020-12","2021-04",
                            "2021-08","2021-12","2022-04","2022-08"),
                   labels=c("Apr\n2020","Aug\n2020","Dec\n2020","Apr\n2021",
                            "Aug\n2021","Dec\n2021", "Apr\n2022","Aug\n2022"))+
  scale_colour_manual(values=cols,"")+
  scale_fill_manual(values=cols,"")+
  mytheme+
  theme(panel.grid=element_line(linetype=2,color='grey'))+
  #theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
  #panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+
  #theme(legend.position = c(0.1,0.78))+
  guides(color=guide_legend(nrow=3,byrow = T))

p1  

ggsave(p1,filename = file.path(paste0(outpath,"fig1-C.pdf")),
       width = 12.9,height =7)






dat3 <- dat2[-which(dat2$month=="2021-11"|dat2$month=="2021-12"|dat2$month=="2022-01"|dat2$month=="2022-02"|
                      dat2$month=="2022-03"|dat2$month=="2022-04"|dat2$month=="2022-05"|dat2$month=="2022-06"|
                      dat2$month=="2022-07"|dat2$month=="2022-08"|dat2$month=="2022-09"), ]

df.summary <- dat3 %>%
  group_by(Income_group) %>%
  summarise(
    sd = sd(dnds_mean_focal, na.rm = TRUE),
    dnds_mean_focal = mean(dnds_mean_focal, na.rm = TRUE)
  )
df.summary


cols <- c("#d4b9da","#9370ae","#494379")

cols <- c("#849f84","#849f84","#849f84")
#cols <- c("#ABD9E9","#5480b9","#1a3a5c")


p2 <- ggplot(dat3, aes(x=Income_group, y=dnds_mean_focal,fill=Income_group,color=Income_group)) + 
      #geom_jitter( position = position_jitter(0.2),alpha=0.35) + 
      #geom_line(aes(group = 1), data = df.summary,color = "black") +
      geom_violin(alpha = 0.55,trim=F,width=0.6,color="#849f84") +
      #geom_boxplot(notch=F, outlier.shape = NA,outlier.size = 2,color="black",
      #         outlier.fill = "black",alpha=0.45, width = 0.15)+
      geom_errorbar(aes(ymin = dnds_mean_focal-sd, ymax = dnds_mean_focal+sd),
                data = df.summary, width = 0.15,color="black",size=1) +
      geom_point(aes(y = dnds_mean_focal),size=3,data = df.summary,  color = "black") +  
      labs(x="",
           y= "Nonsynonymous to synonymous\n divergence ratio")+
      #geom_crossbar(aes(y = dnds_mean_focal, ymin = dnds_mean_focal, ymax = dnds_mean_focal),
      #          data = df.summary, width = 0.2, fill = "white", color = "black") +  # Add this line
      scale_y_continuous(limits=c(0,10),breaks = seq(0,10,5))+
      scale_fill_manual(name = '',values = cols)+
      scale_color_manual(name = '',values = cols)+
      geom_signif(comparisons = list(c("Low coverage region", "High coverage region"),
                                 c("Low coverage region", "Medium coverage region"),
                                 c("Medium coverage region", "High coverage region")),
              map_signif_level=F, textsize=6, y_position = c(9.5,8.5,7.5),test = "wilcox.test",
              color="black")+
      mytheme+
      theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+
      theme(axis.text.x = element_text(angle = 30, hjust = 0.75,vjust = 0.75))+
      theme(plot.margin=unit(c(0.5,1,0.5,1),'lines'))+
      guides(fill = "none",color="none")
p2   
 


fig1<-ggarrange(p1,p2,ncol =2, nrow = 1,widths = c(0.68,0.35),
                labels = c("D","E"),  vjust=1,align = "hv",
                font.label = list(size = 24, face = "bold"))
fig1
#ggsave(fig1,filename = file.path(paste0(outpath,"fig1-C-D.pdf")),
#       width = 19,height =7)





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
#colors <- c("#b4b763","#7eb874","#3b8894")
colors <- c("#ABD9E9","#5480b9","#1a3a5c")

p3 <- ggplot(dat4_S1_202110) +
  geom_point(aes(x = F_Vero, y = total_death_per_million,color=type), size=6,shape=19,alpha=0.75) +
  xlab("Full-dose vaccine coverage (%)") +
  ylab("Total reported deaths per million\n (× 10,000)") +
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
  scale_colour_manual("",values=colors)+
  scale_fill_manual("",values=colors)+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+
  geom_hline(aes(yintercept=1),lty='dashed',colour ="black",size=1)+
  #theme(legend.position = c(0.75,0.78))+
  theme(plot.margin=unit(c(0.5,1,0.5,1),'lines'))+
  guides(colour =guide_legend(nrow=6))

p3


fig1<-ggarrange(p2,p3,ncol =2, nrow = 1,widths = c(1,1),
                labels = c("",""),  vjust=1,align = "hv",
                font.label = list(size = 24, face = "bold"))
fig1
ggsave(fig1,filename = file.path(paste0(outpath,"fig1-D-E.pdf")),
       width = 11.5,height =7)


#fig1<-ggarrange(p1,p2,ncol =2, nrow = 1,widths = c(1,1),
#                labels = c("C","D"),  vjust=1,align = "hv",
#                font.label = list(size = 24, face = "bold"))
#fig1
#ggsave(fig1,filename = file.path(paste0(outpath,"fig1-v1.pdf")),
#       width = 12,height =7)

