setwd("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine")
rm(list = ls())
library(geofacet);library(sf);library(scatterpie)
dat<-read.csv("01_data/subsample_adaptation_vaccine_low.csv")
dat1 <- dat

source("load_package.R")
source("plot_theme.R")
source("load_data_clean.R")

dir.create("03_output/figure_1")
outpath <-  "03_output/figure_1/"


country1 <- as.data.frame(unique(dat1$country))
colnames(country1) <- "Country"
country1$country <- as.character(country1$Country)



#############################sequence###################
vero = subset(dat1,dat1$month== "2021-10"&dat1$gene=="S")
vero <- vero[,c("country","F_Vero")]
vero <- as.data.frame(vero)
vero$country <- as.character(vero$country)


######https://gist.github.com/tadast/8827699#file-countries_codes_and_coordinates-csv#######
grid <- read.csv("02_sup-data/country_latitude_longitude.csv")
grid[which(grid$Country =="Korea, Republic of"),"Country"] <- "South Korea"
grid[which(grid$Country =="Czech Republic"),"Country"] <- "Czechia"
grid[which(grid$Country =="Iran, Islamic Republic of"),"Country"] <- "Iran"
grid[which(grid$Country =="Russian Federation"),"Country"] <- "Russia"
country_grid <- merge(grid,country1,by=c("Country"),all.y = T)
colnames(country_grid) <- c("country","code1","code2","code3","lat","long","note")


vero_prop <- read.csv("02_sup-data/country_dnds_vero_prop_low_202110.csv")
x2 <- which(vero_prop$country=="Curacao"|vero_prop$country=="North Macedonia"|
              vero_prop$country=="Slovakia")
vero_prop <- vero_prop[-c(x2),]


vero_prop1 <- merge(country_grid,vero_prop,by=c("country"),all.y=T)
vero_prop1 <- merge(vero,vero_prop1,by=c("country"),all.y=T)
#vero_prop1$radius <-(vero_prop1$F_Vero)*0.05
vero_prop1$radius <- 3.5



colnames(vero) <- c("region","freq")
vero_region<-cut(vero$freq,breaks=c(0,25,75,1000),include.lowest=T,
              labels = c("Low coverage region","Medium coverage region","High coverage region"))#
vero<-cbind(vero,vero_region)
vero[which(vero$region =="United States"),"region"] <- "USA"
vero[which(vero$region =="United Kingdom"),"region"] <- "UK"


worldMap <- map_data("world")
worldMap[which(worldMap$region =="Czech Republic"),"region"] <- "Czechia"
world.seq.num <- left_join( worldMap,vero, by = "region")
world.seq.num$freq <- as.numeric(world.seq.num$freq)
world.seq.num1 <- world.seq.num
world.seq.num <- world.seq.num[-which(world.seq.num$region =="Antarctica"),] 





theme1 <-  theme(axis.ticks = element_blank(),
                   axis.line = element_blank(),
                   axis.text = element_blank(),
                   legend.justification=c(0.5,0),
                   legend.position = c("bottom"),
                   legend.background = element_blank(),
                   legend.key.size = unit(0.5,"cm"),
                   legend.key.width = unit(1,"cm"),
                   legend.key.height  = unit(0.35,"cm"),
                   legend.title = element_text(size = 18,face = "bold",vjust = 0.85,hjust=2),
                   legend.text = element_text(size = 18,lineheight=6),
                   legend.spacing = unit(6, "cm"), legend.spacing.x = NULL, # the spacing between legends (unit)
                   legend.spacing.y = NULL,#the spacing between legends (unit)
                   axis.title = element_blank(),
                   plot.margin =  margin(0, 0, 0, 0, "cm"),
                   panel.background = element_rect(fill = "white"),
                   legend.box.background =element_blank(),
                   legend.box.margin=  margin(0, 0, 0, 0, "cm"),
                   panel.spacing = unit(0,"cm"),
                   panel.border = element_rect(fill='transparent',colour="transparent"),
                   legend.key = element_blank(),
                   plot.title = element_text(hjust = 0.06, size = 18, vjust = 0))
  
#colors <- c( '#4575b4','#fdae61', '#d73027', '#f46d43','#DCDCDC','#fee090', '#e0f3f8' )  
colors <- c("#5f1c83","#6cbb73","#425da2","#B2182B","#f7b455" ,"#d06a6c")  
colors <- c("#1a3a5c","#77afd0","#425da2","#B2182B","#f7b455" ,"#d06a6c")  
colors <- c("#1a3a5c","#ABD9E9","#5480b9","#d06a6c","#ffdc75","#f29f5f")  
colors <- c("#516084","#cbe9f1","#79b0d2",'#1b9e77',"#ea7449",'#7570b3')  

p1 <- ggplot(world.seq.num, aes(map_id = region,fill=vero_region))+
    geom_map(map = world.seq.num,  color = "azure4",alpha=1,size=0.25)+
    coord_quickmap()+
    expand_limits(x = world.seq.num$long, y = world.seq.num$lat,expand = c(0, 0))+
    geom_scatterpie(aes(x=long, y=lat, group=country, r=radius), data=vero_prop1,
                    cols=c("Vero_mRNA" ,"Vero_In", "Vero_AD"), color="#3B3B3BFF")+
    #geom_scatterpie_legend(vero_prop1$radius, x = -160, y = -55, n = 4, 
    #                     labeller = function(x) x/0.05 )+
    geom_rect(aes(xmin = -16, xmax = 45.5, ymin = 36, ymax = 72), color = "black",size=0.35, fill = NA)  +
    scale_fill_manual("",values =  colors, na.value="grey88")+
    theme1+
    #guides(color = guide_legend())+
    theme(legend.key.size = unit(0.5,"cm"),legend.key.width = unit(1,"cm"),legend.key.height  = unit(0.35,"cm"))+
    theme(panel.border = element_rect(color = "black", size = 1, fill = NA))
    #scale_fill_manual("",values =  colors,
    #                na.value="grey88",labels = c("Low coverage region","Medium coverage region",
    #                                             "High coverage region","Not included/unknown"))+
    #scale_fill_manual(name = "",  values = colors[1:6],labels=c('mRNA vaccine','Inactivated vaccine',"Adenovirus vector vaccine"))
p1
  
  
#ggsave(p1,filename = file.path(outpath, "fig_1_a.pdf"),
#       width = 25,height =15)


  
  
colors <- c("#cbe9f1","#79b0d2","#516084","#eb4c1c","#f7b455","#6cbb73" )  

p2 <- ggplot(world.seq.num, aes(map_id = region,fill=vero_region))+
    geom_map(map = world.seq.num,  color = "azure4",alpha=1,size=0.25)+
    expand_limits(x = world.seq.num$long, y = world.seq.num$lat,expand = c(0, 0))+
    theme1+
    guides(color = guide_legend())+
    theme(legend.position = "none")+
    geom_rect(aes(xmin = -16, xmax = 45.5, ymin = 36, ymax = 72), color = "black",size=0.6, fill = NA)  +
    scale_fill_manual("",values =  colors, na.value="grey88")+
    coord_sf(xlim = c(-16,46), ylim = c(36,72), expand = FALSE)
p2

ggsave(p2,filename = file.path(outpath, "fig_1_A_s.pdf"),
         width = 3.8,height =3.8)
  
 
vero$region <- factor(vero$region,
                      levels = vero$region[order(vero$freq, decreasing = T)],ordered=TRUE)

source("plot_theme.R")


p3 <- ggplot(vero,aes(x=region,y=freq,group=1))+
  geom_bar(aes(x=region,y=freq,fill=vero_region),stat="identity",alpha=1)+
  scale_fill_manual("",values=colors)+
  xlab("") +
  ylab("Full-dose vaccine coverage (%)") +
  scale_y_continuous(expand = c(0,0))+
  guides(fill = "none",color="none")+
  geom_hline(aes(yintercept=25),lty='dashed',colour ="grey25",size=0.5)+
  geom_hline(aes(yintercept=50),lty='dashed',colour ="grey25",size=0.5)+
  geom_hline(aes(yintercept=75),lty='dashed',colour ="grey25",size=0.5)+
  #Rotating labels
  coord_flip()+
  mytheme+  
  theme(axis.text.y.left = element_blank(),axis.text.y.right = element_blank(),
        axis.ticks.y.left = element_blank(),axis.ticks.y.right = element_blank(),
        axis.title.y.left = element_blank(),axis.title.y.right  = element_blank())+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())
#ggtitle("Natural immunity")+
p3 


fig2<-ggarrange(p1,p3,ncol =2, nrow = 1,widths = c(0.85,0.15),
                labels = c("",""),  vjust=1,align = "hv",
                font.label = list(size = 24, face = "bold"))
fig2

ggsave(fig2,filename = file.path(outpath, "fig_1_A_B.pdf"),
       width = 24,height =8)


library(rstatix)




