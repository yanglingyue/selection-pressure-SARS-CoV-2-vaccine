setwd("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine")
rm(list = ls())
library(geofacet);library(sf);library(scatterpie)
dat<-read.csv("01_data/subsample_dnds_vaccine_diversity_low_death.csv")
dat1 <- dat

source("load_package.R")
source("plot_theme.R")
source("load_data_clean.R")

dir.create("04_extended_analysis/plot_output/Figure_S5")
outpath <-  "04_extended_analysis/plot_output/Figure_S5/"



country1 <- as.data.frame(unique(dat1$country))
colnames(country1) <- "Country"
country1$Country <- as.character(country1$Country)


#############################sequence###################
#alldata <- read.csv("01_data/data_process/dnds_BW_NPI_V_low_input_2209.csv")
#alldata<- alldata[order(alldata$country,alldata$month,decreasing=F),]
vero = subset(dat1,dat1$month== "2021-10"&dat1$gene=="S")
vero <- vero[,c("country","F_Vero")]
vero$country <- as.character(vero$country)
colnames(vero) <- c("country","freq")
vero[which(vero$country =="United States"),"country"] <- "USA"
vero[which(vero$country =="United Kingdom"),"country"] <- "UK"
vero$region <- vero$country


######https://gist.github.com/tadast/8827699#file-countries_codes_and_coordinates-csv#######
grid <- read.csv("02_sup-data/country_latitude_longitude.csv")
grid[which(grid$Country =="Korea, Republic of"),"Country"] <- "South Korea"
grid[which(grid$Country =="Czech Republic"),"Country"] <- "Czechia"
grid[which(grid$Country =="Iran, Islamic Republic of"),"Country"] <- "Iran"
grid[which(grid$Country =="Russian Federation"),"Country"] <- "Russia"
country_grid <- merge(grid,country1,by=c("Country"),all.y = T)
country_grid$month <- "Average"
colnames(country_grid) <- c("Country","code1","code2","code3","latitude","longitude","note")

vero_prop <- read.csv("02_sup-data/country_dnds_vero_prop_low_202110.csv")
x2 <- which(vero_prop$country=="Curacao"|vero_prop$country=="North Macedonia"|
              vero_prop$country=="Slovakia")
vero_prop <- vero_prop[-c(x2),]


#colnames(vero_prop)[1] <- "Country"
colnames(grid) <- c("country","code1","code2","code3","lat","long")

vero_prop1 <- merge(grid,vero_prop,by=c("country"),all.y=T)
#vero_prop1 <- vero_prop1[,-c(2:4,7,8,9:17)]
vero_prop1 <- merge(vero,vero_prop1,by=c("country"),all.y=T)
#vero_prop1$radius <-vero_prop1$F_Vero*0.05
vero_prop1$radius <- 3.5


#world_map1 <- world_map[,c("featurecla","geometry")]
#colnames(world_map1)
worldMap <- map_data("world")
#str(shannon_country_grid)

world.seq.num <- left_join( worldMap,vero, by = "region")
world.seq.num$freq <- as.numeric(world.seq.num$freq)
world.seq.num1 <- world.seq.num
world.seq.num <- world.seq.num[-which(world.seq.num$region =="Antarctica"),] 
worldMap[which(worldMap$region =="Czech Republic"),"region"] <- "Czechia"




pal <- rev(brewer.pal(11, "PRGn"))
levels <- pretty(c(0, 3), 5)
col1 <- colorRampPalette(pal[1:6])
col2 <- colorRampPalette(pal[8:11])
cols <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))

#worldMap <- spTransform(worldMap, CRS("+proj=robin"))

  
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
colors <- c( '#7570b3','#1b9e77', '#f46d43','#DCDCDC','#fee090', '#e0f3f8' )  
  
p1 <- ggplot(world.seq.num, aes(map_id = region))+
    geom_map(map = world.seq.num,  color = "azure4", fill ='floralwhite',alpha=1,size=0.25)+
    coord_quickmap()+
    expand_limits(x = world.seq.num$long, y = world.seq.num$lat,expand = c(0, 0))+
    #geom_point(data = sequence_country_grid,aes(x=longitude, y=latitude),alpha = 0.7,color = "black",  size = 2) +
    scale_color_distiller("Vaccine coverage",palette = "GnBu",limits=c(0,100), breaks=c(0,25,50,75,100),
                          labels=c("0","25","50","75","100"),
                          direction = 1,na.value="grey88") + 
    theme1+
    guides(color = guide_legend())+
    theme(legend.key.size = unit(0.5,"cm"),legend.key.width = unit(1,"cm"),legend.key.height  = unit(0.35,"cm"))+
    geom_scatterpie(aes(x=long, y=lat, group=country, r=radius), data=vero_prop1,
                    cols=c("Vero_mRNA" ,"Vero_In", "Vero_AD"), color="#3B3B3BFF")+
    geom_rect(aes(xmin = -16, xmax = 45.5, ymin = 36, ymax = 72), color = "black",size=0.35, fill = NA)  +
    scale_fill_manual(name = "",  values = colors,labels=c('mRNA vaccine','Inactivated vaccine',"Adenovirus vector vaccine"))
p1
  
  
ggsave(p1,filename = file.path(outpath, "fig_S5_a.pdf"),
       width = 15,height =7)


  
  
  
p2 <- ggplot(world.seq.num, aes(map_id = region))+
    geom_map(map = world.seq.num,  color = "azure4", fill ='floralwhite',alpha=1,size=0.25)+
    #coord_quickmap()+
    expand_limits(x = world.seq.num$long, y = world.seq.num$lat,expand = c(0, 0))+
    #geom_point(data = sequence_country_grid,aes(x=longitude, y=latitude),alpha = 0.7,color = "black",  size = 2) +
    scale_color_distiller("Vaccine coverage",palette = "GnBu",limits=c(0,100), breaks=c(0,25,50,75,100),
                          labels=c("0","25","50","75","100"),
                          direction = 1,na.value="grey88") + 
    theme1+
    guides(color = guide_legend())+
    theme(legend.position = "none")+
    geom_scatterpie(aes(x=long, y=lat, group=country, r=1.75), data=vero_prop1,
                    cols=c("Vero_mRNA" ,"Vero_In", "Vero_AD"), color="#3B3B3BFF")+
    geom_rect(aes(xmin = -16, xmax = 45.5, ymin = 36, ymax = 72), color = "black",size=0.6, fill = NA)  +
    coord_sf(xlim = c(-16,46), ylim = c(36,72), expand = FALSE)+
    scale_fill_manual(name = "",  values = colors,labels=c('mRNA vaccine','Inactivated vaccine',"Adenovirus vector vaccine"))
p2

ggsave(p2,filename = file.path(outpath, "fig_S5_b.pdf"),
         width = 3.5,height =3.5)
  
  
library(cowplot)
  gg_inset_map = ggdraw() +
    draw_plot(p1) +
    draw_plot(p2, x = -0.00, y = 0.2, width = 0.23, height = 0.23)
  gg_inset_map
  

 
