setwd("C:/Users/Yanglingyue/Desktop/mutation_dnds/low_coverage_sequence/Covid-19_selection_vaccine")
rm(list = ls())
library(geofacet);library(sf);library(dplyr)
library(ggplot2);library(tidyverse);library(RColorBrewer)
library(showtext)

#############################sequence###################
dir.create("04_extended_analysis/plot_output/Figure_S2")
outpath <-  "04_extended_analysis/plot_output/Figure_S2/"

seq_dat <- read.csv("02_sup-data/sequence_number_metadata_13243751.csv")
colnames(seq_dat) <- c("region","freq")

population <- read.csv("02_sup-data/population_latest.csv")
colnames(population)[1] <- "region"

population[which(population$region =="Bonaire Sint Eustatius and Saba"),"region"] <- "Bonaire"
population[which(population$region =="Cape Verde"),"region"] <- "Cabo Verde"
population[which(population$region =="Czechia"),"region"] <- "Czech Republic"
population[which(population$region =="Democratic Republic of Congo"),"region"] <- "Democratic Republic of the Congo"
population[which(population$region =="Faeroe Islands"),"region"] <- "Faroe Islands"
population[which(population$region =="Congo"),"region"] <- "Republic of the Congo"
population[which(population$region =="Saint Martin (French part)"),"region"] <- "Saint Martin"
population[which(population$region =="Sint Maarten (Dutch part)"),"region"] <- "Sint Maarten"
population[which(population$region =="Bahamas"),"region"] <- "The Bahamas"
population[which(population$region =="Timor"),"region"] <- "Timor-Leste"
population[which(population$region =="United States Virgin Islands"),"region"] <- "U.S. Virgin Islands"
population[which(population$region =="United States"),"region"] <- "USA"
population[which(population$region =="Wallis and Futuna"),"region"] <- "Wallis and Futuna Islands"
population[which(population$region =="Micronesia (country)"),"region"] <- "Micronesia"


sequence1 <- merge(population,seq_dat,by=c("region"),all.y = T)
sequence1$prop <- (sequence1$freq/sequence1$population)*100


######https://gist.github.com/tadast/8827699#file-countries_codes_and_coordinates-csv#######
grid <- read.csv("02_sup-data/country_latitude_longitude.csv")
colnames(grid)[1] <- "region"

grid[which(grid$region =="Korea, Republic of"),"region"] <- "South Korea"
grid[which(grid$region =="Czech Republic"),"region"] <- "Czechia"
grid[which(grid$region =="Iran, Islamic Republic of"),"region"] <- "Iran"
grid[which(grid$region =="Russian Federation"),"region"] <- "Russia"
grid[which(grid$region =="United States"),"region"] <- "USA"
grid[which(grid$region =="Cape Verde"),"region"] <- "Cabo Verde"
grid[which(grid$region =="Bolivia, Plurinational State of"),"region"] <- "Bolivia"
grid[which(grid$region =="Côte d'Ivoire"),"region"] <- "Cote d'Ivoire"
grid[which(grid$region =="Czechia"),"region"] <- "Czech Republic"
grid[which(grid$region =="Swaziland"),"region"] <- "Eswatini"
grid[which(grid$region =="Congo, the Democratic Republic of the"),"region"] <- "Democratic Republic of the Congo"
grid[which(grid$region =="Libyan Arab Jamahiriya"),"region"] <- "Libya"
grid[which(grid$region =="Micronesia, Federated States of"),"region"] <- "Micronesia"
grid[which(grid$region =="Palestinian Territory, Occupied"),"region"] <- "Palestine"
grid[which(grid$region =="Congo"),"region"] <- "Republic of the Congo"
grid[which(grid$region =="Syrian Arab Republic"),"region"] <- "Syria"
grid[which(grid$region =="Taiwan, Province of China"),"region"] <- "Taiwan"
grid[which(grid$region =="Tanzania, United Republic of"),"region"] <- "Tanzania"
grid[which(grid$region =="Bahamas"),"region"] <- "The Bahamas"
grid[which(grid$region =="Virgin Islands, U.S."),"region"] <- "U.S. Virgin Islands"
grid[which(grid$region =="Venezuela, Bolivarian Republic of"),"region"] <- "Venezuela"
grid[which(grid$region =="Wallis and Futuna"),"region"] <- "Wallis and Futuna Islands"
grid[which(grid$region =="Virgin Islands, British"),"region"] <- "British Virgin Islands"
grid[which(grid$region =="Moldova, Republic of"),"region"] <- "Moldova"
grid[which(grid$region =="Brunei Darussalam"),"region"] <- "Brunei"
grid[which(grid$region =="Réunion"),"region"] <- "Reunion"
grid[which(grid$region =="Lao People's Democratic Republic"),"region"] <- "Laos"


sequence <- merge(grid,sequence1,by=c("region"),all.y = T)
colnames(sequence)[2:6] <- c("code1","code2","code3","latitude","longitude")



worldMap <- map_data("world")
#str(shannon_country_grid)


#######TEST#####
sequence[which(sequence$region =="Wallis and Futuna Islands"),"region"] <- "Wallis and Futuna"
#sequence[which(sequence$region =="United States"),"region"] <- "USA"
sequence[which(sequence$region =="United Kingdom"),"region"] <- "UK"
sequence[which(sequence$region =="The Bahamas"),"region"] <- "Bahamas"
sequence[which(sequence$region =="Republic of the Congo"),"region"] <- "Republic of Congo"
#sequence[which(sequence$region =="Hong Kong"),"region"] <- "??????"
#sequence[which(sequence$region =="Gibraltar"),"region"] <- "??????"
sequence[which(sequence$region =="Eswatini"),"region"] <- "Swaziland"
#sequence[which(sequence$region =="Crimea"),"region"] <- "??????"
sequence[which(sequence$region =="Cote d'Ivoire"),"region"] <- "Ivory Coast"
sequence[which(sequence$region =="Cabo Verde"),"region"] <- "Cape Verde"
##sequence[which(sequence$region =="U.S. Virgin Islands"),"region"] <- "Virgin Islands"
sequence[which(sequence$region =="British Virgin Islands"),"freq"] <- sum(sequence[which(sequence$region =="British Virgin Islands"),"freq"]+
                                                                            sequence[which(sequence$region =="U.S. Virgin Islands"),"freq"])
sequence[which(sequence$region =="British Virgin Islands"),"region"] <- "Virgin Islands"
worldMap[which(worldMap$region =="Trinidad"),"region"] <- "Trinidad and Tobago"
worldMap[which(worldMap$region =="Tobago"),"region"] <- "Trinidad and Tobago"
worldMap[which(worldMap$region =="Saint Vincent"),"region"] <- "Saint Vincent and the Grenadines"
worldMap[which(worldMap$region =="Grenadines"),"region"] <- "Saint Vincent and the Grenadines"

#worldMap[which(worldMap$region =="Czech Republic"),"region"] <- "Czechia"

#country_name <- as.data.frame(unique(worldMap$region))
#colnames(country_name) <- "region"
#country_name$type <- "A"
#seq11_test <-  merge(country_name,sequence,by=c("region"),all.y = T)




sequence$region <- as.character(sequence$region)


# seq_num<-cut(sequence$freq,breaks=c(0,100,1000,10000,100000,max(sequence$freq)),include.lowest=T,
#              labels = c("0-99","100-999","1000-9999","10000-99999",">100000"))
seq_num<-cut(sequence$freq,breaks=c(0,100,1000,10000,100000,1000000,max(sequence$freq)),include.lowest=T,
             labels = c("0-99","100-999","1000-9999","10000-99999","100000-9999999","≥1000000"))
num<-cut(sequence$freq,breaks=c(0,100,1000,10000,100000,1000000,max(sequence$freq,na.rm = T)),include.lowest=T,
         labels = c("10","100","1000","10000","100000","1000000"))#
sequence<-cbind(sequence,seq_num,num)
world.seq.num <- left_join(worldMap,sequence, by = "region")
#####Do not include HongKong, Gibraltar and Crimea#######################
world.seq.num$freq <- as.numeric(world.seq.num$freq)
world.seq.num[which(world.seq.num$region =="Taiwan"),"seq_num"] <-"1000-9999"
world.seq.num1 <- world.seq.num

world.seq.num <- world.seq.num[-which(world.seq.num$region =="Antarctica"),] 




pal <- rev(brewer.pal(11, "PRGn"))
levels <- pretty(c(0, 3), 5)
col1 <- colorRampPalette(pal[1:6])
col2 <- colorRampPalette(pal[8:11])
cols <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))

mytheme <-  theme(panel.grid.major.x=element_blank(),
                  panel.grid.minor.x=element_blank(),
                  panel.grid.major.y=element_blank(),
                  panel.grid.minor.y=element_blank(),
                  axis.ticks = element_blank(),
                  axis.line = element_blank(),
                  axis.text = element_blank(),
                  legend.justification=c(0,0),
                  # legend.position = c(0.00,0.06),
                  legend.position = "bottom",
                  legend.background = element_blank(),
                  #legend.key.size = unit(0.5,"cm"),
                  #legend.key.width = unit(0.5,"cm"),
                  #legend.key.height  = unit(0.5,"cm"),
                  legend.title = element_text(size = 10,face = "bold"),
                  legend.text = element_text(size = 10),
                  #legend.spacing = unit(6, "cm"), legend.spacing.x = NULL, # the spacing between legends (unit)
                  # legend.spacing.y = unit(1,"cm"),#the spacing between legends (unit)
                  axis.title = element_blank(),
                  plot.margin =  margin(0.5, 0, 0, 0, "cm"),
                  panel.background = element_rect(fill = "white"),
                  legend.box.background =element_blank(),
                  legend.box.margin=  margin(0, 0, 0, 0, "cm"),
                  panel.spacing = unit(0,"cm"),
                  panel.border = element_rect(fill='transparent',colour="transparent"),
                  legend.key = element_blank(),
                  plot.title = element_text(hjust = 0,size=12,face="bold"))


# c1 = rgb(165,0,38, maxColorValue = 255)
# c2 = rgb(215,48,39, maxColorValue = 255)
# c3 = rgb(244,109,67, maxColorValue = 255)
# c4 = rgb(253,174,97, maxColorValue = 255)
# c5 = rgb(254,224,139, maxColorValue = 255)
# c6 = rgb(255,255,191, maxColorValue = 255)

c1 = rgb(158,1,66, maxColorValue = 255)
c2 = rgb(213,62,79, maxColorValue = 255)
c3 = rgb(244,109,67, maxColorValue = 255)
c4 = rgb(253,174,97, maxColorValue = 255)
c5 = rgb(254,224,139, maxColorValue = 255)
c6 = rgb(255,255,191, maxColorValue = 255)
c0 = rgb(101,10,10, maxColorValue = 255)

sequence$freq <- as.numeric(sequence$freq)
colnames(sequence)[13] <- "Sequences"

p1 <- ggplot(world.seq.num, aes(map_id = region, fill = prop))+
  geom_map(map = world.seq.num,  color = "azure4",size=0.25)+
  coord_quickmap()+
  expand_limits(x = world.seq.num$long, y = world.seq.num$lat,expand = c(0, 0))+
  scale_fill_distiller("Numbers of sequence per 100 people",palette = "OrRd",limits=c(0,10),
                       breaks=c(0,2.5,5,7.5,10),
                       labels=c("0","2.5","5","7.5","10"),
                       direction = 1,na.value="grey88") + 
  geom_point(data = sequence,aes(x=longitude, y=latitude,group=region,size=Sequences,
                                 color=Sequences),alpha=0.7,
            color='#4575b4') +
  mytheme +
  theme(legend.key.size = unit(0.5,"cm"),legend.key.width = unit(1.35,"cm"),legend.key.height  = unit(0.35,"cm"))
  #guides(color = guide_legend())

p1





ggsave(p1,filename = file.path(outpath, "Fig_S2.pdf"),
       width = 9,height =4.2)



