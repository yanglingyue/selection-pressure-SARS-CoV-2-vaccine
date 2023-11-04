dat1$date <- parse_date_time(dat1$month, "ym")
dat1$date <-as.Date(as.POSIXct(dat1$date,tz="Hongkong"))
dat1$month1 <- month(dat1$date)
dat1$year <- year(dat1$date)

##Remove January and February 2020, as the submitted sequences of most countries are unavailable.
y1 <- which(dat1$month=="2020-01"|dat1$month=="2020-02") 
dat1 <-dat1[-y1,]



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

dat1$travel2 <- log(dat1$travel_route)
dat1$travel2<-as.numeric(dat1$travel2)

