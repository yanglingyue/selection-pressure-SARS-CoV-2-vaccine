data$case1 <- log(data$cases_IHME_month)
data$case2 <- log(data$new_cases_month)
data$travel2 <- log(data$travel_route)
data<- data[order(data$country,data$month,decreasing=F),]

# Set maximum lag
nlag = 3

colnames(data)
# Setting the lag for public health and social measure (PHSM) #
lag_npi <- tsModel::Lag(data$NPI_NV, group = data$country_index, k = 0:nlag)  

# Setting the lag for international travel#  
lag_travel <- tsModel::Lag(data$travel2, group = data$country_index, k = 0:nlag)     


# Setting the lag for adjusted vaccine coverage (Consider the waning of booster)#  
lag_FWBW_vero <- tsModel::Lag(data$FW_BW_P, group = data$country_index, k = 0:nlag)

# Setting the lag for adjusted vaccine coverage 
#Consider the difference in vaccine efficacy between Omicron and other VOCs (Menegale et al., 2023)#  
lag_FWBW_jama_vero <- tsModel::Lag(data$FW_BW_P_jama, group = data$country_index, k = 0:nlag)

# Setting the lag for the monthly new confirmed cases# 
lag_cases <- tsModel::Lag(data$case2, group = data$country_index, k = 0:nlag) 

# Setting the lag for the monthly new estimated infections based on IHME# 
lag_cases_IHME <- tsModel::Lag(data$case1, group = data$country_index, k = 0:nlag) 


# Setting the lag for natural immunity, which is calculated based on the estimated infections from IHME
# Considering the average effectiveness of previous infection in preventing reinfection for the Alpha, Beta, Delta, and Omicron variants
lag_NIW1 <- tsModel::Lag(data$NI_IHME_P, group = data$country_index, k = 0:nlag) 

# Setting the lag for natural immunity, which is calculated based on the estimated infections from IHME
# Considering the differences in effectiveness of previous infection in preventing reinfection between omicron and other VOCs 
lag_NIW2 <- tsModel::Lag(data$NI_IHME_P_omicron, group = data$country_index, k = 0:nlag) 



# Define cross basis matrix (combining nonlinear exposure and hysteresis function)
# Set the lag section

var <- lag_npi
npi_cb <- crossbasis(var,
                        argvar = list(fun = "lin"),
                        arglag = list(fun = "ns",df=2))
summary(npi_cb)


var <- lag_FWBW_vero
FBW_cb <- crossbasis(var,
                     argvar = list(fun = "poly",degree=3),
                        arglag = list(fun = "ns",df=2))
summary(FBW_cb)


var <- lag_FWBW_jama_vero
FBW_jama_cb <- crossbasis(var,
                     argvar = list(fun = "poly",degree=3),
                     arglag = list(fun = "ns",df=2))
summary(FBW_jama_cb)




var <- lag_NIW1
NIW1_cb <- crossbasis(var,
                     argvar = list(fun = "poly",degree=3),
                     arglag = list(fun = "ns",df=2))
summary(NIW1_cb)


var <- lag_NIW2
NIW2_cb <- crossbasis(var,
                     argvar = list(fun = "poly",degree=3),
                     arglag = list(fun = "ns",df=2))
summary(NIW2_cb)


var <- lag_travel
travel_cb <- crossbasis(var,
                         argvar = list(fun = "lin"),
                         arglag = list(fun = "ns",df=2))
summary(travel_cb)


# Specifies a unique column name for the gam() model
# Note: GLM (), GAM () or GLM.nb () models are not required
#colnames(npi_cb) = paste0("npi_cb.", colnames(npi_cb))
#colnames(vero_cb) = paste0("vero_cb.", colnames(vero_cb))
#colnames(T_vero_cb) = paste0("T_vero_cb.", colnames(T_vero_cb))
#colnames(F_vero_cb) = paste0("F_vero_cb.", colnames(F_vero_cb))
#colnames(shannon_cb) = paste0("shannon_cb.", colnames(npi_cb))



