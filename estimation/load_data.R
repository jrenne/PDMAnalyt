# ==============================================================================
# Get yield data
# ==============================================================================

library(fredr)
library(Hmisc)

indic.download <- 1
#indic.download <- 0

start.date <- "1968-10-01"
end.date   <- Sys.Date()
freq       <- "a"
frequency  <- 4

fredr_set_key("df65e14c054697a52b4511e77fcfa1f3")
start_date <- as.Date(start.date)
end_date   <- as.Date(end.date)
f <- function(ticker){
  fredr(series_id = ticker,
        observation_start = start_date,observation_end = end_date,
        frequency = freq,aggregation_method = "avg")
}

list.variables <- c("DTB3","GDP","M318501Q027NBEA",
                    "THREEFY1","THREEFY2","THREEFY3","THREEFY4","THREEFY5",
                    "THREEFY6","THREEFY7","THREEFY8","THREEFY9","THREEFY10",
                    "THREEFYTP5","THREEFYTP10",
                    "CPIAUCSL","GDPPOT","GDPC1","DFII5","DFII10","DGS30")

for(i in 1:length(list.variables)){
  data.var <- f(list.variables[i])
  eval(parse(text = gsub(" ","",paste("data.var.frame = data.frame(date=data.var$date,",
                                      list.variables[i],"=data.var$value)",
                                      sep=""))))
  if(i==1){
    DATA = data.var.frame
  }else{
    DATA = merge(DATA,data.var.frame,by="date",all=TRUE)
  }
}

# 4-quarter MA of surpluses:
DATA$surpl <- DATA$M318501Q027NBEA
T_temp <- length(DATA$surpl)
DATA$surpl[4:T_temp] <- DATA$surpl[4:T_temp] + DATA$surpl[3:(T_temp-1)] + 
  DATA$surpl[2:(T_temp-2)] + DATA$surpl[1:(T_temp-3)]
DATA$surpl <- DATA$surpl/DATA$GDP

DATA$z <- log(DATA$GDPC1/DATA$GDPPOT)
DATA$z <- log(DATA$GDPC1/DATA$GDPPOT) - mean(DATA$z,na.rm = TRUE)

lag <- 1
DATA$pi <- NaN
DATA$pi[(lag+1):dim(DATA)[1]] <- log(DATA$CPIAUCSL[(lag+1):dim(DATA)[1]]/
                                       DATA$CPIAUCSL[1:(dim(DATA)[1]-lag)])

DATA$dy <- NaN
DATA$dy[(lag+1):dim(DATA)[1]] <- log(DATA$GDPC1[(lag+1):dim(DATA)[1]]/
                                       DATA$GDPC1[1:(dim(DATA)[1]-lag)])

plot(DATA$date,DATA$z,type="l");grid()



# 
# 
# # Download Survey of Professional Forecasters data: ============================
# 
# # 10 year GDP growth------------------------------------------------------------
# if(indic.download==1){
#   download.file("https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/survey-of-professional-forecasters/data-files/files/mean_rgdp10_level.xlsx",
#                 "Data/mean_rgdp10_level.xlsx")
# }
# SPF <- readxl::read_xlsx(path="Data/mean_rgdp10_level.xlsx",na="#N/A")
# SPF$date <- as.Date(paste(SPF$YEAR,"-",1+3*(SPF$QUARTER-1),"-01",sep=""))
# indic1st <- which(SPF$date==start.date)
# SPF      <- SPF[indic1st:dim(SPF)[1],]
# SPF.GDP10<- data.frame(date=SPF$date,RGDP10=SPF$RGDP10)
# 
# # 1 year GDP growth-------------------------------------------------------------
# if(indic.download==1){
#   download.file("https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/survey-of-professional-forecasters/data-files/files/mean_rgdp_level.xlsxs",
#                 "Data/mean_rgdp_level.xlsxs")
# }
# SPF <- readxl::read_xlsx(path="Data/mean_rgdp_level.xlsxs")
# SPF$date <- as.Date(paste(SPF$YEAR,"-",1+3*(SPF$QUARTER-1),"-01",sep=""))
# SPF <- data.frame(date=SPF$date,RGDP1=100*log(as.numeric(SPF$RGDP6)/as.numeric(SPF$RGDP2)))
# indic1st <- which(SPF$date==start.date)
# SPF <- SPF[indic1st:dim(SPF)[1],]
# SPF.GDP<- data.frame(date=SPF$date,RGDP1=SPF$RGDP1)
# 
# # 1 year CPI -------------------------------------------------------------------
# if(indic.download==1){
#   download.file("https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/survey-of-professional-forecasters/data-files/files/mean_cpi_level.xlsx",
#                 "Data/mean_cpi_level.xlsx")
# }
# SPF <- readxl::read_xlsx(path="Data/mean_cpi_level.xlsx")
# SPF$date <- as.Date(paste(SPF$YEAR,"-",1+3*(SPF$QUARTER-1),"-01",sep=""))
# SPF <- data.frame(date=SPF$date,CPI1=as.numeric(SPF$CPI6))
# indic1st <- which(SPF$date==start.date)
# SPF <- SPF[indic1st:dim(SPF)[1],]
# SPF.CPI <- data.frame(date=SPF$date,CPI1=SPF$CPI1)
# 
# # 10 year CPI-------------------------------------------------------------------
# if(indic.download==1){
#   download.file("https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/survey-of-professional-forecasters/data-files/files/mean_cpi10_level.xlsx",
#                 "Data/mean_cpi10_level.xlsx")
# }
# SPF <- readxl::read_xlsx(path="Data/mean_cpi10_level.xlsx")
# SPF$date <- as.Date(paste(SPF$YEAR,"-",1+3*(SPF$QUARTER-1),"-01",sep=""))
# SPF <- data.frame(date=SPF$date,CPI10=as.numeric(SPF$CPI10))
# indic1st <- which(SPF$date==start.date)
# SPF <- SPF[indic1st:dim(SPF)[1],]
# SPF.CPI10<- data.frame(date=SPF$date,CPI10=SPF$CPI10)
# 
# # 10 year TBILL-----------------------------------------------------------------
# if(indic.download==1){
#   download.file("https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/survey-of-professional-forecasters/data-files/files/mean_bill10_level.xlsx",
#                 "Data/mean_bill10_level.xlsx")
# }
# SPF <- readxl::read_xlsx(path="Data/mean_bill10_level.xlsx")
# SPF$date <- as.Date(paste(SPF$YEAR,"-",1+3*(SPF$QUARTER-1),"-01",sep=""))
# SPF <- data.frame(date=SPF$date,BILL10=as.numeric(SPF$BILL10))
# indic1st <- which(SPF$date==start.date)
# SPF <- SPF[indic1st:dim(SPF)[1],]
# SPF.BILL10<- data.frame(date=SPF$date,BILL10=SPF$BILL10)
# 
# 



# Nominal yield data from GSW 2006:
if(indic.download==1){
  download.file("https://www.federalreserve.gov/data/yield-curve-tables/feds200628.csv",
                "Data/feds200628.csv")
}
DAT <- csv.get(file = "Data/feds200628.csv", skip = 9)
library(tidyverse)
library(zoo)
DAT$Date <- as.Date(DAT$Date)
DAT <- arrange(DAT,Date)
DAT$year <- format(DAT$Date,"%Y")
DAT_yearly <- DAT %>%
  group_by(year) %>%
  summarise_all(function(x){mean(x,na.rm=TRUE)})
days <- as.numeric(format(DAT_yearly$Date,"%d"))
DAT_yearly$Date <- as.Date(paste(format(DAT_yearly$Date,"%Y"),"-01-01",sep=""))
DAT_GSW_nom <- data.frame(date=DAT_yearly$Date,
                          SVENY01=DAT_yearly$SVENY01,
                          SVENY02=DAT_yearly$SVENY02,
                          SVENY03=DAT_yearly$SVENY03,
                          SVENY04=DAT_yearly$SVENY04,
                          SVENY05=DAT_yearly$SVENY05,
                          SVENY06=DAT_yearly$SVENY06,
                          SVENY07=DAT_yearly$SVENY07,
                          SVENY08=DAT_yearly$SVENY08,
                          SVENY09=DAT_yearly$SVENY09,
                          SVENY10=DAT_yearly$SVENY10,
                          SVENY11=DAT_yearly$SVENY11,
                          SVENY12=DAT_yearly$SVENY12,
                          SVENY13=DAT_yearly$SVENY13,
                          SVENY14=DAT_yearly$SVENY14,
                          SVENY15=DAT_yearly$SVENY15,
                          SVENY16=DAT_yearly$SVENY16,
                          SVENY17=DAT_yearly$SVENY17,
                          SVENY18=DAT_yearly$SVENY18,
                          SVENY19=DAT_yearly$SVENY19,
                          SVENY20=DAT_yearly$SVENY20,
                          SVENY21=DAT_yearly$SVENY21,
                          SVENY22=DAT_yearly$SVENY22,
                          SVENY23=DAT_yearly$SVENY23,
                          SVENY24=DAT_yearly$SVENY24,
                          SVENY25=DAT_yearly$SVENY25,
                          SVENY26=DAT_yearly$SVENY26,
                          SVENY27=DAT_yearly$SVENY27,
                          SVENY28=DAT_yearly$SVENY28,
                          SVENY29=DAT_yearly$SVENY29,
                          SVENY30=DAT_yearly$SVENY30)

# Real yield data from GSW 2008:
if(indic.download==1){
  download.file("https://www.federalreserve.gov/data/yield-curve-tables/feds200805.csv",
                "Data/feds200805.csv")
}
DAT <- csv.get(file = "Data/feds200805.csv", skip = 18)
library(tidyverse)
library(zoo)
DAT$Date <- as.Date(DAT$Date)
DAT <- arrange(DAT,Date)
DAT$year <- format(DAT$Date,"%Y")
DAT_yearly <- DAT %>%
  group_by(year) %>%
  summarise_all(function(x){mean(x,na.rm=TRUE)})
days <- as.numeric(format(DAT_yearly$Date,"%d"))
DAT_yearly$Date <- as.Date(paste(format(DAT_yearly$Date,"%Y"),"-01-01",sep=""))
DAT_GSW_real <- data.frame(date=DAT_yearly$Date,
                           TIPSY02=DAT_yearly$TIPSY02,
                           TIPSY03=DAT_yearly$TIPSY03,
                           TIPSY04=DAT_yearly$TIPSY04,
                           TIPSY05=DAT_yearly$TIPSY05,
                           TIPSY06=DAT_yearly$TIPSY06,
                           TIPSY07=DAT_yearly$TIPSY07,
                           TIPSY08=DAT_yearly$TIPSY08,
                           TIPSY09=DAT_yearly$TIPSY09,
                           TIPSY10=DAT_yearly$TIPSY10,
                           TIPSY11=DAT_yearly$TIPSY11,
                           TIPSY12=DAT_yearly$TIPSY12,
                           TIPSY13=DAT_yearly$TIPSY13,
                           TIPSY14=DAT_yearly$TIPSY14,
                           TIPSY15=DAT_yearly$TIPSY15,
                           TIPSY16=DAT_yearly$TIPSY16,
                           TIPSY17=DAT_yearly$TIPSY17,
                           TIPSY18=DAT_yearly$TIPSY18,
                           TIPSY19=DAT_yearly$TIPSY19,
                           TIPSY20=DAT_yearly$TIPSY20)


# Merge all data frames:

# DATA <- merge(DATA,SPF.GDP10,by="date",all=TRUE)
# DATA <- merge(DATA,SPF.GDP,by="date",all=TRUE)
# DATA <- merge(DATA,SPF.CPI10,by="date",all=TRUE)
# DATA <- merge(DATA,SPF.CPI,by="date",all=TRUE)
# DATA <- merge(DATA,SPF.BILL10,by="date",all=TRUE)
DATA <- merge(DATA,DAT_GSW_nom,by="date",all=TRUE)
DATA <- merge(DATA,DAT_GSW_real,by="date",all=TRUE)
# 
# par(mfrow=c(1,2))
# plot(DATA$date,DATA$z,type="l",ylim=c(-.1,.04))
# points(DATA$date,DATA$RGDP1)
# plot(DATA$z,DATA$RGDP1)
# abline(lm(DATA$RGDP1~DATA$z))
# 
# # plot(DATA$date,DATA$z,type="l")
# # lines(DATA$date,DATA$DFII10/100,col="red")
# # lines(DATA$date,DATA$DFII5/100,col="blue")
# 
# # plot(DATA$date,DATA$RGDP1,type="l",ylim=c(-.02,.07))
# # lines(DATA$date,DATA$DFII10/100,col="red")
# # lines(DATA$date,DATA$DFII5/100,col="blue")
# 
# plot(DATA$date,DATA$pi,type="l")
# lines(DATA$date,DATA$CPI1,col="blue")
# lines(DATA$date,DATA$CPI10,col="red")
# 
# plot(DATA$date,DATA$dy,type="l")
# lines(DATA$date,DATA$RGDP1,col="blue")
# points(DATA$date,DATA$RGDP10,col="red")
# 
# DATA$TRP10 <- DATA$THREEFY10 - DATA$BILL10
# 
# DATA$IRP10 <- (DATA$THREEFY10 - DATA$TIPSY10) - DATA$CPI10
# plot(DATA$date,DATA$THREEFYTP10,type="l")
# lines(DATA$date,DATA$IRP10,col="red")
# points(DATA$date,DATA$TRP10)
# 
# save(DATA,file="Data/data_yds.Rda")
# 
# 
# 

save(DATA,file="Data/data.Rda")


