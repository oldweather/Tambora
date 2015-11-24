# Arrange station temperature anomalies by date for adding to the reanalysis videos

library(grid)
library(chron)

# Useful subset of stations for temperature comparisons
good.stations<-c(2,3,4,7,8,13,16,18,19,20,22,24,25,26,27,28,33,36,46,47,49,50,51)

# Load station data and Reanalysis data
load('../Yuri/0_all_station_and_travel_pressure_data.Rdata')
load('../Reanalysis.Rdata')
Dates.356<-Dates
load('../Reanalysis.366.Rdata')
Dates.366<-Dates
load('../Reanalysis.364.Rdata')
Dates.364<-Dates

# The model does not have an adequate diurnal cycle - so average each data source over +-12 hours
strip.diurnal<-function(x,x.dates) {
  result<-x
  for(i in seq_along(x)) {
     w<-which(abs(x.dates-x.dates[i])<=0.5)
     if(length(w)>1) result[i]<-mean(x[w],na.rm=TRUE)
  }
  return(result)
}

Station.anomalies<-list()
Station.anomalies$t2m<-list()
Station.anomalies$longitude<-list()
Station.anomalies$latitude<-list()
Station.anomalies$elevation<-list()

normalise.station<-function(station,Station.anomalies) { # station is an integer in range 1:51

    # Store the fixed data
    Station.anomalies$longitude[[station]]<-PP[[station]]$Longitude[1]
    Station.anomalies$latitude[[station]]<-PP[[station]]$Latitude[1]
    Station.anomalies$elevation[[station]]<-PP[[station]]$Elevation[1]

    # Reformat the station data into lists
    dates<-rep(Dates[[1]]$date,length(Dates))
    t2m<-rep(NA,length(Dates))
    v364<-list()
    v364$t2m.mean<-rep(NA,length(Dates))
    v364$t2m.spread<-rep(NA,length(Dates))
    v366<-list()
    v366$t2m.mean<-rep(NA,length(Dates))
    v366$t2m.spread<-rep(NA,length(Dates))
    v354<-list()
    v354$t2m.mean<-rep(NA,length(Dates))
    v354$t2m.spread<-rep(NA,length(Dates))
    v354$t2m.normal<-rep(NA,length(Dates))
    v356<-list()
    v356$t2m.mean<-rep(NA,length(Dates))
    v356$t2m.spread<-rep(NA,length(Dates))
    for(i in seq_along(Dates)) {
       dates[i]<-Dates[[i]]$date
       t2m[i]<-Dates[[i]]$station[[station]]$TA.orig
       v364$t2m.mean[i]<-Dates.364[[i]]$station[[station]]$v364$t2m.mean-273.15
       v364$t2m.spread[i]<-Dates.364[[i]]$station[[station]]$v364$t2m.spread
       v366$t2m.mean[i]<-Dates.366[[i]]$station[[station]]$v366$t2m.mean-273.15
       v366$t2m.spread[i]<-Dates.366[[i]]$station[[station]]$v366$t2m.spread
       v354$t2m.mean[i]<-Dates.356[[i]]$station[[station]]$v354$t2m.mean-273.15
       v354$t2m.spread[i]<-Dates.356[[i]]$station[[station]]$v354$t2m.spread
       v354$t2m.normal[i]<-Dates.356[[i]]$station[[station]]$v354$t2m.normal-273.15
       v356$t2m.mean[i]<-Dates.356[[i]]$station[[station]]$v356$t2m.mean-273.15
       v356$t2m.spread[i]<-Dates.356[[i]]$station[[station]]$v356$t2m.spread
    }
    t2m<-strip.diurnal(t2m,dates)
    v364$t2m.mean<-strip.diurnal(v364$t2m.mean,dates)
    v364$t2m.spread<-strip.diurnal(v364$t2m.spread,dates)
    v366$t2m.mean<-strip.diurnal(v366$t2m.mean,dates)
    v366$t2m.spread<-strip.diurnal(v366$t2m.spread,dates)
    v354$t2m.mean<-strip.diurnal(v354$t2m.mean,dates)
    v354$t2m.spread<-strip.diurnal(v354$t2m.spread,dates)
    v354$t2m.normal<-strip.diurnal(v354$t2m.normal,dates)
    v356$t2m.mean<-strip.diurnal(v356$t2m.mean,dates)
    v356$t2m.spread<-strip.diurnal(v356$t2m.spread,dates)
   # Make anomalies
    t2m<-t2m-v354$t2m.normal
    v354$t2m.mean<-v354$t2m.mean-v354$t2m.normal
    v356$t2m.mean<-v356$t2m.mean-v354$t2m.normal
    v366$t2m.mean<-v366$t2m.mean-v354$t2m.normal
    v364$t2m.mean<-v364$t2m.mean-v354$t2m.normal
    # Match median anomaly of obs to 366
    t2m<-t2m+median(v366$t2m.mean,na.rm=TRUE)-median(t2m,na.rm=TRUE)

    # Now put the temperature anomalies into the right format for adding to the videos.
    Station.anomalies$dates<-dates
    Station.anomalies$t2m[[station]]<-list()
    Station.anomalies$t2m[[station]]$obs<-t2m
    Station.anomalies$t2m[[station]]$v354.mean<-v354$t2m.mean
    Station.anomalies$t2m[[station]]$v354.spread<-v354$t2m.mean
    Station.anomalies$t2m[[station]]$v356.mean<-v356$t2m.mean
    Station.anomalies$t2m[[station]]$v356.spread<-v356$t2m.mean
    Station.anomalies$t2m[[station]]$v364.mean<-v364$t2m.mean
    Station.anomalies$t2m[[station]]$v364.spread<-v364$t2m.mean
    Station.anomalies$t2m[[station]]$v366.mean<-v354$t2m.mean
    Station.anomalies$t2m[[station]]$v366.spread<-v354$t2m.mean
   
    return(Station.anomalies)

}

for(station in good.stations) { # 51
   Station.anomalies<-normalise.station(station,Station.anomalies)
}

save(Station.anomalies,file='Station.anomalies.Rdata')