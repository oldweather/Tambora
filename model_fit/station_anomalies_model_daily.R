#  Calculate linear model relating obs and model versions

library(grid)
library(chron)
library(pracma)

# Divide stations into high-quality and 'other'
good.stations<-c(2,3,4,7,8,13,16,18,19,20,22,24,25,26,27,28,33,36,46,47,49,50,51)
bad.stations<-c(1,5,6,9,10,11,12,14,15,17,21,23,29,30,31,32,34,35,37,38,39,40,41,42,43,44,45,48)


# Load station data and Reanalysis data
load('../Yuri/0_all_station_and_travel_pressure_data.Rdata')
load('../Reanalysis.364.Rdata')
Dates.364<-Dates
load('../Reanalysis.366.Rdata')
Dates.366<-Dates
load('../Reanalysis.Rdata')


# The model does not have an adequate diurnal cycle - so average each data source over +-12 hours
strip.diurnal<-function(x,x.dates) {
  result<-x
  for(i in seq_along(x)) {
     w<-which(abs(x.dates-x.dates[i])<=0.5)
     if(length(w)>1) result[i]<-mean(x[w],na.rm=TRUE)
  }
  return(result)
}
model.station<-function(station) { # station is an integer in range 1:51

    # Reformat the station data into lists
    dates<-rep(Dates[[1]]$date,length(Dates))
    t2m<-rep(NA,length(Dates))
    prmsl<-rep(NA,length(Dates))
    v354<-list()
    v354$t2m.mean<-rep(NA,length(Dates))
    v354$t2m.spread<-rep(NA,length(Dates))
    v354$t2m.normal<-rep(NA,length(Dates))
    v354$prmsl.mean<-rep(NA,length(Dates))
    v354$prmsl.spread<-rep(NA,length(Dates))
    v354$prmsl.normal<-rep(NA,length(Dates))
    v355<-list()
    v355$t2m.mean<-rep(NA,length(Dates))
    v355$t2m.spread<-rep(NA,length(Dates))
    v355$prmsl.mean<-rep(NA,length(Dates))
    v355$prmsl.spread<-rep(NA,length(Dates))
    v356<-list()
    v356$t2m.mean<-rep(NA,length(Dates))
    v356$t2m.spread<-rep(NA,length(Dates))
    v356$prmsl.mean<-rep(NA,length(Dates))
    v356$prmsl.spread<-rep(NA,length(Dates))
    v364<-list()
    v364$t2m.mean<-rep(NA,length(Dates))
    v364$t2m.spread<-rep(NA,length(Dates))
    v364$prmsl.mean<-rep(NA,length(Dates))
    v364$prmsl.spread<-rep(NA,length(Dates))
    v366<-list()
    v366$t2m.mean<-rep(NA,length(Dates))
    v366$t2m.spread<-rep(NA,length(Dates))
    v366$prmsl.mean<-rep(NA,length(Dates))
    v366$prmsl.spread<-rep(NA,length(Dates))
    for(i in seq_along(Dates)) {
       dates[i]<-Dates[[i]]$date
       t2m[i]<-Dates[[i]]$station[[station]]$TA.orig
       prmsl[i]<-Dates[[i]]$station[[station]]$QFF
       v354$t2m.mean[i]<-Dates[[i]]$station[[station]]$v354$t2m.mean-273.15
       v354$t2m.spread[i]<-Dates[[i]]$station[[station]]$v354$t2m.spread
       v354$t2m.normal[i]<-Dates[[i]]$station[[station]]$v354$t2m.normal-273.15
       v355$t2m.mean[i]<-Dates[[i]]$station[[station]]$v355$t2m.mean-273.15
       v355$t2m.spread[i]<-Dates[[i]]$station[[station]]$v355$t2m.spread
       v356$t2m.mean[i]<-Dates[[i]]$station[[station]]$v356$t2m.mean-273.15
       v356$t2m.spread[i]<-Dates[[i]]$station[[station]]$v356$t2m.spread
       v364$t2m.mean[i]<-Dates.364[[i]]$station[[station]]$v364$t2m.mean-273.15
       v364$t2m.spread[i]<-Dates.364[[i]]$station[[station]]$v364$t2m.spread
       v366$t2m.mean[i]<-Dates.366[[i]]$station[[station]]$v366$t2m.mean-273.15
       v366$t2m.spread[i]<-Dates.366[[i]]$station[[station]]$v366$t2m.spread
       v354$prmsl.mean[i]<-Dates[[i]]$station[[station]]$v354$prmsl.mean/100
       v354$prmsl.spread[i]<-Dates[[i]]$station[[station]]$v354$prmsl.spread/100
       v354$prmsl.normal[i]<-Dates[[i]]$station[[station]]$v354$prmsl.normal/100
       v355$prmsl.mean[i]<-Dates[[i]]$station[[station]]$v355$prmsl.mean/100
       v355$prmsl.spread[i]<-Dates[[i]]$station[[station]]$v355$prmsl.spread/100
       v356$prmsl.mean[i]<-Dates[[i]]$station[[station]]$v356$prmsl.mean/100
       v356$prmsl.spread[i]<-Dates[[i]]$station[[station]]$v356$prmsl.spread/100
       v364$prmsl.mean[i]<-Dates.364[[i]]$station[[station]]$v364$prmsl.mean/100
       v364$prmsl.spread[i]<-Dates.364[[i]]$station[[station]]$v364$prmsl.spread/100
       v366$prmsl.mean[i]<-Dates.366[[i]]$station[[station]]$v366$prmsl.mean/100
       v366$prmsl.spread[i]<-Dates.366[[i]]$station[[station]]$v366$prmsl.spread/100
   }
    t2m<-strip.diurnal(t2m,dates)
    v354$t2m.mean<-strip.diurnal(v354$t2m.mean,dates)
    v354$t2m.spread<-strip.diurnal(v354$t2m.spread,dates)
    v354$t2m.normal<-strip.diurnal(v354$t2m.normal,dates)
    v355$t2m.mean<-strip.diurnal(v355$t2m.mean,dates)
    v355$t2m.spread<-strip.diurnal(v355$t2m.spread,dates)
    v356$t2m.mean<-strip.diurnal(v356$t2m.mean,dates)
    v356$t2m.spread<-strip.diurnal(v356$t2m.spread,dates)
    v364$t2m.mean<-strip.diurnal(v364$t2m.mean,dates)
    v364$t2m.spread<-strip.diurnal(v364$t2m.spread,dates)
    v366$t2m.mean<-strip.diurnal(v366$t2m.mean,dates)
    v366$t2m.spread<-strip.diurnal(v366$t2m.spread,dates)
    # Make anomalies
    t2m<-t2m-v354$t2m.normal
    v354$t2m.mean<-v354$t2m.mean-v354$t2m.normal
    v355$t2m.mean<-v355$t2m.mean-v354$t2m.normal
    v356$t2m.mean<-v356$t2m.mean-v354$t2m.normal
    v364$t2m.mean<-v364$t2m.mean-v354$t2m.normal
    v366$t2m.mean<-v366$t2m.mean-v354$t2m.normal
    prmsl<-prmsl-v354$prmsl.normal
    v354$prmsl.mean<-v354$prmsl.mean-v354$prmsl.normal
    v355$prmsl.mean<-v355$prmsl.mean-v354$prmsl.normal
    v356$prmsl.mean<-v356$prmsl.mean-v354$prmsl.normal
    v364$prmsl.mean<-v364$prmsl.mean-v354$prmsl.normal
    v366$prmsl.mean<-v366$prmsl.mean-v354$prmsl.normal
  
    aerosol.diff<-v356$t2m.mean-v354$t2m.mean
    sst.diff<-v364$t2m.mean-v354$t2m.mean
    forced.diff<-v366$t2m.mean-v354$t2m.mean

    # fit models by orthogonal regression
    w<-which(!is.na(t2m))
    lm.0f<-odregress(v354$t2m.mean[w],t2m[w])
    w<-which(!is.na(t2m))
    lm.cf<-odregress(v366$t2m.mean[w],t2m[w])
    w<-which(!is.na(t2m)&!is.na(forced.diff))
    lm.1f<-odregress(cbind(v354$t2m.mean[w],forced.diff[w]),t2m[w])
    w<-which(!is.na(t2m)&!is.na(sst.diff))
    lm.1s<-odregress(cbind(v354$t2m.mean[w],sst.diff[w]),t2m[w])
    w<-which(!is.na(t2m)&!is.na(aerosol.diff))
    lm.1a<-odregress(cbind(v354$t2m.mean[w],aerosol.diff[w]),t2m[w])
    w<-which(!is.na(t2m)&!is.na(aerosol.diff)&!is.na(sst.diff))
    lm.2f<-odregress(cbind(v354$t2m.mean[w],aerosol.diff[w],sst.diff[w]),t2m[w])

    return(list(lm.0f=lm.0f,lm.1f=lm.1f,lm.2f=lm.2f,lm.1s=lm.1s,lm.1a=lm.1a,lm.cf=lm.cf,
                dates=dates,obs=t2m,
                v354=v354$t2m.mean,v366=v366$t2m.mean,spread<-v366$t2m.spread,
                aerosol.diff=aerosol.diff,sst.diff=sst.diff,forced.diff=forced.diff))
}

models<-list()
for(station in good.stations) {
   models[[station]]<-model.station(station)
}

save(models,file='models.Rdata')

