#Calculate station:reanalysis correlations for each station

library(grid)
library(chron)

# Load station data and Reanalysis data
load('../Yuri/0_all_station_and_travel_pressure_data.Rdata')
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
cor.station<-function(station) { # station is an integer in range 1:51

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
       v354$prmsl.mean[i]<-Dates[[i]]$station[[station]]$v354$prmsl.mean/100
       v354$prmsl.spread[i]<-Dates[[i]]$station[[station]]$v354$prmsl.spread/100
       v354$prmsl.normal[i]<-Dates[[i]]$station[[station]]$v354$prmsl.normal/100
       v355$prmsl.mean[i]<-Dates[[i]]$station[[station]]$v355$prmsl.mean/100
       v355$prmsl.spread[i]<-Dates[[i]]$station[[station]]$v355$prmsl.spread/100
       v356$prmsl.mean[i]<-Dates[[i]]$station[[station]]$v356$prmsl.mean/100
       v356$prmsl.spread[i]<-Dates[[i]]$station[[station]]$v356$prmsl.spread/100
    }
    t2m<-strip.diurnal(t2m,dates)
    v354$t2m.mean<-strip.diurnal(v354$t2m.mean,dates)
    v354$t2m.spread<-strip.diurnal(v354$t2m.spread,dates)
    v354$t2m.normal<-strip.diurnal(v354$t2m.normal,dates)
    v355$t2m.mean<-strip.diurnal(v355$t2m.mean,dates)
    v355$t2m.spread<-strip.diurnal(v355$t2m.spread,dates)
    v356$t2m.mean<-strip.diurnal(v356$t2m.mean,dates)
    v356$t2m.spread<-strip.diurnal(v356$t2m.spread,dates)
    # Make anomalies
    t2m<-t2m-v354$t2m.normal
    v354$t2m.mean<-v354$t2m.mean-v354$t2m.normal
    v355$t2m.mean<-v355$t2m.mean-v354$t2m.normal
    v356$t2m.mean<-v356$t2m.mean-v354$t2m.normal
    prmsl<-prmsl-v354$prmsl.normal
    v354$prmsl.mean<-v354$prmsl.mean-v354$prmsl.normal
    v355$prmsl.mean<-v355$prmsl.mean-v354$prmsl.normal
    v356$prmsl.mean<-v356$prmsl.mean-v354$prmsl.normal


    cor.slp.p<-cor(v356$prmsl.mean,prmsl,use='na.or.complete',method='pearson')
    cor.at.p<-cor(v356$t2m.mean,t2m,use='na.or.complete',method='pearson')
    cor.slp.s<-cor(v356$prmsl.mean,prmsl,use='na.or.complete',method='spearman')
    cor.at.s<-cor(v356$t2m.mean,t2m,use='na.or.complete',method='spearman')

    print(sprintf("%-20s %4.2f %4.2f %4.2f %4.2f",PP[[station]]$Station[1],
     cor.slp.p,cor.slp.s,cor.at.p,cor.at.s))
}

for(station in seq(1,51)) {
   cor.station(station)
}
