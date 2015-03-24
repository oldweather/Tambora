# Sort Yuri's observations by date, and fill in the gaps between them

library(chron)

# Get Yuri's observations
load('Yuri/0_all_station_and_travel_pressure_data_1815-17.Rdata')
Last.lats<-rep(NA,51) # Use to interpolate positions for fixed stations
Last.lons<-rep(NA,51)
for(station in seq(1,51)) {
  PP[[station]]$chron<-as.chron(PP[[station]]$UTC.date)
  if((max(PP[[station]]$Latitude,na.rm=T)-
     min(PP[[station]]$Latitude,na.rm=T))<0.1) {
    Last.lats[station]<-max(PP[[station]]$Latitude,na.rm=T)
  }
  if((max(PP[[station]]$Longitude,na.rm=T)-
     min(PP[[station]]$Longitude,na.rm=T))<0.1) {
    Last.lons[station]<-max(PP[[station]]$Longitude,na.rm=T)
  }
}
# Make a list of all the dates we are going to use - every 6 hours
#  in 1815-17
c.date<-chron(dates="1815/01/01",times="00:00:00",
          format=c(dates='y/m/d',times='h:m:s'))
e.date<-chron(dates="1817/12/31",times="23:00:00",
          format=c(dates='y/m/d',times='h:m:s'))
Dates = list()
count=1
while(c.date<e.date) {
  for(hour in c(0,6,12,18)) {
     Dates[[count]]<-list(date=as.chron(c.date+hour/24),stations=list())
     for(station in seq(1,51)) {
       Dates[[count]]$stations[[station]]<-list(name=PP[[station]]$Station.name[1],
                                                Longitude=NA,
                                                Latitude=NA,
                                                QFF=NA,
                                                TA.orig=NA)
       w<-which(abs(Dates[[count]]$date-PP[[station]]$chron)<=0.125)
       if(length(w)>0) {
          for(var in c('Longitude','Latitude','QFF','TA.orig')) {
             Dates[[count]]$stations[[station]][[var]]<-PP[[station]][[var]][w[1]]
          }       
       }
       if(is.na(Dates[[count]]$stations[[station]]$Latitude)) {
          Dates[[count]]$stations[[station]]$Latitude<-Last.lats[station]
       }
       if(is.na(Dates[[count]]$stations[[station]]$Longitude)) {
          Dates[[count]]$stations[[station]]$Longitude<-Last.lons[station]
       }
     }       
     count<-count+1
   }
  c.date<-c.date+1
}


save(Dates,file='By.date.Rdata')
