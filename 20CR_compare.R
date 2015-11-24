# Compare obs and reanalysis data

library(GSDF.TWCR)
library(parallel)
library(chron)

# Get Yuri's observations
load('Yuri/0_all_station_and_travel_pressure_data_1815-17.Rdata')

# Will probably run this more than once, cache the field accesses.
GSDF.cache.dir<-sprintf("%s/GSDF.cache",Sys.getenv('GSCRATCH'))

# Get all the lats and lons for which Yuri has data for a given time
get.ll<-function(year,month,day,hour) {
   indices<-list()
   lats<-numeric()
   lons<-numeric()
   for(i in seq(1,length(PP))) {
      ob.dates<-as.chron(PP[[i]]$UTC.date)
      w<-which(year==as.integer(as.character(years(ob.dates))) &
               month == as.integer(months(ob.dates)) &
               day == days(ob.dates) &
               hour==as.integer(hours(ob.dates)))
      if(length(w)>0) {
         indices[[i]]<-w
         lats<-c(lats,PP[[i]]$Latitude[w])
         lons<-c(lons,PP[[i]]$Longitude[w])
      }
   }
   return(list(indices=indices,lats=lats,lons=lons))
}

# Get means and spreads for each position, for the given time.
get.comparisons<-function(year,month,day,hour,lats,lons,version) {
  
  t2m<-TWCR.get.slice.at.hour('air.2m',year,month,day,hour,
                              version=version)
  t2m.mean<-GSDF.interpolate.ll(t2m,lats,lons)  
  t2m<-TWCR.get.slice.at.hour('air.2m',year,month,day,hour,
                              type='spread',
                              version=version)
  t2m.spread<-GSDF.interpolate.ll(t2m,lats,lons)  
  prmsl<-TWCR.get.slice.at.hour('prmsl',year,month,day,hour,
                              version=version)
  prmsl.mean<-GSDF.interpolate.ll(prmsl,lats,lons)
  prmsl<-TWCR.get.slice.at.hour('prmsl',year,month,day,hour,
                              type='spread',
                              version=version)
  prmsl.spread<-GSDF.interpolate.ll(prmsl,lats,lons)
  return(list(t2m.mean=t2m.mean,t2m.spread=t2m.spread,
              prmsl.mean=prmsl.mean,prmsl.spread=prmsl.spread))
}

#  Get all the necessary reanalysis data for an hour
get.reanalysis<-function(dte) {

  year<-dte$year
  month<-dte$month
  day<-dte$day
  hour<-dte$hour

  li<-get.ll(year,month,day,hour)
  if(length(li$lats)==0) return(NULL)
  v354=get.comparisons(year,month,day,hour,li$lats,li$lons,'3.5.4')
  v355=get.comparisons(year,month,day,hour,li$lats,li$lons,'3.5.5')
  v356=get.comparisons(year,month,day,hour,li$lats,li$lons,'3.5.6')
  return(list(v354=v354,v355=v355,v356=v356,indices=li$indices))
}

# Merge the reanalysis data with the original frame
merge.reanalysis<-function(Reanalysis,res.list) {
   count<-0
   if(!is.list(res.list)) return(Reanalysis)
   if(length(Reanalysis)<length(res.list$indices)) {
    Reanalysis[[length(res.list$indices)]]<-list()
   }
   for(i in seq(1,length(res.list$indices))) {
     if(is.null(Reanalysis[[i]])) Reanalysis[[i]]<-list()
     if(is.null(res.list$indices[[i]])) next
     count<-count+1
     for(ensda in c('v354','v355','v356')) {
        if(is.null(Reanalysis[[i]][[ensda]])) Reanalysis[[i]][[ensda]]<-list()
        for(var in c('t2m.mean','t2m.spread','prmsl.mean','prmsl.spread')) {
           if(is.null(Reanalysis[[i]][[ensda]][[var]])) {
              Reanalysis[[i]][[ensda]][[var]]<-rep(NA,length(PP[[i]]$UTC.date))
           }
           Reanalysis[[i]][[ensda]][[var]][res.list$indices[[i]]]<-res.list[[ensda]][[var]][count]
        }
     }
   }
   return(Reanalysis)
}

c.date<-chron(dates="1817/12/20",times="00:00:00",
          format=c(dates='y/m/d',times='h:m:s'))
e.date<-chron(dates="1817/12/30",times="00:00:00",
          format=c(dates='y/m/d',times='h:m:s'))
Dates = list()
count=1
while(c.date<e.date) {
  for(hour in seq(0,23)) {
     Dates[[count]]<-list(
             year=as.numeric(as.character(years(c.date))),
             month=as.integer(months(c.date)),
             day=as.integer(days(c.date)),
             hour=hour)
     count<-count+1
   }
  c.date<-c.date+1
}

Results<-mclapply(Dates,get.reanalysis,mc.cores=8)

#Reanalysis<-list()
load('Reanalysis.Rdata')
for(i in seq(1:length(Results))) {
   if(is.null(Results[[i]])) next
   Reanalysis<-merge.reanalysis(Reanalysis,Results[[i]])
}
save(Reanalysis,file='Reanalysis.Rdata')
    