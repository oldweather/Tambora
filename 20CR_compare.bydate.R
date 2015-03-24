# Compare obs and reanalysis data

library(GSDF.TWCR)
library(parallel)
library(chron)

# Get data by date
load('By.date.Rdata')

# Will probably run this more than once, cache the field accesses.
GSDF.cache.dir<-sprintf("%s/GSDF.cache",Sys.getenv('GSCRATCH'))

# Get all the lats and lons for which Yuri has data for a given time
get.ll<-function(index) {
   indices<-list()
   lats<-numeric()
   lons<-numeric()
   for(i in seq(1,length(Dates[[index]]$stations))) {
      if(!is.na(Dates[[index]]$stations[[i]]$Latitude) &
         !is.na(Dates[[index]]$stations[[i]]$Longitude)) {
         indices<-c(indices,i)
         lats<-c(lats,Dates[[index]]$stations[[i]]$Latitude)
         lons<-c(lons,Dates[[index]]$stations[[i]]$Longitude)
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
  t2m.normal<-rep(NA,length(t2m.mean))
  if(version=='3.5.4') { 
     t2m<-TWCR.get.slice.at.hour('air.2m',year,month,day,hour,
                                 version='3.5.4',type='normal')
     t2m.normal<-GSDF.interpolate.ll(t2m,lats,lons)
  }  
  prmsl<-TWCR.get.slice.at.hour('prmsl',year,month,day,hour,
                              version=version)
  prmsl.mean<-GSDF.interpolate.ll(prmsl,lats,lons)
  prmsl<-TWCR.get.slice.at.hour('prmsl',year,month,day,hour,
                              type='spread',
                              version=version)
  prmsl.spread<-GSDF.interpolate.ll(prmsl,lats,lons)
  prmsl.normal<-rep(NA,length(prmsl.mean))
  if(version=='3.5.4') { 
    prmsl<-TWCR.get.slice.at.hour('prmsl',year,month,day,hour,
                              version='3.5.4',type='normal')
     prmsl.normal<-GSDF.interpolate.ll(prmsl,lats,lons)
  }
  return(list(t2m.mean=t2m.mean,t2m.spread=t2m.spread,
              t2m.normal=t2m.normal,
              prmsl.mean=prmsl.mean,prmsl.spread=prmsl.spread,
              prmsl.normal=prmsl.normal))
}

#  Get all the necessary reanalysis data for an hour
get.reanalysis<-function(index) {

  dte<-Dates[[index]]$date

  year=as.numeric(as.character(years(dte)))
  month=as.integer(months(dte))
  day=as.integer(days(dte))
  hour=as.integer(hours(dte))

  li<-get.ll(index)
  if(length(li$lats)==0) return(NULL)
  v354=get.comparisons(year,month,day,hour,li$lats,li$lons,'3.5.4')
  v355=get.comparisons(year,month,day,hour,li$lats,li$lons,'3.5.5')
  v356=get.comparisons(year,month,day,hour,li$lats,li$lons,'3.5.6')
  return(list(v354=v354,v355=v355,v356=v356,indices=li$indices))
}

# Merge the reanalysis data with the original frame
merge.reanalysis<-function(Dates,index,res.list) {
   count<-0
   if(!is.list(res.list)) return(Dates)
   for(i in seq(1,51)) {
     for(ensda in c('v354','v355','v356')) {
        if(is.null(Dates[[index]]$stations[[i]][[ensda]])) {
           Dates[[index]]$stations[[i]][[ensda]]<-list()
        }
        for(var in c('t2m.mean','t2m.spread','t2m.normal',
                     'prmsl.mean','prmsl.spread','prmsl.normal')) {
           if(is.null(Dates[[index]]$stations[[i]][[ensda]][[var]])) {
              Dates[[index]]$stations[[i]][[ensda]][[var]]<-NA
           }
        }
      }
   }
   count<-0
   for(i in res.list$indices ) {
      count<-count+1
      for(ensda in c('v354','v355','v356')) {
        for(var in c('t2m.mean','t2m.spread','t2m.normal',
                     'prmsl.mean','prmsl.spread','prmsl.normal')) {
           Dates[[index]]$stations[[i]][[ensda]][[var]]<-res.list[[ensda]][[var]][count]
        }
     }
   }
   return(Dates)
}


Results<-mclapply(seq(1,length(Dates)),get.reanalysis,mc.cores=8)
save(Results,file='Results.Rdata')

for(i in seq(1:length(Results))) {
   print(i)
   if(is.null(Results[[i]])) next
   Dates<-merge.reanalysis(Dates,i,Results[[i]])
   if(i%%100==0) save(Dates,file='Reanalysis.Rdata')
}
save(Dates,file='Reanalysis.Rdata')
    
