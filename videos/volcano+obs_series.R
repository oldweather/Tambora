#!/usr/common/graphics/R/R-3.1.1/bin/Rscript --no-save

# Pressure fields from 20CR2c for 1815/7
#  Add the temperature anomalies from Yri's data
# Show the time-series beside the plot

library(GSDF.TWCR)
library(GSDF.WeatherMap)
library(parallel)
library(getopt)

opt = getopt(c(
  'year',   'y', 2, "integer",
  'month',   'm', 2, "integer",
  'version',   'v', 2, "character",
  'nmonths','n', 2, "integer",
  'obs', 'o', 2, 'character'
));
if ( is.null(opt$year  ) )   { stop("Year not specified") }
if ( is.null(opt$month  ) )   { stop("Month not specified") }
if ( is.null(opt$nmonths ) ) { opt$nmonths  = 6 }
if ( is.null(opt$version ) ) { opt$version  = '3.5.4' }
if ( is.null(opt$obs ) ) { opt$obs='F' }

Year<-opt$year
Month<-opt$month
Day<-1
Hour<-0
n.total<-opt$nmonths*31*24 # Total number of hours to be rendered
version<-opt$version
fog.threshold<-exp(1)

GSDF.cache.dir<-sprintf("%s/GSDF.cache",Sys.getenv('GSCRATCH'))
Imagedir<-sprintf("%s/images/Tambora+obs.%s",Sys.getenv('GSCRATCH'),opt$version)
if(!file.exists(Imagedir)) dir.create(Imagedir,recursive=TRUE)

c.date<-chron(dates=sprintf("%04d/%02d/%02d",Year,Month,Day),
          times=sprintf("%02d:00:00",Hour),
          format=c(dates='y/m/d',times='h:m:s'))

Options<-WeatherMap.set.option(NULL)

Options<-WeatherMap.set.option(Options,'lat.min',-40)
Options<-WeatherMap.set.option(Options,'lat.max',40)
Options<-WeatherMap.set.option(Options,'lon.min',-35)
Options<-WeatherMap.set.option(Options,'lon.max',105)
Options<-WeatherMap.set.option(Options,'pole.lon',10)
Options<-WeatherMap.set.option(Options,'pole.lat',40)
Options$label.yp<-0.03
Options$mslp.base=0                    # Base value for anomalies
#Options$mslp.base=101325                    # Base value for anomalies
Options$mslp.range=50000                    # Anomaly for max contour
Options$mslp.step=500                       # Smaller -> more contours
Options$mslp.tpscale=500                    # Smaller -> contours less transparent
Options$mslp.lwd=2
Options<-WeatherMap.set.option(Options,'obs.size',0.5)
Options<-WeatherMap.set.option(Options,'obs.colour',rgb(255,215,0,255,
                                                       maxColorValue=255))

Options$ice.points<-50000
land<-WeatherMap.get.land(Options)

load('../Reanalysis.Rdata') # obs data and comparisons
dates<-rep(Dates[[1]]$date,length(Dates))
for(i in seq_along(Dates)) {
   dates[i]<-Dates[[i]]$date
} # Dates vector for selection
T2M<-list()
for(station in seq(1,51)) {
   T2M[[station]]<-rep(NA,length(Dates))
   for(i in seq_along(Dates)) {
      T2M[[station]][i]<-Dates[[i]]$station[[station]]$TA.orig
   }
}

Draw.obs<-function(year,month,day,hour,Options,Trange=10,size=0.015) {
   od<-chron(dates=sprintf("%04d/%02d/%02d",year,month,day),
          times=sprintf("%02d:00:00",hour),
          format=c(dates='y/m/d',times='h:m:s'))
   for(station in c(18,24,36,47,51)) {
      w<-which(abs(as.numeric(dates)-as.numeric(od))<0.125 &
               !is.na(T2M[[station]]))
      if(length(w)<1) {
         w<-which(abs(as.numeric(dates)-as.numeric(od))<0.25 &
                  !is.na(T2M[[station]]))
      }
      if(length(w)<1) {
         w<-which(abs(as.numeric(dates)-as.numeric(od))<0.5 &
                  !is.na(T2M[[station]]))
      }
      if(length(w)==0) next
      t2m<-Dates[[w[1]]]$station[[station]]$TA.orig
      Lat<-Dates[[w[1]]]$station[[station]]$Latitude
      Lon<-Dates[[w[1]]]$station[[station]]$Longitude
      if(is.na(t2m) | is.na(Lat) | is.na(Lon)) next
      if(Options$pole.lon!=0 || Options$pole.lat!=90) {
           l2<-GSDF.ll.to.rg(Lat,Lon,Options$pole.lat,Options$pole.lon)
           Lon<-l2$lon
           Lat<-l2$lat
      }
      t2m<-t2m-(Dates[[w[1]]]$station[[station]]$v354$t2m.normal-273.15)
      spread<-Dates[[w[1]]]$station[[station]]$v354$t2m.spread
      if(abs(t2m/spread)>3) col=rgb(1,1,1,1)
      else col=rgb(0,0,0,1)
      if(t2m>0) fill=rgb(1,0,0,min(0.99,t2m/Trange))
      else fill=rgb(0,0,1,min(0.99,-1*t2m/Trange))
      gp<-gpar(col=col,fill=fill,lwd=2)
      grid.points(x=unit(Lon,'native'),
                  y=unit(Lat,'native'),
                  pch=21,
                  size=unit(size,'snpc'),
                  gp=gp)
   }
}


Draw.temperature<-function(temperature,Options,Trange=10) {
  
  Options.local<-Options
  Options.local$fog.min.transparency<-0.5
  tplus<-temperature
  tplus$data[]<-pmax(1,pmin(Trange,tplus$data))/Trange
  Options.local$fog.colour<-c(1,0,0)
  WeatherMap.draw.fog(tplus,Options.local)
  tminus<-temperature
  tminus$data[]<-tminus$data*-1
  tminus$data[]<-pmax(1,pmin(Trange,tminus$data))/Trange
  Options.local$fog.colour<-c(0,0,1)
  WeatherMap.draw.fog(tminus,Options.local)
}



Draw.pressure<-function(mslp,Options,colour=c(0,0,0)) {

  M<-GSDF.WeatherMap:::WeatherMap.rotate.pole(mslp,Options)
  lats<-M$dimensions[[GSDF.find.dimension(M,'lat')]]$values
  longs<-M$dimensions[[GSDF.find.dimension(M,'lon')]]$values
    # Need particular data format for contourLines
  if(lats[2]<lats[1] || longs[2]<longs[1] || max(longs) > 180 ) {
    if(lats[2]<lats[1]) lats<-rev(lats)
    if(longs[2]<longs[1]) longs<-rev(longs)
    longs[longs>180]<-longs[longs>180]-360
    longs<-sort(longs)
    M2<-M
    M2$dimensions[[GSDF.find.dimension(M,'lat')]]$values<-lats
    M2$dimensions[[GSDF.find.dimension(M,'lon')]]$values<-longs
    M<-GSDF.regrid.2d(M,M2)
  }
  z<-matrix(data=M$data,nrow=length(longs),ncol=length(lats))
  contour.levels<-seq(Options$mslp.base-Options$mslp.range,
                      Options$mslp.base+Options$mslp.range,
                      Options$mslp.step)
  lines<-contourLines(longs,lats,z,
                       levels=contour.levels)
  if(!is.na(lines) && length(lines)>0) {
     for(i in seq(1,length(lines))) {
         tp<-min(1,(abs(lines[[i]]$level-Options$mslp.base)/
                    Options$mslp.tpscale))
         lt<-2
         lwd<-1
         if(lines[[i]]$level<=Options$mslp.base) {
             lt<-1
             lwd<-1
         }
         gp<-gpar(col=rgb(colour[1],colour[2],colour[3],tp),
                             lwd=Options$mslp.lwd*lwd,lty=lt)
         grid.xspline(x=unit(lines[[i]]$x,'native'),
                    y=unit(lines[[i]]$y,'native'),
                    shape=1,
                    gp=gp)
     }
  }
}

plot.station<-function(station) { # station is an integer in range 1:51

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
    # Make anomalies
    t2m<-t2m-v354$t2m.normal
    v354$t2m.mean<-v354$t2m.mean-v354$t2m.normal
    v355$t2m.mean<-v355$t2m.mean-v354$t2m.normal
    v356$t2m.mean<-v356$t2m.mean-v354$t2m.normal
    prmsl<-prmsl-v354$prmsl.normal
    v354$prmsl.mean<-v354$prmsl.mean-v354$prmsl.normal
    v355$prmsl.mean<-v355$prmsl.mean-v354$prmsl.normal
    v356$prmsl.mean<-v356$prmsl.mean-v354$prmsl.normal

   date.range<-chron(dates=c("1815/01/01","1817/12/31"),
                        times=c("00:0:01","23:59:59"),
                        format=c(dates = "y/m/d", times = "h:m:s"))
    tics=pretty(date.range)
    ticl=attr(tics,'labels')

    gpe<-list()
    gpe$v354<-gpar(col=rgb(0.8,0.8,1,1),fill=rgb(0.8,0.8,1,1))
    gpe$v355<-gpar(col=rgb(0.6,0.6,1,1),fill=rgb(0.6,0.6,1,1))
    gpe$v356<-gpar(col=rgb(0.4,0.4,1,1),fill=rgb(0.4,0.4,1,1))

    # Pressure along the bottom with x axis
    pmin<- -20
    for(ensda in c('v354','v355','v356')) {
        v<-get(ensda)
        rmin<-min(v$prmsl.mean-v$prmsl.spread*2,na.rm=TRUE)
        pmin<-min(pmin,rmin,na.rm=TRUE)
    }
    omin<-quantile(prmsl,probs=0.05,na.rm=TRUE)
    pmin<-min(pmin,omin,na.rm=TRUE)
    pmax<-10
    for(ensda in c('v354','v355','v356')) {
        v<-get(ensda)
        rmax<-max(v$prmsl.mean+v$prmsl.spread*2,na.rm=TRUE)
        pmax<-max(pmax,rmax,na.rm=TRUE)
    }
   omax<-quantile(prmsl,probs=0.95,na.rm=TRUE)
    pmax<-max(pmax,omax,na.rm=TRUE)
    pushViewport(viewport(width=1.0,height=0.55,x=0.0,y=0.0,
                          just=c("left","bottom"),name="Page",clip='off'))
       pushViewport(plotViewport(margins=c(4,6,0,0)))
          pushViewport(dataViewport(date.range,c(pmin,pmax),clip='off'))
             grid.xaxis(at=as.numeric(tics),label=ticl,main=T)
             grid.text('Date',y=unit(-3,'lines'))
             grid.yaxis(main=T)
             grid.text('Sea-level pressure (hPa)',x=unit(-4,'lines'),rot=90)

          # Mark the eruption
             gp=gpar(col=rgb(0,0,0,1),fill=rgb(0,0,0,1))

             p.y<-c(pmin,2040)
             p.x<-chron(dates=c("1815/04/10","1815/04/10"),
                        times=rep("12:00:00",2),
                        format=c(dates = "y/m/d", times = "h:m:s"))
             grid.lines(x=unit(p.x,'native'),
                          y=unit(p.y,'native'),
                          gp=gp)


             # Analysis spreads

             for(ensda in c('v354','v356','v355')) { # Order matters
                 v<-get(ensda)
                 for(i in seq_along(dates)) {
                    x<-c(dates[i]-0.125,
                         dates[i]+0.125,
                         dates[i]+0.125,
                         dates[i]-0.125)
                    y<-c(v$prmsl.mean[i]-
                           v$prmsl.spread[i]*2,
                         v$prmsl.mean[i]-
                          v$prmsl.spread[i]*2,
                         v$prmsl.mean[i]+
                          v$prmsl.spread[i]*2,
                         v$prmsl.mean[i]+
                          v$prmsl.spread[i]*2)
                    grid.polygon(x=unit(x,'native'),
                                 y=unit(y,'native'),
                              gp=gpe[[ensda]])
                  }
             }
           # Observation
             gp=gpar(col=rgb(0,0,0,1),fill=rgb(0,0,0,1))
             grid.points(x=unit(dates,'native'),
                         y=unit(prmsl,'native'),
                         size=unit(0.005,'npc'),
                         pch=20,
                         gp=gp)
          popViewport()
       popViewport()
    popViewport()

    # AT on top
    tmin<- -5
    for(ensda in c('v354','v355','v356')) {
        v<-get(ensda)
        rmin<-min(v$t2m.mean-v$t2m.spread*2,na.rm=TRUE)
        tmin<-min(tmin,rmin,na.rm=TRUE)
    }
    omin<-quantile(t2m,probs=0.05,na.rm=TRUE)
    tmin<-min(tmin,omin,na.rm=TRUE)
    tmax<- 5
    for(ensda in c('v354','v355','v356')) {
        v<-get(ensda)
        rmax<-max(v$t2m.mean+v$t2m.spread*2,na.rm=TRUE)
        tmax<-max(tmax,rmax,na.rm=TRUE)
    }
    omax<-quantile(t2m,probs=0.95,na.rm=TRUE)
    tmax<-max(tmax,omax,na.rm=TRUE)
    pushViewport(viewport(width=1.0,height=0.45,x=0.0,y=0.55,
                          just=c("left","bottom"),name="Page",clip='off'))
       pushViewport(plotViewport(margins=c(0,6,0,0)))
          pushViewport(dataViewport(date.range,c(tmin,tmax),clip='off'))
             grid.yaxis(main=T)
             grid.text('Air temperature (C)',x=unit(-4,'lines'),rot=90)

          # Mark the eruption
             gp=gpar(col=rgb(0,0,0,1),fill=rgb(0,0,0,1))

             p.y<-c(-200,300)
             p.x<-chron(dates=c("1815/04/10","1815/04/10"),
                        times=rep("12:00:00",2),
                        format=c(dates = "y/m/d", times = "h:m:s"))
             grid.lines(x=unit(p.x,'native'),
                          y=unit(p.y,'native'),
                          gp=gp)

             grid.yaxis(main=T)
             grid.text('Air temperature (C)',x=unit(-4,'lines'),rot=90)

             # Analysis spreads

             for(ensda in c('v354','v356','v355')) { # Order matters
                 v<-get(ensda)
                 for(i in seq_along(dates)) {
                    x<-c(dates[i]-0.125,
                         dates[i]+0.125,
                         dates[i]+0.125,
                         dates[i]-0.125)
                    y<-c(v$t2m.mean[i]-v$t2m.spread[i]*2,
                         v$t2m.mean[i]-v$t2m.spread[i]*2,
                         v$t2m.mean[i]+v$t2m.spread[i]*2,
                         v$t2m.mean[i]+v$t2m.spread[i]*2)
                    grid.polygon(x=unit(x,'native'),
                                 y=unit(y,'native'),
                              gp=gpe[[ensda]])
                  }
              }
            # Observation
             gp=gpar(col=rgb(0,0,0,1),fill=rgb(0,0,0,1))
             grid.points(x=unit(dates,'native'),
                         y=unit(t2m,'native'),
                         size=unit(0.005,'npc'),
                         pch=20,
                         gp=gp)
          popViewport()
       popViewport()
    popViewport()

    # Name the station
    grid.text(PP[[station]]$Station.name[1],
              x=unit(0.03,'npc'),y=unit(0.03,'npc'),just<-c('left','bottom'))

}


plot.hour<-function(year,month,day,hour) {    

    image.name<-sprintf("%04d-%02d-%02d:%02d.png",year,month,day,hour)

    ifile.name<-sprintf("%s/%s",Imagedir,image.name)
    if(file.exists(ifile.name) && file.info(ifile.name)$size>0) return()

    prmsl<-TWCR.get.slice.at.hour('prmsl',year,month,day,hour,version=version)
    prmsl.spread<-TWCR.get.slice.at.hour('prmsl',year,month,day,hour,version=version,
                                              type='spread')
    prmsl.sd<-TWCR.get.slice.at.hour('prmsl',year,month,day,hour,
                                         version='3.4.1',type='standard.deviation')
    prmsl.normal<-TWCR.get.slice.at.hour('prmsl',year,month,day,hour,version='3.5.4',
                                             type='normal')
    fog<-TWCR.relative.entropy(prmsl.normal,prmsl.sd,prmsl,prmsl.spread)
    fog$data[]<-1-pmin(fog.threshold,pmax(0,fog$data))/fog.threshold
    if(opt$obs=='T') {
       obs<-TWCR.get.obs(year,month,day,hour,version=version)
       w<-which(obs$Longitude>180)
       obs$Longitude[w]<-obs$Longitude[w]-360
    }
    t2m<-TWCR.get.slice.at.hour('air.2m',year,month,day,hour,version=version)
    t2m.normal<-TWCR.get.slice.at.hour('air.2m',year,month,day,hour,version='3.5.4',
                                             type='normal')
    prmsl$data[]<-prmsl$data-prmsl.normal$data
    t2m$data[]<-t2m$data-t2m.normal$data

     png(ifile.name,
             width=1080*16/9,
             height=1080,
             bg=Options$sea.colour,
             pointsize=24,
             type='cairo')
       pushViewport(dataViewport(c(Options$lon.min,Options$lon.max),
                                 c(Options$lat.min,Options$lat.max),
                                  extension=0))
          ip<-WeatherMap.rectpoints(Options$ice.points,Options)
          #WeatherMap.draw.ice(ip$lat,ip$lon,icec,Options)
          WeatherMap.draw.land(land,Options)

        Draw.temperature(t2m,Options)
        Draw.pressure(prmsl,Options,colour=c(0,0,0))
        WeatherMap.draw.fog(fog,Options)

        Draw.obs(year,month,day,hour,Options) 
        if(opt$obs=='T') WeatherMap.draw.obs(obs,Options)

        Options$label<-sprintf("%04d-%02d-%02d:%02d",year,month,day,hour)
        WeatherMap.draw.label(Options)

   # Add the data plots
    pushViewport(viewport(width=0.5,height=1.0,x=0.5,y=0.0,
                          just=c("left","bottom"),name="Page",clip='off'))

    # Plain background
    grid.polygon(x=unit(c(0,1,1,0),'npc'),y=unit(c(0,0,1,1),'npc'),
                 gp=gpar(col='white',fill='white'))


     upViewport()

    dev.off()
}

for(n.count in seq(0,n.total)) {

    n.date<-c.date+n.count/24
    year<-as.numeric(as.character(years(n.date)))
    month<-months(n.date)
    day<-days(n.date)
    hour<-n.count%%24

    image.name<-sprintf("%04d-%02d-%02d:%02d.png",year,month,day,hour)
    ifile.name<-sprintf("%s/%s",Imagedir,image.name)
    if(file.exists(ifile.name) && file.info(ifile.name)$size>0) next
    # Each plot in a seperate parallel process
    mcparallel(plot.hour(year,month,day,hour))
    if(hour==20) mccollect(wait=TRUE)

}
mccollect()
