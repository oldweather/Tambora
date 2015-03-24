#!/usr/common/graphics/R/R-3.1.1/bin/Rscript --no-save

# Pressure fields from 20CR2c for 1815/6

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
if ( is.null(opt$obs ) ) { opt$obs='T' }

Year<-opt$year
Month<-opt$month
Day<-1
Hour<-0
n.total<-opt$nmonths*31*24 # Total number of hours to be rendered
version<-opt$version
fog.threshold<-exp(1)

GSDF.cache.dir<-sprintf("%s/GSDF.cache",Sys.getenv('GSCRATCH'))
Imagedir<-sprintf("%s/images/Tambora.%s",Sys.getenv('GSCRATCH'),opt$version)
if(!file.exists(Imagedir)) dir.create(Imagedir,recursive=TRUE)

c.date<-chron(dates=sprintf("%04d/%02d/%02d",Year,Month,Day),
          times=sprintf("%02d:00:00",Hour),
          format=c(dates='y/m/d',times='h:m:s'))

Options<-WeatherMap.set.option(NULL)

Options<-WeatherMap.set.option(Options,'lat.min',-40)
Options<-WeatherMap.set.option(Options,'lat.max',40)
Options<-WeatherMap.set.option(Options,'lon.min',-70)
Options<-WeatherMap.set.option(Options,'lon.max',70)
Options<-WeatherMap.set.option(Options,'pole.lon',150)
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
        if(opt$obs=='T') WeatherMap.draw.obs(obs,Options)
        Draw.pressure(prmsl,Options,colour=c(0,0,0))
        WeatherMap.draw.fog(fog,Options)

        Options$label<-sprintf("%04d-%02d-%02d:%02d",year,month,day,hour)
        WeatherMap.draw.label(Options)
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
