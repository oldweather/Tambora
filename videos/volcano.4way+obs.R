#!/usr/common/graphics/R/R-3.1.1/bin/Rscript --no-save

# Pressures and temperatures from 20CR2c for 1815-17
# Show 24-hour running mean temperatures (avoids diurnal cycle problem).

# Compares versions 3.5.4 3.5.6 3.6.4 and 3.6.6 (Europe only)

library(GSDF.TWCR)
library(GSDF.WeatherMap)
library(parallel)
library(getopt)
library(chron)

opt = getopt(c(
  'year',   'y', 2, "integer",
  'month',   'm', 2, "integer",
  'nmonths','n', 2, "integer",
  'obs', 'o', 2, 'character'
));
if ( is.null(opt$year  ) )   { stop("Year not specified") }
if ( is.null(opt$month  ) )   { stop("Month not specified") }
if ( is.null(opt$nmonths ) ) { opt$nmonths  = 6 }
if ( is.null(opt$obs ) ) { opt$obs='T' }

# Useful subset of stations for temperature comparisons
good.stations<-c(2,3,4,7,8,13,16,18,19,20,22,24,25,26,27,28,33,36,46,47,49,50,51)
#load('../Yuri/0_all_station_and_travel_pressure_data.Rdata')
load('Station.anomalies.Rdata')

Year<-opt$year
Month<-opt$month
Day<-1
Hour<-0
n.total<-opt$nmonths*31*24 # Total number of hours to be rendered
fog.threshold<-exp(1)

GSDF.cache.dir<-sprintf("%s/GSDF.cache",Sys.getenv('SCRATCH'))
if(!file.exists(GSDF.cache.dir)) dir.create(GSDF.cache.dir,recursive=TRUE)
Imagedir<-sprintf("%s/images/Tambora.4way+obs",Sys.getenv('SCRATCH'))
if(!file.exists(Imagedir)) dir.create(Imagedir,recursive=TRUE)

c.date<-chron(dates=sprintf("%04d/%02d/%02d",Year,Month,Day),
          times=sprintf("%02d:00:00",Hour),
          format=c(dates='y/m/d',times='h:m:s'))

Options<-WeatherMap.set.option(NULL)

Options<-WeatherMap.set.option(Options,'lat.min',-20)
Options<-WeatherMap.set.option(Options,'lat.max',20)
Options<-WeatherMap.set.option(Options,'lon.min',-35)
Options<-WeatherMap.set.option(Options,'lon.max',35)
Options<-WeatherMap.set.option(Options,'pole.lon',195)
Options<-WeatherMap.set.option(Options,'pole.lat',35)
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

Draw.obs<-function(year,month,day,hour,Options,Trange=10,size=0.015,version=version) {
   od<-chron(dates=sprintf("%04d/%02d/%02d",year,month,day),
          times=sprintf("%02d:00:00",hour),
          format=c(dates='y/m/d',times='h:m:s'))
   for(station in good.stations) {
      w<-which(abs(as.numeric(Station.anomalies$dates)-as.numeric(od))<0.125 &
               !is.na(Station.anomalies$t2m[[station]]$obs))
      if(length(w)<1) {
         w<-which(abs(as.numeric(Station.anomalies$dates)-as.numeric(od))<0.25 &
                  !is.na(Station.anomalies$t2m[[station]]$obs))
      }
      if(length(w)<1) {
         w<-which(abs(as.numeric(Station.anomalies$dates)-as.numeric(od))<0.5 &
                  !is.na(Station.anomalies$t2m[[station]]$obs))
      }
      if(length(w)==0) next
      t2m<-Station.anomalies$t2m[[station]]$obs[w[1]]
      Lat<-Station.anomalies$latitude[[station]]
      Lon<-Station.anomalies$longitude[[station]]
      if(is.na(t2m) | is.na(Lat) | is.na(Lon)) next
      if(Options$pole.lon!=0 || Options$pole.lat!=90) {
           l2<-GSDF.ll.to.rg(Lat,Lon,Options$pole.lat,Options$pole.lon)
           Lon<-l2$lon
           Lat<-l2$lat
      }
      spread<-Station.anomalies$t2m[[station]]$v354.spread[w[1]]
      if(version=='3.5.6') spread<-Station.anomalies$t2m[[station]]$v356.spread[w[1]]
      if(version=='3.6.4') spread<-Station.anomalies$t2m[[station]]$v364.spread[w[1]]
      if(version=='3.6.6') spread<-Station.anomalies$t2m[[station]]$v366.spread[w[1]]
      t.diff<-t2m-Station.anomalies$t2m[[station]]$v354.mean[w[1]]
      if(version=='3.5.6') t.diff<-t2m-Station.anomalies$t2m[[station]]$v356.mean[w[1]]
      if(version=='3.6.4') t.diff<-t2m-Station.anomalies$t2m[[station]]$v364.mean[w[1]]
      if(version=='3.6.6') t.diff<-t2m-Station.anomalies$t2m[[station]]$v366.mean[w[1]]
      #if(abs(t.diff/spread)>3) col=rgb(1,1,1,1)
      col=rgb(0,0,0,1)
      if(t2m>0) fill=rgb(1,0,0,min(0.99,t2m/Trange))
      else fill=rgb(0,0,1,min(0.99,-1*t2m/Trange))
      gp<-gpar(col=col,fill=fill,lwd=1)
      grid.points(x=unit(Lon,'native'),
                  y=unit(Lat,'native'),
                  pch=21,
                  size=unit(size,'snpc'),
                  gp=gp)
   }
}



Options$ice.points<-50000
land<-WeatherMap.get.land(Options)

next.day<-function(year,month,day) {
  c.date<-chron(dates=sprintf("%04d/%02d/%02d",year,month,day),
          times=sprintf("%02d:00:00",12),
          format=c(dates='y/m/d',times='h:m:s'))+1
  return(c(as.integer(as.character(years(c.date))),
           as.integer(months(c.date)),
           as.integer(days(c.date))))
}
previous.day<-function(year,month,day) {
  c.date<-chron(dates=sprintf("%04d/%02d/%02d",year,month,day),
          times=sprintf("%02d:00:00",12),
          format=c(dates='y/m/d',times='h:m:s'))-1
  return(c(as.integer(as.character(years(c.date))),
           as.integer(months(c.date)),
           as.integer(days(c.date))))
}
  

get.temperature.daily<-function(year,month,day,hour,version,type='mean') {

   date.range<-c('','')
   if(hour<12) {
      pd<-previous.day(year,month,day)
      date.range<-c(sprintf("%04d-%02d-%02d:%02d",pd[1],pd[2],pd[3],hour+12),
                    sprintf("%04d-%02d-%02d:%02d",year,month,day,hour+12))
   }
   if(hour>=12) {
      nd<-next.day(year,month,day)
      date.range<-c(sprintf("%04d-%02d-%02d:%02d",year,month,day,hour-12),
                    sprintf("%04d-%02d-%02d:%02d",nd[1],nd[2],nd[3],hour-12))
   }

   t<-TWCR.get.slab.from.hourly('air.2m',date.range=date.range,version=version,type=type) 
   t$data<-apply(t$data,c(1,2),mean) # average over the times
   t$dimensions[[3]]$values<-mean(t$dimensions[[3]]$values)
   return(t)
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

plot.hour<-function(year,month,day,hour) {    

    image.name<-sprintf("%04d-%02d-%02d:%02d.png",year,month,day,hour)

    ifile.name<-sprintf("%s/%s",Imagedir,image.name)
    if(file.exists(ifile.name) && file.info(ifile.name)$size>0) return()

    t2m.normal<-get.temperature.daily(year,month,day,hour,version='3.5.4',
                                             type='normal')
    prmsl.sd<-TWCR.get.slice.at.hour('prmsl',year,month,day,hour,
                                         version='3.4.1',type='standard.deviation')
    prmsl.normal<-TWCR.get.slice.at.hour('prmsl',year,month,day,hour,version='3.5.4',
                                             type='normal')
    if(opt$obs=='T') {
       obs<-TWCR.get.obs(year,month,day,hour,version='3.5.4')
       w<-which(obs$Longitude>180)
       obs$Longitude[w]<-obs$Longitude[w]-360
    }

    versions<-c('3.5.4','3.5.6','3.6.4','3.6.6')
    captions<-c('Fixed SST, No aerosols',
                'Fixed SST, Crowley aerosols',
                'SODAsi SST, No aerosols',
                'SODAsi SST, Crowley Aerosols')
    x.min<-c(0.02,0.51,0.02,0.51)
    y.min<-c(0.05,0.05,0.525,0.525)

     png(ifile.name,
             width=1080*16/9,
             height=1080,
             bg='white',
             pointsize=24,
             type='cairo')

    for(vn in seq(1,4)) {

	prmsl<-TWCR.get.slice.at.hour('prmsl',year,month,day,hour,version=versions[vn])
	prmsl.spread<-TWCR.get.slice.at.hour('prmsl',year,month,day,hour,version=versions[vn],
						  type='spread')
	fog<-TWCR.relative.entropy(prmsl.normal,prmsl.sd,prmsl,prmsl.spread)
	fog$data[]<-1-pmin(fog.threshold,pmax(0,fog$data))/fog.threshold
	t2m<-get.temperature.daily(year,month,day,hour,version=versions[vn])
	prmsl$data[]<-prmsl$data-prmsl.normal$data
	t2m$data[]<-t2m$data-t2m.normal$data

	   pushViewport(viewport(width=0.47,height=0.45,x=x.min[vn],y=y.min[vn],
			      just=c("left","bottom"),clip='on'))

	   pushViewport(dataViewport(c(Options$lon.min,Options$lon.max),
				     c(Options$lat.min,Options$lat.max),
				      extension=0))
              grid.polygon(x=unit(c(0,1,1,0),'npc'),y=unit(c(0,0,1,1),'npc'),
                           gp=gpar(col=Options$sea.colour,fill=Options$sea.colour))
	      ip<-WeatherMap.rectpoints(Options$ice.points,Options)
	      #WeatherMap.draw.ice(ip$lat,ip$lon,icec,Options)
	      WeatherMap.draw.land(land,Options)

	    Draw.temperature(t2m,Options)
	    if(opt$obs=='T') Draw.obs(year,month,day,hour,Options,Trange=10,size=0.025,version=version)
	    Draw.pressure(prmsl,Options,colour=c(0,0,0))
	    WeatherMap.draw.fog(fog,Options)

	    Options$label<-captions[vn]
	    WeatherMap.draw.label(Options)
	 upViewport()
	 upViewport()

    }
	 Options$label<-sprintf("%04d-%02d-%02d:%02d",year,month,day,hour)
	 #WeatherMap.draw.label(Options)
         grid.text(Options$label,y=unit(0.01,'npc'),x=unit(0.99,'npc'),
                   just=c('right','bottom'))

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
