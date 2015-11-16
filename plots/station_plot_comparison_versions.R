# Compare temperature anomalies differences from obs for all the runs

library(grid)
library(chron)

# Useful subset of stations for temperature comparisons
good.stations<-c(2,3,4,7,8,13,16,18,19,20,22,24,25,26,27,28,33,36,46,47,49,50,51)

# Load station data and Reanalysis data
load('../Yuri/0_all_station_and_travel_pressure_data.Rdata')
load('../Reanalysis.Rdata')
Dates.356<-Dates
load('../Reanalysis.366.Rdata')

# The model does not have an adequate diurnal cycle - so average each data source over +-12 hours
strip.diurnal<-function(x,x.dates) {
  result<-x
  for(i in seq_along(x)) {
     w<-which(abs(x.dates-x.dates[i])<=0.5)
     if(length(w)>1) result[i]<-mean(x[w],na.rm=TRUE)
  }
  return(result)
}

# Get a consistent loess span adjusted for missing data quantity
get.span<-function(dta) {
  w<-which(is.na(dta))
  return(0.5/(1-length(w)/length(dta)))
}

plot.station<-function(station) { # station is an integer in range 1:51

    # Reformat the station data into lists
    dates<-rep(Dates[[1]]$date,length(Dates))
    t2m<-rep(NA,length(Dates))
    prmsl<-rep(NA,length(Dates))
    v366<-list()
    v366$t2m.mean<-rep(NA,length(Dates))
    v366$t2m.spread<-rep(NA,length(Dates))
    v366$prmsl.mean<-rep(NA,length(Dates))
    v366$prmsl.spread<-rep(NA,length(Dates))
    v354<-list()
    v354$t2m.mean<-rep(NA,length(Dates))
    v354$t2m.spread<-rep(NA,length(Dates))
    v354$t2m.normal<-rep(NA,length(Dates))
    v354$prmsl.mean<-rep(NA,length(Dates))
    v354$prmsl.normal<-rep(NA,length(Dates))
    v356<-list()
    v356$t2m.mean<-rep(NA,length(Dates))
    v356$t2m.spread<-rep(NA,length(Dates))
    v356$prmsl.mean<-rep(NA,length(Dates))
    v356$prmsl.spread<-rep(NA,length(Dates))
    for(i in seq_along(Dates)) {
       dates[i]<-Dates[[i]]$date
       t2m[i]<-Dates[[i]]$station[[station]]$TA.orig
       prmsl[i]<-Dates[[i]]$station[[station]]$QFF
       v366$t2m.mean[i]<-Dates[[i]]$station[[station]]$v366$t2m.mean-273.15
       v366$t2m.spread[i]<-Dates[[i]]$station[[station]]$v366$t2m.spread
       v354$t2m.mean[i]<-Dates.356[[i]]$station[[station]]$v354$t2m.mean-273.15
       v354$t2m.spread[i]<-Dates.356[[i]]$station[[station]]$v354$t2m.spread
       v354$t2m.normal[i]<-Dates.356[[i]]$station[[station]]$v354$t2m.normal-273.15
       v356$t2m.mean[i]<-Dates.356[[i]]$station[[station]]$v356$t2m.mean-273.15
       v356$t2m.spread[i]<-Dates.356[[i]]$station[[station]]$v356$t2m.spread
       v366$prmsl.mean[i]<-Dates[[i]]$station[[station]]$v366$prmsl.mean/100
       v366$prmsl.spread[i]<-Dates[[i]]$station[[station]]$v366$prmsl.spread/100
       v354$prmsl.mean[i]<-Dates.356[[i]]$station[[station]]$v354$prmsl.mean/100
       v354$prmsl.spread[i]<-Dates.356[[i]]$station[[station]]$v354$prmsl.spread/100
       v354$prmsl.normal[i]<-Dates.356[[i]]$station[[station]]$v354$prmsl.normal/100
       v356$prmsl.mean[i]<-Dates.356[[i]]$station[[station]]$v356$prmsl.mean/100
       v356$prmsl.spread[i]<-Dates.356[[i]]$station[[station]]$v356$prmsl.spread/100
    }
    t2m<-strip.diurnal(t2m,dates)
    v366$t2m.mean<-strip.diurnal(v366$t2m.mean,dates)
    v366$t2m.spread<-strip.diurnal(v366$t2m.spread,dates)
    v354$t2m.mean<-strip.diurnal(v354$t2m.mean,dates)
    v354$t2m.spread<-strip.diurnal(v354$t2m.spread,dates)
    v354$t2m.normal<-strip.diurnal(v354$t2m.normal,dates)
    v356$t2m.mean<-strip.diurnal(v356$t2m.mean,dates)
    v356$t2m.spread<-strip.diurnal(v356$t2m.spread,dates)
   # Make anomalies
   # Match median anomaly of obs to 354
    t2m<-t2m+median(v354$t2m.mean,na.rm=TRUE)-median(t2m,na.rm=TRUE)
    t2m<-t2m-v354$t2m.normal
    v354$t2m.mean<-v354$t2m.mean-v354$t2m.normal-t2m
    v356$t2m.mean<-v356$t2m.mean-v354$t2m.normal-t2m
    v366$t2m.mean<-v366$t2m.mean-v354$t2m.normal-t2m
    prmsl<-prmsl-v354$prmsl.normal
    v354$prmsl.mean<-v354$prmsl.mean-v354$prmsl.normal
    v356$prmsl.mean<-v356$prmsl.mean-v354$prmsl.normal
    v366$prmsl.mean<-v366$prmsl.mean-v354$prmsl.normal

    # Smooth the observations to show the low-frequency behaviour
    
    lm<-loess(t2m~dates,span=get.span(t2m))
    t2m<-predict(lm,dates)
    lm<-loess(prmsl~dates,span=get.span(prmsl))
    prmsl<-predict(lm,dates)
    lm<-loess(v354$t2m.mean~dates,span=get.span(v354$t2m.mean),na.action = na.exclude)
    v354$t2m.mean<-predict(lm,dates)
    lm<-loess(v354$prmsl.mean~dates,span=get.span(v354$prmsl.mean),na.action = na.exclude)
    v354$prmsl.mean<-predict(lm,dates)
    lm<-loess(v356$t2m.mean~dates,span=get.span(v356$t2m.mean),na.action = na.exclude)
    v356$t2m.mean<-predict(lm,dates)
    lm<-loess(v356$prmsl.mean~dates,span=get.span(v356$prmsl.mean),na.action = na.exclude)
    v356$prmsl.mean<-predict(lm,dates)
    lm<-loess(v366$t2m.mean~dates,span=get.span(v366$t2m.mean),na.action = na.exclude)
    v366$t2m.mean<-predict(lm,dates)
    lm<-loess(v366$prmsl.mean~dates,span=get.span(v366$prmsl.mean),na.action = na.exclude)
    v366$prmsl.mean<-predict(lm,dates)



            # Analysis means
	     gp=gpar(col=rgb(1,0,0,0.1),fill=rgb(1,0,0,0.1))
	     grid.lines(x=unit(dates,'native'),
			 y=unit(v366$t2m.mean,'native'),
			 gp=gp)
	     gp=gpar(col=rgb(0,0,1,0.1),fill=rgb(0,0,1,0.1))
	     grid.lines(x=unit(dates,'native'),
			 y=unit(v354$t2m.mean,'native'),
			 gp=gp)



}

pdf(file="Stations.sst_comparison_all.pdf",
    width=10*sqrt(2),height=10,family='Helvetica',
    paper='special',pointsize=18,bg='white')

    date.range<-chron(dates=c("1815/01/01","1817/12/31"),
			times=c("00:0:01","23:59:59"),
			format=c(dates = "y/m/d", times = "h:m:s"))
    tics=pretty(date.range)
    ticl=attr(tics,'labels')
    
    gpe<-list()
    gpe$v354<-gpar(col=rgb(0.8,0.8,1,1),fill=rgb(0.8,0.8,1,1))
    gpe$v356<-gpar(col=rgb(0.6,0.6,1,1),fill=rgb(0.6,0.6,1,1))
    gpe$v366<-gpar(col=rgb(0.4,0.4,1,1),fill=rgb(0.4,0.4,1,1))
                     
    pushViewport(viewport(width=1.0,height=1.0,x=0.0,y=0.0,
			  just=c("left","bottom"),name="Page",clip='off'))
       pushViewport(plotViewport(margins=c(4,6,0,0)))
	  pushViewport(dataViewport(date.range,c(-5,5),clip='off'))
	     grid.xaxis(at=as.numeric(tics),label=ticl,main=T)
	     grid.text('Date',y=unit(-3,'lines'))
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


             for(station in good.stations) { # 23
                plot.station(station)
             }

	  popViewport()
       popViewport()
    popViewport()
