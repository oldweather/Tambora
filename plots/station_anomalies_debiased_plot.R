# Plot Pressure and temperature anomalies for one station
#  Apply simple bias adjustment to obs - match median anomaly
#  to 354 median anomaly.

library(grid)
library(chron)

# Load station data and Reanalysis data
load('../Yuri/0_all_station_and_travel_pressure_data.Rdata')
load('../Reanalysis.Rdata')

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
    # Match median anomaly of obs to 354
    t2m<-t2m+median(v354$t2m.mean,na.rm=TRUE)-median(t2m,na.rm=TRUE)
    prmsl<-prmsl+median(v354$prmsl.mean,na.rm=TRUE)-median(prmsl,na.rm=TRUE)


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

pdf(file="Station.anomalies.debiased.pdf",
    width=10*sqrt(2),height=10,family='Helvetica',
    paper='special',pointsize=18)

for(station in seq(1,51)) {
   grid.newpage()
   plot.station(station)
}
