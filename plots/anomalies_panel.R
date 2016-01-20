# Plot Pressure and temperature anomalies for one station
# Scatter plot and time-series

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



    date.range<-chron(dates=c("1816/01/01","1816/12/31"),
			times=c("00:0:01","23:59:59"),
			format=c(dates = "y/m/d", times = "h:m:s"))
    tics=pretty(date.range,min.n=5)
    ticl=attr(tics,'labels')
    
    gpe<-list()
    gpe$v354<-gpar(col=rgb(0.8,0.8,1,1),fill=rgb(0.8,0.8,1,1))
    gpe$v355<-gpar(col=rgb(0.6,0.6,1,1),fill=rgb(0.6,0.6,1,1))
    gpe$v356<-gpar(col=rgb(0.4,0.4,1,1),fill=rgb(0.4,0.4,1,1))
                     
    # Pressure along the bottom with x axis
    pmin<- -20
    for(ensda in c('v356')) {
        v<-get(ensda)
	rmin<-min(v$prmsl.mean-v$prmsl.spread*2,na.rm=TRUE)
	pmin<-min(pmin,rmin,na.rm=TRUE)
    }
    omin<-quantile(prmsl,probs=0.05,na.rm=TRUE)
    pmin<-min(pmin,omin,na.rm=TRUE)
    pmax<-10
    for(ensda in c('v356')) {
        v<-get(ensda)
	rmax<-max(v$prmsl.mean+v$prmsl.spread*2,na.rm=TRUE)
	pmax<-max(pmax,rmax,na.rm=TRUE)
    }
    omax<-quantile(prmsl,probs=0.95,na.rm=TRUE)
    pmax<-max(pmax,omax,na.rm=TRUE)
    pushViewport(viewport(width=1.0,height=0.3,x=0.0,y=0.0,
			  just=c("left","bottom"),name="Page",clip='off'))
       pushViewport(plotViewport(margins=c(4,6,0,0)))
	  pushViewport(dataViewport(date.range,c(pmin,pmax),clip='off'))
	     grid.xaxis(at=as.numeric(tics),label=ticl,main=T)
	     grid.text('Date',y=unit(-3,'lines'))
	     grid.yaxis(main=T)
	     grid.text('Sea-level\npressure (hPa)',x=unit(-4,'lines'),rot=90)

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

	     for(ensda in c('v356')) { # Order matters 
                 v<-get(ensda)
		 for(i in seq_along(dates)) {
                    if(dates[i]<date.range[1] ||
                       dates[i]>date.range[2]) next
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
             w<-which(dates>=date.range[1]&
                      dates<=date.range[2])
	     grid.points(x=unit(dates[w],'native'),
			 y=unit(prmsl[w],'native'),
			 size=unit(0.005,'npc'),
			 pch=20,
			 gp=gp)
	  popViewport()
       popViewport()
    popViewport()

    # AT on top
    tmin<- -5
    for(ensda in c('v356')) {
        v<-get(ensda)
	rmin<-min(v$t2m.mean-v$t2m.spread*2,na.rm=TRUE)
	tmin<-min(tmin,rmin,na.rm=TRUE)
    }
    omin<-quantile(t2m,probs=0.05,na.rm=TRUE)
    tmin<-min(tmin,omin,na.rm=TRUE)
    tmax<- 5
    for(ensda in c('v356')) {
        v<-get(ensda)
	rmax<-max(v$t2m.mean+v$t2m.spread*2,na.rm=TRUE)
	tmax<-max(tmax,rmax,na.rm=TRUE)
    }
    omax<-quantile(t2m,probs=0.95,na.rm=TRUE)
    tmax<-max(tmax,omax,na.rm=TRUE)
    pushViewport(viewport(width=1.0,height=0.2,x=0.0,y=0.3,
			  just=c("left","bottom"),name="Page",clip='off'))
       pushViewport(plotViewport(margins=c(0,6,0,0)))
	  pushViewport(dataViewport(date.range,c(tmin,tmax),clip='off'))
	     grid.yaxis(main=T)
	     grid.text('Air temperature\n(C)',x=unit(-4,'lines'),rot=90)

	  # Mark the eruption
	     gp=gpar(col=rgb(0,0,0,1),fill=rgb(0,0,0,1))

	     p.y<-c(-200,300)
	     p.x<-chron(dates=c("1815/04/10","1815/04/10"),
			times=rep("12:00:00",2),
			format=c(dates = "y/m/d", times = "h:m:s"))
	     grid.lines(x=unit(p.x,'native'),
			  y=unit(p.y,'native'),
			  gp=gp)


	     # Analysis spreads

	     for(ensda in c('v356')) { # Order matters
                 v<-get(ensda)
		 for(i in seq_along(dates)) {
                    if(dates[i]<date.range[1] ||
                       dates[i]>date.range[2]) next
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
             w<-which(dates>=date.range[1]&
                      dates<=date.range[2])
	     grid.points(x=unit(dates[w],'native'),
			 y=unit(t2m[w],'native'),
			 size=unit(0.005,'npc'),
			 pch=20,
			 gp=gp)
	  popViewport()
       popViewport()
    popViewport()

# SLP scatter
    pushViewport(viewport(width=0.5,height=0.5,x=0.0,y=0.5,
			  just=c("left","bottom"),name="Page",clip='off'))
       pushViewport(plotViewport(margins=c(4,6,0,0)))
	  pushViewport(dataViewport(c(pmin,pmax),c(pmin,pmax),clip='off'))
	     grid.xaxis()
	     grid.text('Observed SLP',y=unit(-3,'lines'))
	     grid.yaxis(main=T)
	     grid.text('Reanalysis SLP',x=unit(-4,'lines'),rot=90)

	     for(ensda in c('v356')) { # Order matters 
                 v<-get(ensda)
		 for(i in seq_along(dates)) {
                    if(dates[i]<date.range[1] ||
                       dates[i]>date.range[2]) next
		    x<-c(prmsl[i],prmsl[i])
		    y<-c(v$prmsl.mean[i]-
			   v$prmsl.spread[i]*2,
			 v$prmsl.mean[i]+
			  v$prmsl.spread[i]*2)
		    grid.lines(x=unit(x,'native'),
				 y=unit(y,'native'),
			      gp=gpe[[ensda]])
		  }
             }
              grid.lines(x=unit(c(pmin,pmax),'native'),
                         y=unit(c(pmin,pmax),'native'),
                         gp=gpar(col=rgb(.8,.8,.8,1)))
            
	  popViewport()
       popViewport()
    popViewport()

# AT scatter
    pushViewport(viewport(width=0.5,height=0.5,x=0.5,y=0.5,
			  just=c("left","bottom"),name="Page",clip='off'))
       pushViewport(plotViewport(margins=c(4,6,0,0)))
	  pushViewport(dataViewport(c(tmin,tmax),c(tmin,tmax),clip='off'))

	     grid.xaxis()
	     grid.text('Observed AT',y=unit(-3,'lines'))
	     grid.yaxis(main=T)
	     grid.text('Reanalysis AT',x=unit(-4,'lines'),rot=90)

	     for(ensda in c('v356')) { # Order matters
                 v<-get(ensda)
		 for(i in seq_along(dates)) {
                    if(dates[i]<date.range[1] ||
                       dates[i]>date.range[2]) next
		    x<-c(t2m[i],t2m[i])
		    y<-c(v$t2m.mean[i]-
			   v$t2m.spread[i]*2,
			 v$t2m.mean[i]+
			  v$t2m.spread[i]*2)
		    grid.lines(x=unit(x,'native'),
				 y=unit(y,'native'),
			      gp=gpe[[ensda]])
		  }
              }

              grid.lines(x=unit(c(tmin,tmax),'native'),
                         y=unit(c(tmin,tmax),'native'),
                         gp=gpar(col=rgb(.8,.8,.8,1)))
    
	  popViewport()
       popViewport()
    popViewport()

}

pdf(file="anomalies.panel.pdf",
    width=10*sqrt(2),height=10,family='Helvetica',
    paper='special',pointsize=18,bg='white')

for(station in seq(18,18)) { # 18 is Geneva
   grid.newpage()
   plot.station(station)
}
