# Plot Pressure and temperature anomalies for one station

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


    date.range<-chron(dates=c("1815/01/01","1817/12/31"),
			times=c("00:0:01","23:59:59"),
			format=c(dates = "y/m/d", times = "h:m:s"))
    tics=pretty(date.range)
    ticl=attr(tics,'labels')
    
    gpe<-list()
    gpe$v354<-gpar(col=rgb(0.4,0.4,1,.5),fill=rgb(0.4,0.4,1,0.5))
    gpe$v355<-gpar(col=rgb(0.6,0.6,1,1),fill=rgb(0.6,0.6,1,0.5))
    gpe$v356<-gpar(col=rgb(0.4,0.4,1,1),fill=rgb(0.4,0.4,1,0.5))
                     
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
    pushViewport(viewport(width=1.0,height=0.5,x=0.0,y=0.0,
			  just=c("left","bottom"),name="Page",clip='off'))
       pushViewport(plotViewport(margins=c(4,6,0,0)))
	  pushViewport(dataViewport(c(pmin,pmax),c(pmin,pmax),clip='off'))
	     grid.xaxis()
	     grid.text('Observed SLP',y=unit(-3,'lines'))
	     grid.yaxis(main=T)
	     grid.text('Reanalysis SLP',x=unit(-4,'lines'),rot=90)

	     for(ensda in c('v354')) { # Order matters 
                 v<-get(ensda)
		 for(i in seq_along(dates)) {
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
    pushViewport(viewport(width=1.0,height=0.5,x=0.0,y=0.5,
			  just=c("left","bottom"),name="Page",clip='off'))
       pushViewport(plotViewport(margins=c(4,6,0,0)))
	  pushViewport(dataViewport(c(tmin,tmax),c(tmin,tmax),clip='off'))

	     grid.xaxis()
	     grid.text('Observed AT',y=unit(-3,'lines'))
	     grid.yaxis(main=T)
	     grid.text('Reanalysis AT',x=unit(-4,'lines'),rot=90)

	     for(ensda in c('v354')) { # Order matters
                 v<-get(ensda)
		 for(i in seq_along(dates)) {
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

    # Name the station
    grid.text(PP[[station]]$Station.name[1],
	      x=unit(0.03,'npc'),y=unit(0.01,'npc'),just<-c('left','bottom'))

}

pdf(file="Station.anomalies.scatter.pdf",
    width=10,height=10*sqrt(2),family='Helvetica',
    paper='special',pointsize=18)

for(station in seq(1,51)) {
   grid.newpage()
   plot.station(station)
}
