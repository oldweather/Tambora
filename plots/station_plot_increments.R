# Plot Pressure and temperature increments for one station

library(grid)
library(chron)

good.stations<-c(2,3,4,7,8,13,16,18,19,20,22,24,25,26,27,28,33,36,46,47,49,50,51)

# Load station data and Reanalysis data
load('../Yuri/0_all_station_and_travel_pressure_data.Rdata')
load('../Reanalysis.Rdata')
Dates.356<-Dates
load('../Reanalysis.366.Rdata')
Dates.366<-Dates
load('../Reanalysis.354.fd.Rdata')
Dates.354.fd<-Dates
load('../Reanalysis.366.fd.Rdata')
Dates.366.fd<-Dates


plot.station<-function(station) { # station is an integer in range 1:51

    # Reformat the station data into lists
    dates<-rep(Dates[[1]]$date,length(Dates))
    t2m<-rep(NA,length(Dates))
    prmsl<-rep(NA,length(Dates))
    v354<-list()
    v354$t2m.mean<-rep(NA,length(Dates))
    v354$t2m.fd<-rep(NA,length(Dates))
    v354$prmsl.mean<-rep(NA,length(Dates))
    v354$prmsl.fd<-rep(NA,length(Dates))
    v366<-list()
    v366$t2m.mean<-rep(NA,length(Dates))
    v366$t2m.fd<-rep(NA,length(Dates))
    v366$prmsl.mean<-rep(NA,length(Dates))
    v366$prmsl.fd<-rep(NA,length(Dates))
    for(i in seq_along(Dates)) {
       dates[i]<-Dates[[i]]$date
       v354$t2m.mean[i]<-Dates.356[[i]]$station[[station]]$v354$t2m.mean
       v354$t2m.fd[i]<-Dates.354.fd[[i]]$station[[station]]$v354$t2m.mean
       v354$prmsl.mean[i]<-Dates.356[[i]]$station[[station]]$v354$prmsl.mean/100
       v354$prmsl.fd[i]<-Dates.354.fd[[i]]$station[[station]]$v354$prmsl.mean/100
       v366$t2m.mean[i]<-Dates.366[[i]]$station[[station]]$v366$t2m.mean
       v366$t2m.fd[i]<-Dates.366.fd[[i]]$station[[station]]$v366$t2m.mean
       v366$prmsl.mean[i]<-Dates.366[[i]]$station[[station]]$v366$prmsl.mean/100
       v366$prmsl.fd[i]<-Dates.366.fd[[i]]$station[[station]]$v366$prmsl.mean/100
    }

    date.range<-chron(dates=c("1815/01/01","1817/12/31"),
			times=c("00:0:01","23:59:59"),
			format=c(dates = "y/m/d", times = "h:m:s"))
    tics=pretty(date.range)
    ticl=attr(tics,'labels')
    
    gpe<-list()
    gpe$v354<-gpar(col=rgb(1,0.4,0.4,1),fill=rgb(0.8,0.8,1,1))
    gpe$v366<-gpar(col=rgb(0.4,0.4,1,1),fill=rgb(0.4,0.4,1,1))
                     
    # Pressure along the bottom with x axis
    pmin<- -1
    for(ensda in c('v354','v366')) {
        v<-get(ensda)
	rmin<-min(v$prmsl.mean-v$prmsl.fd,na.rm=TRUE)
	pmin<-min(pmin,rmin,na.rm=TRUE)
    }
    omin<-quantile(prmsl,probs=0.05,na.rm=TRUE)
    pmin<-min(pmin,omin,na.rm=TRUE)
    pmax<-1
    for(ensda in c('v354','v366')) {
        v<-get(ensda)
	rmax<-max(v$prmsl.mean-v$prmsl.fd,na.rm=TRUE)
	pmax<-max(pmax,rmax,na.rm=TRUE)
    }
    pushViewport(viewport(width=1.0,height=0.55,x=0.0,y=0.0,
			  just=c("left","bottom"),name="Page",clip='off'))
       pushViewport(plotViewport(margins=c(4,6,0,0)))
	  pushViewport(dataViewport(date.range,c(pmin,pmax),clip='off'))
	     grid.xaxis(at=as.numeric(tics),label=ticl,main=T)
	     grid.text('Date',y=unit(-3,'lines'))
	     grid.yaxis(main=T)
	     grid.text('Pressure increments (hPa)',x=unit(-4,'lines'),rot=90)

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

	     for(ensda in c('v354','v366')) { # Order matters 
             v<-get(ensda)
	     grid.points(x=unit(dates,'native'),
			 y=unit(v$prmsl.mean-v$prmsl.fd,'native'),
                         size=unit(0.005,'npc'),
                         pch=20,
			 gp=gpe[[ensda]])
                 
              }
            
	  popViewport()
       popViewport()
    popViewport()

    # AT on top
    tmin<- -1
    for(ensda in c('v354','v366')) {
        v<-get(ensda)
	rmin<-min(v$t2m.mean-v$t2m.fd,na.rm=TRUE)
	tmin<-min(tmin,rmin,na.rm=TRUE)
    }
    omin<-quantile(t2m,probs=0.05,na.rm=TRUE)
    tmin<-min(tmin,omin,na.rm=TRUE)
    tmax<- 1
    for(ensda in c('v354','v366')) {
        v<-get(ensda)
	rmax<-max(v$t2m.mean-v$t2m.fd,na.rm=TRUE)
	tmax<-max(tmax,rmax,na.rm=TRUE)
    }
    pushViewport(viewport(width=1.0,height=0.45,x=0.0,y=0.55,
			  just=c("left","bottom"),name="Page",clip='off'))
       pushViewport(plotViewport(margins=c(0,6,0,0)))
	  pushViewport(dataViewport(date.range,c(tmin,tmax),clip='off'))
	     grid.yaxis(main=T)
	     grid.text('Air temperature increment (C)',x=unit(-4,'lines'),rot=90)

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

	     for(ensda in c('v354','v366')) { # Order matters
             v<-get(ensda)
	     grid.points(x=unit(dates,'native'),
			 y=unit(v$t2m.mean-v$t2m.fd,'native'),
                         size=unit(0.005,'npc'),
                         pch=20,
			 gp=gpe[[ensda]])
               }
	  popViewport()
       popViewport()
    popViewport()

    # Name the station
    grid.text(PP[[station]]$Station.name[1],
	      x=unit(0.03,'npc'),y=unit(0.03,'npc'),just<-c('left','bottom'))

}

pdf(file="Increments.daily.pdf",
    width=10*sqrt(2),height=10,family='Helvetica',
    paper='special',pointsize=18,bg='white')

for(station in good.stations) {
   grid.newpage()
   plot.station(station)
}
