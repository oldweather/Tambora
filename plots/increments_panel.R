#  Calculate reanalysis increments at each station
#   plot them for 354 and 366.

library(grid)
library(chron)

# Divide stations into high-quality and 'other'
good.stations<-c(2,3,4,7,8,13,16,18,19,20,22,24,25,26,27,28,33,36,46,47,49,50,51)
bad.stations<-c(1,5,6,9,10,11,12,14,15,17,21,23,29,30,31,32,34,35,37,38,39,40,41,42,43,44,45,48)

# Load station data and Reanalysis data
load('../Yuri/0_all_station_and_travel_pressure_data.Rdata')
load('../Reanalysis.Rdata')
Dates.356<-Dates
load('../Reanalysis.366.Rdata')
Dates.366<-Dates
load('../Reanalysis.354.fd.Rdata')
Dates.354.fd<-Dates
load('../Reanalysis.356.fd.Rdata')
Dates.356.fd<-Dates

cor.station<-function(station) { # station is an integer in range 1:51

    # Reformat the station data into lists
    dates<-rep(Dates[[1]]$date,length(Dates))
    t2m<-rep(NA,length(Dates))
    prmsl<-rep(NA,length(Dates))
    v354<-list()
    v354$t2m.mean<-rep(NA,length(Dates))
    v354$t2m.fd<-rep(NA,length(Dates))
    v354$prmsl.mean<-rep(NA,length(Dates))
    v354$prmsl.fd<-rep(NA,length(Dates))
    v356<-list()
    v356$t2m.mean<-rep(NA,length(Dates))
    v356$t2m.fd<-rep(NA,length(Dates))
    v356$prmsl.mean<-rep(NA,length(Dates))
    v356$prmsl.fd<-rep(NA,length(Dates))
    for(i in seq_along(Dates)) {
       dates[i]<-Dates[[i]]$date
       v354$t2m.mean[i]<-Dates.356[[i]]$station[[station]]$v354$t2m.mean
       v354$t2m.fd[i]<-Dates.354.fd[[i]]$station[[station]]$v354$t2m.mean
       v354$prmsl.mean[i]<-Dates.356[[i]]$station[[station]]$v354$prmsl.mean/100
       v354$prmsl.fd[i]<-Dates.354.fd[[i]]$station[[station]]$v354$prmsl.mean/100
       v356$t2m.mean[i]<-Dates.356[[i]]$station[[station]]$v356$t2m.mean
       v356$t2m.fd[i]<-Dates.356.fd[[i]]$station[[station]]$v356$t2m.mean
       v356$prmsl.mean[i]<-Dates.356[[i]]$station[[station]]$v356$prmsl.mean/100
       v356$prmsl.fd[i]<-Dates.356.fd[[i]]$station[[station]]$v356$prmsl.mean/100
    }


    slp.354<-sd(v354$prmsl.mean-v354$prmsl.fd,na.rm=TRUE)
    t2m.354<-sd(v354$t2m.mean-v354$t2m.fd,na.rm=TRUE)
    slp.356<-sd(v356$prmsl.mean-v356$prmsl.fd,na.rm=TRUE)
    t2m.356<-sd(v356$t2m.mean-v356$t2m.fd,na.rm=TRUE)

    gp<-gpar(col=rgb(0,0,0,1),fill=rgb(0,0,0,1))
    if(slp.354>slp.356) gp<-gpar(col=rgb(1,0,0,1),fill=rgb(1,0,0,1))
    grid.lines(x=unit(c(0.05,0.4),'native'),
                y=unit(c(slp.354,slp.356),'native'),gp=gp)
    grid.points(x=unit(c(0.05,0.4),'native'),
	    y=unit(c(slp.354,slp.356),'native'),
	    size=unit(0.01,'npc'),pch=20,gp=gp)
   return(slp.356)

}

pdf(file="Increments.panel.pdf",
    width=10/sqrt(2),height=10,family='Helvetica',
    paper='special',pointsize=18,bg='white')

    pushViewport(viewport(width=1.0,height=1.0,x=0.0,y=0.0,
                          just=c("left","bottom"),name="Page",clip='off'))
       pushViewport(plotViewport(margins=c(3,4,0,0)))
          pushViewport(dataViewport(c(0,1),c(1,6),clip='off'))
             grid.xaxis(at=c(0.05,0.4),label=c('No forcing','With forcing'),main=T)
             grid.yaxis(main=T)
             grid.text('Pressure increment sd (hPa)',x=unit(-3,'lines'),rot=90)

            gp<-gpar(col=rgb(0.7,0.7,0.7,1),fill=rgb(0.7,0.7,0.7,1))
             with.by.name<-list()
             for(station in good.stations) {
                with.by.name[[PP[[station]]$Station.name[1]]]<-cor.station(station)
             }
             w<-order(unlist(with.by.name))
             for(sti in seq_along(with.by.name)) {
                yp<-0.0+(sti/length(w))*0.8
                stnname<-names(with.by.name)[w[sti]]
                if(length(grep('itenice',stnname))>0) stnname<-'Zitenice' # Unprintable
                grid.text(stnname,
                  x=unit(0.5,'npc'),y=unit(yp,'npc'),just<-c('left','center'))
                grid.lines(x=unit(c(0.43,0.49),'npc'),
                           y=unit.c(unit(with.by.name[[w[sti]]],'native'),
                                    unit(yp,'npc')),
                           gp=gp)
             }


          popViewport()
       popViewport()
    popViewport()
dev.off()

