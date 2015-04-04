# Make summer (JJA) 1815-7 Precipitation anomaly composites

library(GSDF.TWCR)

get.anomaly.composite<-function(date.range,version) {

   n<-TWCR.get.slab.from.hourly('prate',
            date.range<-date.range,
            type='normal',version='3.5.4')
   t<-TWCR.get.slab.from.hourly('prate',
               date.range<-date.range,
               type='mean',version=version)
   p<-TWCR.get.slab.from.hourly('prmsl',
               date.range<-date.range,
               type='mean',version=version)
   s<-TWCR.get.slab.from.hourly('prmsl',
               date.range<-date.range,
               type='spread',version=version)
   t$data[]<-t$data/n$data
   c<-GSDF.reduce.1d(t,'time',mean)
   v<-GSDF.reduce.1d(p,'time',var)
   s<-GSDF.reduce.1d(s,'time',mean)
   r<-s
   r$data[]<-v$data/(v$data+s$data**2)
   c2<-c
   c2$meta$pole.lat<-90
   c2$meta$pole.lon<-0
   c<-GSDF.regrid.2d(c,c2)
   r<-GSDF.regrid.2d(r,c2)
   w<-which(r$data<0.5)
   is.na(c$data[w])<-TRUE
   return(c)
}

for(year in c(1815,1816,1817)) {

   c<-get.anomaly.composite(c(
        sprintf("%04d-06-01:00",year),
        sprintf("%04d-08-31:23",year)),'3.5.4')
   c2<-get.anomaly.composite(c(
        sprintf("%04d-06-01:00",year),
        sprintf("%04d-08-31:23",year)),'3.5.6')
   d<-c
   d$data[]<-c2$data/c$data
   c$data[]<-log(c$data)
   d$data[]<-log(d$data)
   png(filename=sprintf("summer_prate_%04d.png",year),width=1024,height=768,pointsize=24)
   GSDF.pplot.2d(c,d,x.range=c(90,220),y.range=c(30,80),levels=seq(-1.2,1.2,.1),
                 x.scale=0,y.scale=0,x.label='',y.label='')
   dev.off()
}