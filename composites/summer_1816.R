# Make summer (JJA) 1816 Temperature anomaly composites

library(GSDF.TWCR)

get.anomaly.composite<-function(date.range,version) {

   n<-TWCR.get.slab.from.hourly('air.2m',
            date.range<-date.range,
            type='normal',version='3.5.4')
   t<-TWCR.get.slab.from.hourly('air.2m',
               date.range<-date.range,
               type='mean',version=version)
   p<-TWCR.get.slab.from.hourly('prmsl',
               date.range<-date.range,
               type='mean',version=version)
   s<-TWCR.get.slab.from.hourly('prmsl',
               date.range<-date.range,
               type='spread',version=version)
   t$data[]<-t$data-n$data
   c<-GSDF.reduce.1d(t,'time',mean)
   v<-GSDF.reduce.1d(p,'time',var)
   s<-GSDF.reduce.1d(s,'time',mean)
   r<-s
   r$data[]<-v$data/(v$data+s$data**2)
   c2<-c
   c2$meta$pole.lat<--40
   c2$meta$pole.lon<-0
   c<-GSDF.regrid.2d(c,c2)
   r<-GSDF.regrid.2d(r,c2)
   w<-which(r$data<0.7)
   is.na(c$data[w])<-TRUE
   return(c)
}

c<-get.anomaly.composite(c('1816-06-01:00','1816-08-31:23'),'3.5.4')
c2<-get.anomaly.composite(c('1816-06-01:00','1816-08-31:23'),'3.5.5'
#GSDF.pplot.2d(c,d,x.range=c(160,210),y.range=c(-20,25),levels=seq(-4,4,.2))
