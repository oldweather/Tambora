# Make summer (JJA) 1815 Monthly temperature anomaly composites

library(GSDF.TWCR)
library(grid)

get.anomaly.composite<-function(date.range,version) {

   n<-TWCR.get.slab.from.hourly('air.2m',
            date.range<-date.range,
            type='normal',version='3.5.6')
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
   c2$meta$pole.lat<-90
   c2$meta$pole.lon<-0
   c<-GSDF.regrid.2d(c,c2)
   r<-GSDF.regrid.2d(r,c2)
   w<-which(r$data<0.5)
   is.na(c$data[w])<-TRUE
   w<-which(c$data>5)
   c$data[w]<-5
   w<-which(c$data< -5)
   c$data[w]<- -5
   return(c)
}

c1<-get.anomaly.composite(c("1815-06-01:00","1815-06-30:23"),'3.5.6')
c2<-get.anomaly.composite(c("1815-07-01:00","1815-07-31:23"),'3.5.6')
c3<-get.anomaly.composite(c("1815-08-01:00","1815-08-31:23"),'3.5.6')

png(filename="1815_monthly.png",width=412,height=768,pointsize=24)
   p1<-GSDF.plot.2d(c1,x.range=c(150,215),y.range=c(30,70),
                    levels=seq(-5,5,.2),draw=FALSE,
                    x.scale=0,y.scale=0,x.label='',y.label='')
   pushViewport(viewport(width=1.0,height=0.33,x=0.0,y=0.67,
                       just=c("left","bottom"),name="Page",clip='off'))
      print(p1,newpage=FALSE)
   popViewport()
   p2<-GSDF.plot.2d(c2,x.range=c(150,215),y.range=c(30,70),
                    levels=seq(-5,5,.2),draw=FALSE,
                    x.scale=0,y.scale=0,x.label='',y.label='')
   pushViewport(viewport(width=1.0,height=0.33,x=0.0,y=0.33,
                       just=c("left","bottom"),name="Page",clip='off'))
      print(p2,newpage=FALSE)
   popViewport()
   p3<-GSDF.plot.2d(c3,x.range=c(150,215),y.range=c(30,70),
                    levels=seq(-5,5,.2),draw=FALSE,
                    x.scale=0,y.scale=0,x.label='',y.label='')
   pushViewport(viewport(width=1.0,height=0.33,x=0.0,y=0.00,
                       just=c("left","bottom"),name="Page",clip='off'))
      print(p3,newpage=FALSE)
   popViewport()
dev.off()
