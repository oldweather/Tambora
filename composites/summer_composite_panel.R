# Make summer (JJA) 1816 Temperature anomaly composites
# Portrait format for publication.

library(GSDF.TWCR)

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
   return(c)
}

for(year in c(1816)) {

   c<-get.anomaly.composite(c(
        sprintf("%04d-06-01:00",year),
        sprintf("%04d-08-31:23",year)),'3.5.4')
   c2<-get.anomaly.composite(c(
        sprintf("%04d-06-01:00",year),
        sprintf("%04d-08-31:23",year)),'3.5.6')
   d<-c
   d$data[]<-c2$data-c$data
   png(filename=sprintf("summer_%04d_panel.png",year),width=600,height=1024,pointsize=24)
   trellis.par.set('axis.text',gpar(cex=2))
    grid.newpage()
    pushViewport(viewport(x = 0, y = 0.66, width = 1, height = 0.33,
        just = c("left", "bottom")))
    m1 <- GSDF.plot.2d(c2,draw = F, x.range=c(165,205),y.range=c(40,65),levels=seq(-5,5,.2),
                 x.scale=0,y.scale=0,x.label='',y.label='')
    print(m1, newpage = F)
    popViewport()
    pushViewport(viewport(x = 0, y = 0.33, width = 1, height = 0.33,
        just = c("left", "bottom")))
    m2 <- GSDF.plot.2d(c,draw = F, x.range=c(165,205),y.range=c(40,65),levels=seq(-5,5,.2),
                 x.scale=0,y.scale=0,x.label='',y.label='')
    print(m2, newpage = F)
    popViewport()
    pushViewport(viewport(x = 0, y = 0, width = 1, height = 0.33,
        just = c("left", "bottom")))
    m3 <- GSDF.plot.2d(d,draw = F, x.range=c(165,205),y.range=c(40,65),levels=seq(-5,5,.2),
                 x.scale=0,y.scale=0,x.label='',y.label='')
    print(m3, newpage = F)
    popViewport()

   dev.off()
}
