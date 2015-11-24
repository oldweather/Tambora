# Make 1831-1850 normals for 20CR 3.5.1 and store them for 3.5.4 (Tambora)
# Special set of normals for the Tambora analysis.

library(GSDF.TWCR)
library(chron)

version<-'3.5.1'
start.year<-1951
end.year<-1980

# The only tricky issue is leap years - follow
#  20CRv2 and have no climatology for Feb 29
# Need a function to take out the Feb 29th data
Unleap<-function(f) {
   nl1<-length(f$dimensions[[1]]$values)
   nl2<-length(f$dimensions[[2]]$values)
   per.day<-length(f$dimensions[[3]]$values)/366
   w<-c(seq(1,nl1*nl2*(31+28)*per.day),
        seq(nl1*nl2*(31+30)*per.day,nl1*nl2*366*per.day))
   f$data<-array(data=f$data[w],dim=c(nl1,nl2,365*per.day))
   w<-seq((31+28)*per.day,(31+29)*per.day)
   f$dimensions[[3]]$values<-f$dimensions[[3]]$values[-w]
   return(f)
}

if(leap.year(start.year)) {
   stop("Can't start climatology period on a leap year")
}

var.clim<-function(parameter) {

    n<-TWCR.get.slab.from.hourly(parameter,
			  date.range=c(sprintf("%04d-01-01:00",start.year),
				       sprintf("%04d-12-31-23",start.year)),
			  version=version)

    for(year in seq(start.year+1,end.year)) {
       f<-TWCR.get.slab.from.hourly(parameter,
			  date.range=c(sprintf("%04d-01-01:00",year),
				       sprintf("%04d-12-31-23",year)),
			  version=version)
       if(leap.year(year)) f<-Unleap(f)
       n$data[]<-n$data+f$data
    }
    n$data[]<-n$data/(end.year-start.year+1)

    # Change the dates to 1981
    s.date<-chron(dates="1981/01/01",times="00:00:00",
          format=c(dates='y/m/d',times='h:m:s'))
    n$dimensions[[3]]$values<-as.chron(as.numeric(n$dimensions[[3]]$values) - 
                                       as.numeric(n$dimensions[[3]]$values[1]) + 
                                       as.numeric(s.date))

    # Output the normal to the appropriate file.
    fn<-sprintf("/project/projectdirs/m958/netCDF.data/20CR_v3.5.6/hourly/normals/%s.nc",parameter)
    s<-GSDF.ncdf.write(n,file.name=fn,name=parameter)
    warnings()
    if(is.null(s)) stop('Normals write failed')

}

lapply(c('air.2m','prmsl'),var.clim)
