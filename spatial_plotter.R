library(rgdal)  #conatins the read/writeOGR for reading shapelies and read/writeRGDAL for reading raster
library(rgeos)  # neccessary for ggplot2::fortify.sp(); serves as a replacement for gpclib
library(maptools)  #Contains the overlay command
library(spdep)  #Contains a number of useful spatial stat functions
library(spatstat)  #Contains functions for generating random points drawn from a specific data generatin
library(raster)  #contains a number of useful functions for raster data, especially extract()
library(ggplot2)
library(reshape2)
library(scales)


plot.map <- function(spatial_file, image_file, xlab, ylab, title){
    d = read.csv(spatial_file, header=F)
    colnames(d) = c("long","lat","order","hole","piece","group","id","name")
    d$pc = d$id

    data = data.frame("id"= d[,c("id")])
    data$random = sample(1:10, nrow(d), replace=T)
    
    p1 = ggplot(data)
    p1 = p1 + geom_map(aes(fill = random, map_id = id), map = d) + scale_fill_gradient(low = "darkred", high ="red")
    p1 <- p1 + expand_limits(x = d$long, y = d$lat) + labs(x=xlab, y=ylab, title=title)
    p1 <- p1 + coord_equal()
    ggsave(p1,file=image_file)
}


