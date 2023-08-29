## ##############################################################
##
## Illusration of the concepts for k-means clustering
##

library(mva)

source(file.path(dir.util, "util_mobile_centers.R"))
    
## ##############################################################
## Main routine

## generate a random data set
x <- rbind(matrix(rnorm(200, sd = 0.3), ncol = 2),
           matrix(rnorm(50, mean = 1, sd = 0.3), ncol = 2)
           )

## number of clusters
k <- 2 

## seed centers
seed.centers <- x[c(1,2),]

## plot initial conditions
plot(x, col = '#AAAAAA',main="initial conditions",pch=19,xlab=NA,ylab=NA,cex.main=2)
points(seed.centers, col='#0000DD',pch=19)
setwd(dir.figures); export.plot(file.prefix=paste("mobile_centers_step",0,sep=""), export.formats=export.formats.plots, width=8,height=8)
print(seed.centers)

## number of iterations
for (i in c(1,2,3,4,5,6,7,8,9,10,15)) {
  #  cl <- kmeans(x, centers=seed.centers, iter.max=i)
  cl <- mobile.centers(x, centers=seed.centers, iter.max=i, plot=T,cex.main=2)
  setwd(dir.figures); export.plot(file.prefix=paste("mobile_centers_step",i,sep=""), export.formats=export.formats.plots, width=8,height=8)
  print(cl$centers)
}


## ##############################################################
## Export a single figure with all iterations
par (mfrow=c(3,3))

## plot initial conditions
plot(x, col = '#AAAAAA',main="initial conditions",pch=19,xlab=NA,ylab=NA,cex.main=2)
points(seed.centers, col='#0000DD',pch=19)

## number of iterations
for (i in c(1,2,3,4,5,6,7,8)) {
  #  cl <- kmeans(x, centers=seed.centers, iter.max=i)
  cl <- mobile.centers(x, centers=seed.centers, iter.max=i,plot=T,cex.main=2)
}
setwd(dir.figures); export.plot(file.prefix="mobile_centers", export.formats=export.formats.plots, width=13,height=9)
par(mfrow=c(1,1))


