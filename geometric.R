################################################################
##
## Illustrations of the geometric distribution
##
## Running this script requires to first run the script config.R
##    source('http://pedagogix-tagc.univ-mrs.fr/courses/stat
##
## Single-command execution:
##   source(file.path(dir.R.files, 'geometric.R'))


## Parameters
p <- 0.25
n <- 30
x <- 0:n
in.colors <- F

test.x <- 3

## color specification (actually I don't use colors for the book
if (in.colors) {
  plot.colors <- c('density'='blue',
              'CDF'='green',
              'dCDF'='violet'
              )
} else {
  plot.colors <- c('density'='black',
              'CDF'='black',
              'dCDF'='black'
              )
}

## drawing parameters
X11(width=7,height=5)
par(font.lab=2)
par(font.axis=2)
par(font.main=2)
par(cex=1.2)
par(cex.main=1.2)
par(cex.axis=1.2)
par(cex.lab=1.2)
grid.color <- '#CCCCCC'

## grid sizes
x.grid <- 0:6*5
y.grid <- 0:20*0.05

## ##############################################################
##
## Waiting time before the first success in a Bernouli schema
##

################################################################
## Plot the probability mass function
plot(x - 0.5,
     dgeom(x,p),
     ylim=c(0,1),
     main=paste('Geometric probability mass function; p=',format(p),', n=', format(n)),
     xlab='Waiting time',
     ylab='probability mass function P(X=x)',
     panel.first=c(abline(v=x.grid,col=grid.color),
       abline(h=y.grid,col=grid.color)),
     col=plot.colors['density'],
     lwd=2,
     type='s'
     )

## Draw arrows to illustrate the relationship between density curve and probability
test.p <- dgeom(test.x,p)
arrows(test.x, 0,test.x,test.p,col='#888888',length=0.1,angle=25,code=2,lwd=2)
arrows(test.x, test.p,0,test.p,col='#888888',length=0.1,angle=25,code=2,lwd=2)
text(test.x,0, test.x,pos=4)
text(0,test.p, round(test.p, digits=3),pos=1)

setwd(dir.figures); export.plot(file.prefix='geometric_pmf', export.formats=export.formats.plots, width=7,height=5)

################################################################
## Plot the probability mass function with logarithmic Y scale. Show an
## extended domain of X, to illustrate the capability of log scales to
## display small probabilities.
ext.x <- 0:(n+10)
plot(ext.x - 0.5,
     dgeom(ext.x,p),
     main=paste('probability mass function (log Y scale)'),
     xlab='Waiting time',
     ylab='probability mass function P(X=x) (log scale)',
     ylim=c(min(dgeom(ext.x,p)),1),
     panel.first=grid(equilog=F,col=grid.color,lty='solid'),
     col=plot.colors['density'],
     lwd=2,
     type='s',
     log='y'
     )

## Draw arrows to illustrate the relationship between CDF and left tail probability
log.test.x = c(test.x, 29)
log.test.p <- round(dgeom(log.test.x, p),digits=c(3,5))
arrows(log.test.x, 1e-5,log.test.x,log.test.p,col='#888888',length=0.1,angle=25,code=2,lwd=2)
arrows(log.test.x, log.test.p,0,log.test.p,col='#888888',length=0.1,angle=25,code=2,lwd=2)
text(log.test.x,1e-5, log.test.x,pos=1)
text(0,log.test.p, log.test.p,pos=1)

setwd(dir.figures); export.plot(file.prefix='geometric_pmf_logY', export.formats=export.formats.plots, width=7,height=5)

################################################################
## Plot the probability mass function and highlight the left tail
left.tail.range <- 0:test.x
plot(x - 0.5,
     dgeom(x,p),
     ylim=c(0,1),
     main=paste('Geometric left tail; p=',format(p),', n=', format(n), ', X<=', test.x),
     xlab='Waiting time',
     ylab='probability mass function P(X=x)',
     panel.first=c(abline(v=x.grid,col=grid.color),
       abline(h=y.grid,col=grid.color),
       rect(left.tail.range-0.5, 0, left.tail.range+0.5, dgeom(left.tail.range,p), density = NULL, angle = 45,
            col = '#BBBBBB', border = '#BBBBBB')
       ),
     col=plot.colors['density'],
     lwd=2,
     type='s'
     )
legend('topright',legend=paste('P(X<=',test.x,') = ',round(pgeom(test.x,p),digits=3),sep=''),bty='o', bg='#FFFFFF')
axis(1, at=test.x, labels=test.x,las=1,col.ticks='#888888')

setwd(dir.figures); export.plot(file.prefix='geometric_pmf_left_tail', export.formats=export.formats.plots, width=7,height=5)

################################################################
## Plot the probability mass function and highlight the right tail
right.tail.range <- test.x:n
plot(x - 0.5,
     dgeom(x,p),
     ylim=c(0,1),
     main=paste('Geometric right tail; p=',format(p),', n=', format(n), ', X>=', test.x),
     xlab='Waiting time',
     ylab='probability mass function P(X=x)',
     panel.first=c(abline(v=x.grid,col=grid.color),
       abline(h=y.grid,col=grid.color),
       rect(right.tail.range-0.5, 0, right.tail.range+0.5, dgeom(right.tail.range,p), density = NULL, angle = 45,
            col = '#BBBBBB', border = '#BBBBBB')
       ),
     col=plot.colors['density'],
     lwd=2,
     type='s'
     )
legend('topright',legend=paste('P(X>=',test.x,') = ',round(pgeom(test.x-1,p,lower.tail=F),digits=3),sep=''),bty='o', bg='#FFFFFF')
axis(1, at=test.x, labels=test.x,las=1,col.ticks='#888888')
setwd(dir.figures); export.plot(file.prefix='geometric_pmf_right_tail', export.formats=export.formats.plots, width=7,height=5)

################################################################
## Plot the CDF function
plot(x - 0.5,
     pgeom(x,p),
     ylim=c(0,1),
     main=paste('cumulative distribution function (CDF)'),
     xlab='Waiting time',
     ylab='CDF = P(X<=x)',
     panel.first=c(abline(v=x.grid,col=grid.color),
       abline(h=y.grid,col=grid.color)),
     col=plot.colors['CDF'],
     lwd=2,
     type='s'
     )

## Draw arrows to illustrate the relationship between CDF and left tail probability
left.tail.p <- pgeom(test.x,p)
arrows(test.x, 0,test.x,left.tail.p,col='#888888',length=0.1,angle=25,code=2,lwd=2)
arrows(test.x, left.tail.p,0,left.tail.p,col='#888888',length=0.1,angle=25,code=2,lwd=2)
text(test.x,0, test.x,pos=4)
text(0,left.tail.p, round(left.tail.p, digits=3),pos=3)

setwd(dir.figures); export.plot(file.prefix='geometric_CDF', export.formats=export.formats.plots, width=7,height=5)

################################################################
## Plot the dCDF function
plot(x - 0.5,
     pgeom(x-1,p,lower.tail = F), ## We use x-1 to obtain the inclusive dCDF rather than the default exclusive one
     ylim=c(0,1),
     main=paste('Decreasing cumulative distribution function (dCDF)'),
     xlab='Waiting time',
     ylab='dCDF = P(X>=x)',
     panel.first=c(abline(v=x.grid,col=grid.color),
       abline(h=y.grid,col=grid.color)),
     col=plot.colors['dCDF'],
     lwd=2,
     type='s'
     )

## Draw arrows to illustrate the relationship between dCDF and right tail probability
right.tail.p <- pgeom(test.x-1,p,lower.tail=F)
arrows(test.x, 0,test.x,right.tail.p,col='#888888',length=0.1,angle=25,code=2,lwd=2)
arrows(test.x, right.tail.p,0,right.tail.p,col='#888888',length=0.1,angle=25,code=2,lwd=2)
text(test.x,0, test.x,pos=4)
text(0,right.tail.p, round(right.tail.p, digits=3),pos=3)

setwd(dir.figures); export.plot(file.prefix='geometric_dCDF', export.formats=export.formats.plots, width=7,height=5)


## Show that the sum of CDF and dCDF differ from 1
print(c(dgeom(test.x, p), left.tail.p, right.tail.p, right.tail.p + left.tail.p),digits=3)
