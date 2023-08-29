################################################################
##
## Draw frequency polygon with class frequencies
##
## Author: Jacques van Helden
##
## Running this script requires to first run the script config.R
##
## Single-command execution:
##   source(file.path(dir.R.files, 'freq_polygon_orf_lengths.R'))

#### read the data
yeast.gene.len	<- scan(file.path(dir.data, 'orf_lengths/yeast_orf_length_enum.txt'))
ecoli.gene.len	<- scan(file.path(dir.data, 'orf_lengths/ecoli_orf_length_enum.txt'))

#### calculate number of classes
max.value <- max(c(yeast.gene.len),c(ecoli.gene.len))
max.value.plot <- 5000 ## Skip the right tail of the distribution, due to an outlier
class.interval <- 100
class.number <- as.integer(max.value/class.interval) + 1
class.breaks <- class.interval*0:class.number
class.min <- class.breaks[1:class.number]
class.max <- class.breaks[2:(class.number+1)]
class.mids <- (class.min + class.max)/2

#### Use the hist() function to compute yeast class frequencies without plotting the histo
yeast.histo <- hist(yeast.gene.len,
                    breaks = class.breaks,
                    plot=F
                    )

#### Use the hist() function to compute ecoli class frequencies  without plotting the histo
ecoli.histo <- hist(ecoli.gene.len,
                    breaks = class.breaks,
                    plot=F
                    )


#### calculate relative frequencies
yeast.freq <- yeast.histo$counts/sum(yeast.histo$counts)
ecoli.freq <- ecoli.histo$counts/sum(ecoli.histo$counts)
max.freq <- max (c(yeast.freq,ecoli.freq))

## Define graphical parameters
colors <- c(yeast="#888888",
            ecoli="#000000")

lty <- c(yeast="solid",
         ecoli="solid")

################################################################
## Draw frequency polygons

#par(cex.main=2)
#par(cex.lab=1.5)
X11(width=7,height=5)

plot(yeast.histo$mid,
     yeast.freq,
     type="l",
     lty=lty["yeast"],
     col=colors["yeast"],
     main='Frequency distributions of gene lengths',
     xlab=paste("gene length (class interval = ", class.interval, ")", sep=""),
     ylab="Frequency",
     lwd=2,
     xlim = range(0,max.value.plot),
     ylim= range(0,max.freq),
     panel.first=c(## grid
       abline(h=0,col='black'),
       abline(v=1000*(0:16),col='#BBBBBB'),
       abline(h=seq(from=0,to=0.1,by=0.01),col='#BBBBBB'))
     )

lines(ecoli.histo$mid,
      ecoli.freq,
      type="l",
      lty=lty["ecoli"],
      col=colors["ecoli"],
      lwd=2
      )

legend("topright", legend=c("Saccharomyces cerevisiae", "Escherichia coli K12"),lty=lty, col=colors, bty="o",bg="white",lwd=3)

## Export the plot
setwd(dir.figures); export.plot(file.prefix="freq_polygon_gene_lengths", export.formats=export.formats.plots, width=7,height=5)

## ##############################################################
## plot the cumulative frequency polygon
yeast.cum.freq <- cumsum(yeast.freq)
plot(class.max, ## For cumulative frequencies, the horizontal axis displays the class max
     yeast.cum.freq,
     type="l",
     lty=lty["yeast"],
     col=colors["yeast"],
     main='Cumulative frequency distributions of gene lengths',
     xlab=paste("Gene length (class interval = ", class.interval, ")", sep=""),
     ylab="Cumulative frequency",
     lwd=2,
     xlim = range(0,max.value.plot),
     ylim= range(0,1),
     panel.first=c(## grid
       abline(h=0.1*(0:10),col='#DDDDDD'),
       abline(v=1000*(0:16),col='#DDDDDD'),
       abline(h=c(0,1),col='black'),
       abline(v=0,col='black'))
     )
ecoli.cum.freq <- cumsum(ecoli.freq)
lines(class.min,
      ecoli.cum.freq,
      type="l",
      lty=lty["ecoli"],
      col=colors["ecoli"],
      lwd=2
      )
legend("bottomright", legend=c("Saccharomyces cerevisiae", "Escherichia coli K12"),lty=lty, col=colors, bty="o",bg="white",lwd=3)

################################################################
## Estimate the frequency below a given quantile and display the
## principle of its identification.
quantile <- 0.3
above.quantile.index <- which(yeast.cum.freq >= quantile)[1]
below.quantile.index <- above.quantile.index -1
quantile.est <- class.max[below.quantile.index] + class.interval * (quantile - yeast.cum.freq[below.quantile.index]) / (yeast.cum.freq[above.quantile.index] - yeast.cum.freq[below.quantile.index])

## Draw arrows to illustrate the principle of quantile idenfication on the graph
arrows(0, quantile,quantile.est,quantile,col=colors['yeast'],length=0.1,angle=25,code=2,lwd=1)
arrows(quantile.est, quantile,quantile.est,0,col=colors['yeast'],length=0.1,angle=25,code=2,lwd=1)
axis(2 , at=quantile, labels=quantile,las=2,col.ticks=colors['yeast'])
axis(1 , at=quantile.est, labels=round(quantile.est),las=2,col.ticks=colors['yeast'])

## Draw arrows showing the way from values (abcsissa) to cumulated frequencies (ordinate)
value <- 1700
above.value.index <- which(class.max>value)[1]
below.value.index <- above.value.index - 1
above.value.freq <- yeast.cum.freq[above.value.index]
below.value.freq <- yeast.cum.freq[below.value.index]
value.freq <- below.value.freq + (above.value.freq - below.value.freq) * (value - class.max[below.value.index]) / (class.max[above.value.index] - class.max[below.value.index])

arrows(value, 0, value, below.value.freq,col=colors['yeast'],length=0.1,angle=25,code=2,lwd=1)
arrows(value, below.value.freq, 0,below.value.freq,col=colors['yeast'],length=0.1,angle=25,code=2,lwd=1)
axis(1 , at=value, labels=value,las=2,col.ticks=colors['yeast'])
axis(2 , at=below.value.freq, labels=round(below.value.freq, digits=2),las=2,col.ticks=colors['yeast'])

## Export the plot
setwd(dir.figures); export.plot(file.prefix="freq_polygon_gene_lengths_cum", export.formats=export.formats.plots, width=7,height=5)

## ##############################################################
## Plot decreasing cumulative frequency polygon
yeast.decr.cum.freq <- rev(cumsum(rev(yeast.freq)))
plot(class.min, ## For decreasing cumulative frequencies, the X axis shows the class min
     yeast.decr.cum.freq,
     type="l",
     lty=lty["yeast"],
     col=colors["yeast"],
     main='Decreasing cumulative frequency distributions of gene lengths',
     xlab=paste("Gene length (class interval = ", class.interval, ")", sep=""),
     ylab="Decreasing cumulative frequency",
     lwd=2,
     xlim = range(0,max.value.plot),
     ylim= range(0,1),
     panel.first=c(## grid
       abline(h=0.1*(0:10),col='#DDDDDD'),
       abline(v=1000*(0:16),col='#DDDDDD'),
       abline(h=c(0,1),col='black'),
       abline(v=0,col='black'))
     )
ecoli.decr.cum.freq <- rev(cumsum(rev(ecoli.freq)))
lines(class.min,
      ecoli.decr.cum.freq,
      type="l",
      lty=lty["ecoli"],
      col=colors["ecoli"],
      lwd=2
      )


legend("topright", legend=c("Saccharomyces cerevisiae", "Escherichia coli K12"),lty=lty, col=colors, bty="o",bg="white",lwd=3)

################################################################
## Estimate the frequency below a given quantile and display the
## principle of its identification.
quantile <- 0.12

## Yeast quantile
above.quantile.index <- which(yeast.decr.cum.freq <= quantile)[1]
below.quantile.index <- above.quantile.index -1
quantile.est <- class.min[below.quantile.index] + class.interval * (quantile - yeast.decr.cum.freq[below.quantile.index]) / (yeast.decr.cum.freq[above.quantile.index] - yeast.decr.cum.freq[below.quantile.index])
arrows(0, quantile,quantile.est,quantile,col=colors['yeast'],length=0.1,angle=25,code=2,lwd=1)
arrows(quantile.est, quantile,quantile.est,0,col=colors['yeast'],length=0.1,angle=25,code=2,lwd=1)
axis(2 , at=quantile, labels=quantile,las=2,col.ticks=colors['yeast'])
axis(1 , at=quantile.est, labels=round(quantile.est),las=2,col.ticks=colors['yeast'])

## Ecoli quantile
above.quantile.index <- which(ecoli.decr.cum.freq <= quantile)[1]
below.quantile.index <- above.quantile.index -1
quantile.est <- class.min[below.quantile.index] + class.interval * (quantile - ecoli.decr.cum.freq[below.quantile.index]) / (ecoli.decr.cum.freq[above.quantile.index] - ecoli.decr.cum.freq[below.quantile.index])
arrows(0, quantile,quantile.est,quantile,col=colors['ecoli'],length=0.1,angle=25,code=2,lwd=1)
arrows(quantile.est, quantile,quantile.est,0,col=colors['ecoli'],length=0.1,angle=25,code=2,lwd=1)
axis(2 , at=quantile, labels=quantile,las=2,col.ticks=colors['ecoli'])
axis(1 , at=quantile.est, labels=round(quantile.est),las=2,col.ticks=colors['ecoli'])

## Export the plot
setwd(dir.figures); export.plot(file.prefix="freq_polygon_gene_lengths_inv_cum", export.formats=export.formats.plots, width=7,height=5)

################################################################
## Summarize the data in a data frame
yeast.distrib <- data.frame(min=class.min,
                            mid=yeast.histo$mid,
                            max=class.max,
                            count=yeast.histo$count,
                            cumcount=cumsum(yeast.histo$count),
                            invcumcount=rev(cumsum(rev(yeast.histo$count))),
                            freq=yeast.freq,
                            cumfreq = cumsum(yeast.freq),
                            invcum = rev(cumsum(rev(yeast.freq)))
                            )
ecoli.distrib <- data.frame(min=class.min,
                            mid=ecoli.histo$mid,
                            max=class.max,
                            count=ecoli.histo$count,
                            cumcount=cumsum(ecoli.histo$count),
                            invcumcount=rev(cumsum(rev(ecoli.histo$count))),
                            freq=ecoli.freq,
                            cumfreq = cumsum(ecoli.freq),
                            invcum = rev(cumsum(rev(ecoli.freq)))
                            )
