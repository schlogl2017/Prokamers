################################################################
#
# Plot a histogram with indication of central tendency
# and dispersion parameters.
#

#### plot a histogram with somt central tendency and dispersion parameters
central.tendency.histo  <-  function(x=rnorm(1000),
                                     main='histogram',
                                     xlab = NA,
                                     ylab = NA,
                                     ymax=NA, ## Height of the histogram
                                     class.number = NA,
                                     class.interval = NA,
                                     show.central=T, ## Display the central tendency parameters
                                     show.dispersion=F, # Display the central dispersion parameters
                                     tick.height = 0.05, ## Height of the ticks for the dispersion bars
                                     stat.colors = c(mean='blue',
                                       median='#00bb00',
                                       mode='#00bbbb',
                                       sd='#0000ee',
                                       iqr='#ff8800',
                                       iqr.norm='#BB6600',
                                       mad='#00DD00',
                                       mad.norm='#008800'
                                       ),
                                     stat.lty = c(mean='solid',
                                       median='solid',
                                       mode='solid',
                                       sd='solid',
                                       iqr='dotted',
                                       iqr.norm='solid',
                                       mad='dotted',
                                       mad.norm='solid'
                                       ),
                                     digits=NULL, ## Number of digits to display on the legends
                                     legend.box=FALSE, ## Draw the legend in a box rather than besides each estimator
                                     side.legend.x = NULL, ## X position of the legend besides each dispersion bar
                                     ... ## Additional parameters are passed to hist()
                                     ) {

  
  #### color palette
  dispersion.bar.height <- c(sd=0.8,
                             iqr.norm=0.65,
                             mad.norm=0.5,
                             mad=0.35,
                             iqr=0.2)

  
  #### class boundaries
  max <- max(x, na.rm=TRUE)
  min <- min(x, na.rm=TRUE)
  if (is.na(class.number)) {
    if (is.na(class.interval)) {
      class.number <- 50
    } else {
      class.number <- as.integer((max - min)/class.interval) + 1
    }
    print(paste("auto class number", class.number))
  }
  if (is.na(class.interval)) {
    breaks <- pretty(x, class.number)
#    print(breaks)
    class.interval <- (max(breaks)-min(breaks))/(length(breaks)-1)
    print(paste("Auto class interval", class.interval))
  } 
  min.class <- as.integer(min/class.interval) - 1
  max.class <- as.integer(max/class.interval) + 1
  
  #### Compue limits for the histogram
  if (is.na(xlab)) { xlab <- paste("class interval = ", class.interval) }
  if (is.na(ylab)) { ylab <- "Occurrences" }
  if (is.na(ymax)) {
    ylim=NULL
  } else {
    ylim=c(0,ymax)
  }

  ## Draw histogram
  histo <- hist(x,
                main=main,
                ylim=ylim,
                breaks = class.interval*min.class:max.class,
                xlab=xlab,
                ylab=ylab, ...
                )

  ## We need ymax for defining the height of the dispersion bars
  if (is.na(ymax)) {
    ymax = max(histo$counts)
  }
  
  ## Calculate central tendency and dispersion parameters
  m <- mean(x,na.rm=TRUE)
  med <- median(x,na.rm=TRUE)
  q1 <- quantile(x,0.25,na.rm=TRUE)
  q3 <- quantile(x,0.75,na.rm=TRUE)
  iqr <- IQR(x, na.rm=T)
  iqr.norm <- IQR(x, na.rm=T)/(qnorm(0.75)-qnorm(0.25))
  mad.norm <- mad(x,na.rm=T)
  mad <- mad(x,constant=1,na.rm=T)
  mode <- histo$mid[which.max(histo$counts)]
#  n <- length(x)
  s <- sqrt(var(x,na.rm=TRUE))  # standard deviation
  legend <- vector()
  legend.colors <- vector()
  legend.lty <- vector()
  
  ## dispersion parameters
  if (show.dispersion) {

    if (is.null(side.legend.x)) {
      side.legend.x <- max(m+s, med+mad.norm, med+iqr.norm) * 1.1
    }
    
    ## IQR
    y <- ymax*dispersion.bar.height['iqr']
    y1 <- y - ymax*tick.height
    y2 <- y + ymax*tick.height
    lines(x=c(q1,q1),y=c(y1,y2),col=stat.colors['iqr'],lty='solid',lwd=1)
    lines(x=c(q3,q3),y=c(y1,y2),col=stat.colors['iqr'],lty='solid',lwd=1)
    lines(c(q1,q3),
          c(y,y),
          col=stat.colors['iqr'],
          lty=stat.lty['iqr'],
          lwd=2
          )
    if (!legend.box) {
      text(x=side.legend.x, y=y, labels=paste("IQR =", format(iqr, digits=digits)), pos=4)
    }

    ## Normalized IQR
    y <- ymax*dispersion.bar.height['iqr.norm']
    y1 <- y - ymax*tick.height
    y2 <- y + ymax*tick.height
    lines(x=c(med,med),y=c(y1,y2),col=stat.colors['iqr.norm'],lty='solid',lwd=1)
    lines(x=c(med+iqr.norm,med+iqr.norm),y=c(y1,y2),col=stat.colors['iqr.norm'],lty='solid',lwd=1)
    lines(c(med,med+iqr.norm),
          c(y,y),
          col=stat.colors['iqr.norm'],
          lty=stat.lty['iqr.norm'],
          lwd=2
          )
    if (!legend.box) {
      text(x=side.legend.x, y=y, labels=paste("IQR.norm =", format(iqr.norm, digits=digits)), pos=4)
    }

    ## MAD
    y <- ymax*dispersion.bar.height['mad']
    y1 <- y - ymax*tick.height
    y2 <- y + ymax*tick.height
    lines(x=c(med,med),y=c(y1,y2),col=stat.colors['mad'],lty='solid',lwd=1)
    lines(x=c(med+mad, med+mad),y=c(y1,y2),col=stat.colors['mad'],lty='solid',lwd=1)
    lines(c(med,med+mad),
          c(y,y),
          col=stat.colors['mad'],
          lty=stat.lty['mad'],
          lwd=2
          )
    if (!legend.box) {
      text(x=side.legend.x, y=y, labels=paste("mad =", format(mad, digits=digits)), pos=4)
    }
    
    ## Normalized MAD
    y <- ymax*dispersion.bar.height['mad.norm']
    y1 <- y - ymax*tick.height
    y2 <- y + ymax*tick.height
    lines(x=c(med,med),y=c(y1,y2),col=stat.colors['mad.norm'],lty='solid',lwd=1)
    lines(x=c(med+mad.norm, med+mad.norm),y=c(y1,y2),col=stat.colors['mad.norm'],lty='solid',lwd=1)
    lines(c(med,med+mad.norm),
          c(y,y),
          col=stat.colors['mad.norm'],
          lty=stat.lty['mad.norm'],
          lwd=2
          )
    if (!legend.box) {
      text(x=side.legend.x, y=y, labels=paste("MAD.norm =", format(mad.norm, digits=digits)), pos=4)
    }
    
    ## Standard deviation
    y <- ymax*dispersion.bar.height['sd']
    y1 <- y - ymax*tick.height
    y2 <- y + ymax*tick.height
    lines(x=c(m,m),y=c(y1,y2),col=stat.colors['sd'],lty='solid',lwd=1)
    lines(x=c(m+s,m+s),y=c(y1,y2),col=stat.colors['sd'],lty='solid',lwd=1)
    lines(c(m,m+s),
          c(y,y),
          col=stat.colors['sd'],
          lty=stat.lty['sd'],
          lwd=2
          )
    if (!legend.box) {
      text(x=side.legend.x, y=y, labels=paste("sd =", format(s, digits=digits)), pos=4)
    }
    
    if (legend.box) {
      legend <- append(legend,
                       c(paste("s =", format(s,digits=digits)),
                         paste("IQR.norm =", format(iqr.norm,digits=digits)),
                         paste("MAD.norm =", format(mad.norm,digits=digits)),
                         paste("MAD =", format(mad,digits=digits)),
                         paste("IQR =", format(iqr,digits=digits))
                         ))

      legend.colors <- append(legend.colors,
                              c(stat.colors['sd'],
                                stat.colors['iqr.norm'],
                                stat.colors['mad.norm'],
                                stat.colors['mad'],
                                stat.colors['iqr']
                                )
                              )
      legend.lty <- append(legend.lty,
                           c(stat.lty['sd'],
                             stat.lty['iqr.norm'],
                             stat.lty['mad.norm'],
                             stat.lty['mad'],
                             stat.lty['iqr']
                             )
                           )
    }
  }

  
  #### central tendency parmeters
  if (show.central) {
    abline(v=mode,col=stat.colors['mode'],lty=stat.lty['mode'],lwd=2)
    abline(v=med,col=stat.colors['median'],lty=stat.lty['median'],lwd=2)
    abline(v=m,col=stat.colors['mean'],lty=stat.lty['mean'],lwd=2)     
    if (legend.box) {
      legend <- append(c(paste("mean =", format(m, digits=digits)),
                         paste("median =",format(med, digits=digits)),
                         paste("mode =", format(mode, digits=digits))
                         ), legend)
      
      legend.colors <- append(c(stat.colors['mean'],
                                stat.colors['median'],
                                stat.colors['mode']
                                ),legend.colors
                              )
      legend.lty <- append(c(stat.lty['mean'],
                             stat.lty['median'],
                             stat.lty['mode']
                             ),legend.lty
                           )
    }
  }

  ## Add a legend
  if (legend.box) {
    legend("topright",
           legend,
           col=legend.colors,
           lty=legend.lty,
           lwd=2)
  }

  ## Generate the result
  result <- list()
  result$mean <- m
  result$mode <- mode
  result$q1 - q1
  result$median <- med
  result$q3 - q3
  result$sd <- s
  result$iqr <- iqr
  result$mad <- mad
  result$mad.norm <- mad.norm
  result$ymax <- ymax
  result$histo <- histo
  return(result)
}


