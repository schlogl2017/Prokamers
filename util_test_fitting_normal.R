## ##############################################################
##
## Fit an enumeration with a normal distribution
## plot the data and fitted curve,
## and test the goodness of fit.
##
## Loading: source(file.path(dir.util, 'util_test_fitting_normal.R'))

library(stats)

test.fitting.normal <- function(x,        # a vector of values
				mean.est=NA,  # impose your own estimate of the mean (if not, it will be computed from input values)
				sd.est=NA,  # impose your own estimate of the standard deviation (if not, it will be computed from input values)
				quartile.est=FALSE, # quartile estimates (more robust)
				plot=TRUE,          # plot the result
				class.interval=NA,  # class interval
                                xlab='values',
                                ylab='frequency',
				xlim=NA,            # x limits for the plot
				ylim=NA,            # y limits for the plot
				tail.group=FALSE,   # plot expected values with tail grouping
				lwd.bars=4,         # plot expected values with tail grouping
				test.ks= F,         # kolmogorov-smirnov fitting test
				test.chi2= F,       # chi2 fitting test
				main="Normal fitting",
                                colors=c(obs="gray",fit="blue"),
                                test.sig=F,          # Test the significance of each value
                                cex.leg=1,
				... # Additional parameters are passsed to the plot function
				) {

  result <- list()
  x <- na.omit(x)
  x <- x[is.finite(x)]  

  ## calculate some parameters on the normal
  n <- length(x)
  m  <- mean(x,na.rm=TRUE)
  s  <- sd(x,na.rm=TRUE)
  s2  <- var(x,na.rm=TRUE)
  min <- min(x,na.rm=TRUE)
  max <- max(x,na.rm=TRUE)
  q <- quantile(x,probs=0.25*1:3,na.rm=TRUE)

  result$n <- n
  result$mean <- m
  result$stdev <- s
  result$min <- min
  result$max <- max
  result$iqr <- IQR(x)
  result$Q <- q
  q.vect <- as.vector(q)

  ## Estimate the mean of the population
  if (is.na(mean.est)) {
    if (quartile.est) {
      mean.est <- q.vect[2]
      mean.est.label <- paste("median=",format(mean.est,digits=3))
    } else {
      mean.est <- m
      mean.est.label <- paste("mean=",format(mean.est,digits=3))
    }
  }
  result$mean.est <- mean.est

  ## Estimate the standard deviation of the population
  if (is.na(sd.est)) {
    if (quartile.est) {
      iqr <- IQR(x, na.rm=T)
      sd.est <- iqr/(qnorm(0.75) - qnorm(0.25))
      sd.est.label <- paste("norm iqr =",format(sd.est,digits=3))
    } else {
      sd.est <- s
      sd.est.label <- paste("sd =",format(sd.est,digits=3))
    }
  }
  result$sd.est <- sd.est
  
  ## Class interval for the histogram
  if (is.na(class.interval)) {
    class.interval <- sd.est/5
  }

  ## calculate the observed distribution of occurrences
  br.min <- floor(min/class.interval)
  br.max <- floor(max/class.interval) + 1
  br  <- class.interval*br.min:br.max
  class.nb <- length(br) -1
  
  histo  <- hist(x,breaks=br,plot=FALSE)
  observed <- histo$counts
  values <- histo$mid
  
  ## calculate expected distribution of occurrences
  pnorm  <- pnorm(br,mean=mean.est,sd=sd.est)
  first <- pnorm[2] # include the left tail of normal proba
  last <- 1 - pnorm[class.nb] # include the right tail of normal proba
  class.max <- pnorm[3:class.nb]
  class.min <- pnorm[2:(class.nb-1)]
  between <- class.max - class.min
  norm <- c(first,between,last)
  expected <- (norm) * n
  
  if (tail.group) {
    expected.plot <- expected
  } else {
    expected.plot <- n * (pnorm[2:(class.nb+1)] - pnorm[1:(class.nb)])
  }

  ## plot expected and observed distributions
  if (plot) {
    if (is.na(xlim)) {
      xlim        <- c(min(x,na.rm=TRUE),max(x,na.rm=TRUE))
    }
    if (is.na(ylim)) {
      ylim        <- c(0,max(c(expected,observed)))
    }
    xmin <- xlim[1]
    xmax <- xlim[2]
    ymin <- ylim[1]
    ymax <- ylim[2]

    plot(values,
         expected.plot,
         type = "l",
         lwd = 2,
         col = colors[2],
         main = main, 
         font.main = 2,
         font.axis = 2,
         font.lab = 2,
         xlab = xlab,
         ylab = ylab,
         xlim = xlim,
         ylim = ylim,
	 ...
         )
    
    lines(values, ## observed is plotted as vertical bars
          observed,
          type = "h",
          lwd = lwd.bars,
          col = colors[1]
          )
    lines(values, ## expected is plotted as lines
          expected.plot,
          type = "l",
          lwd = 2,
          col = colors[2]
          )
    
    legend('topleft',
           c("observed", "fitted"),
           col = colors,
	   bty="n",
           lty = 1,
           lwd = 2,
           cex = cex.leg
           )
    text(xmin,
         ymin + 0.9*(ymax -ymin),
         pos=4,
         font=2,
         paste(mean.est.label,
               sd.est.label,
               sep="\n"),
         )
  }
  
  ## Kolmogorov-Smirnov test
  if (test.ks) {
    result$ks <- ks.test(x,"pnorm",mean=mean.est,sd=sd.est)
  } 

  ## Chi2 fitting test
  if (test.chi2) {
    result$chi2 <- chisq.test(observed,p=norm)

    ## add  the chi2 details  in a table  
    diff  <- expected - observed
    diff2  <- diff^2
    chisq.vect <- diff2/expected
    chisq.obs <- sum(chisq.vect)
    result$table <- data.frame(
			       values		= values,
			       observed	= observed,
			       expected	= expected,
			       diff		= diff,
			       diff.sq		= diff2,
			       chisq.vect	= chisq.vect,
			       chisq.obs
			       )
  } 

  ## Significance test (bilateral)
  if (test.sig) {
    z <-(x - mean.est)/sd.est
    Pval <- pnorm(-abs(z))
    Eval <- Pval * length(z)
    sig <- -log(Eval, base=10)
    result$sig.table <- data.frame(x=x,
                                   z=z,
                                   Pval=Pval,
                                   Eval=Eval,
                                   sig=sig
                                   )
  }
  
  return(result)

}


################################################################
## Draw an histogram of the values and superimpose a normal
## distribution.
plot.normal.fit <- function(values,
                            est="classical", ## supported: "robust", "classical"
                            main=c(est, "normal fit"),
                            palette =c(
                              "obs"="#BBBBBB",
                              "fitted"="black",
                              "mean"="black",
                              "sd"="black",
                              "grid"="#CCCCCC"
                              ),
                            lty=c(
                              "obs"="solid",
                              "fitted"="solid",
                              "sd"="dotted",
                              "mean"="dashed",
                              "grid"="solid"
                              ),
                            cex.leg=1, ## Font size for the legend
                            las=1,
                            breaks=seq(-7, 7, by=0.1), ## breaks between class intervals for the histogram
                            xlab="",
                            ylab="Density",
                            ... ## Additional parameters are passed to plot()
                            ) {

  values.nona <- na.omit(values)
  n <- length(values)
  
  ## Compute robust estimators
  q <- quantile(values.nona,probs=0.25*1:3,na.rm=TRUE)  

  leg.text <- c("observed", "fitted normal")
  
  if (est == "robust") {
    ## Second quartile == median is a robust, non-biaised estimator of the mean
    mean.est <- q[2]  
    
    ## Inter-quartile range is a robust, non-biased estimtator of the standard deviation (with scaling factor)
    sd.est <- sqrt(n/(n-1))*(q[3]-q[1]) / (qnorm(0.75) - qnorm(0.25))

    leg.text <- append(leg.text,
                       c(paste("m.robust=", round(mean.est, digits=2)),
                         paste("sd.robust=", round(sd.est, digits=2))))
  } else {
    ## Second quartile == median is a robust, non-biaised estimator of the mean
    mean.est <- mean(values.nona)
    
    ## Inter-quartile range is a robust, non-biased estimtator of the standard deviation (with scaling factor)
    sd.est <- sd(values.nona)
    leg.text <- append(leg.text,
                       c(paste("mean=", round(mean.est, digits=2)),
                         paste("sd=", round(sd.est, digits=2))))

  }

  h <- hist(values.nona, breaks=breaks,plot=F)

  ## Plot the histogram
  plot(h$mids,
       h$density,
       type="h",
       col=palette["obs"],
       lwd=2,
       main=main,
       xlab=xlab,
       ylab=ylab,
       panel.first=grid(lty=lty["grid"],col=palette["grid"]),
       las=las,
       ...)
  abline(v=(mean.est + sd.est*(-1:1)),col=palette["sd"], lty=lty["sd"])
  abline(v=mean.est,col=palette["mean"], lty=lty["mean"])
  
  ## Overlay the fitted normal curve
  lines(b,
        dnorm(b,mean=mean.est,sd=sd.est),
        type="l",
        col=palette["fitted"],
        lwd=1)

  ## Legend
  legend("topleft", legend=leg.text,
         col=palette[c("obs", "fitted", "mean", "sd")],
         lty=lty[c("obs", "fitted", "mean", "sd")],
         lwd=c(3,1,1,1),
         bty="o", bg="white", cex=cex.leg)

  ## Add a box with the estimator values
  ##    legend("topleft", legend=leg.text, bg="white", bty="o", cex=cex.leg)

  ## Return the results
  result <- list(esimtators=est,
                 mean.est=mean.est,
                 sd.est=sd.est,
                 z=(values - mean.est)/sd.est)

  return(result)
}
