################################################################
#
# Fit a set of observed values with a Poisson distribution plot the
# data and fitted curves, and test the goodness of fit with a
# chi-squared test

## source(file.path(dir.R.files, 'util/util_test_fitting_poisson.R'))


library(stats)
source(file.path(dir.util, 'util_chi2_merged_tails.R'))


## ##############################################################
## The values are provided as an enumeration
test.fitting.poisson <- function(x,  # a vector of observed values
				 lambda=NA,    # expected mean
                                 breaks=NA, ## Support user-specific breaks
				 ...
				 ) {


  #### calculate some parameters on the Poisson
  m  <- mean(x)
  if (is.na(lambda)) {
    lambda <- m
  }

  s  <- sd(x)
  s2  <- var(x)
  
  ## Define the breaks for the distributions
  ## By default, one break per natural value (given the discrete nature of the Poisson distribution)
  ## However, in some cases it might be good to specify a number of breaks, 
  ## for example if the poissson has a very large lambda (Poisson tending towards a normal).
  if (is.na(breaks)) {
    if (max(x) > 50) {
      br  <- pretty(x, n=25)
    } else {
      br  <- seq(-0.5,max(x)+0.5,1)
    }
  } else {
    if (length(breaks)==1) { 
      br  <- pretty(x, n=breaks)
    } else {
      br <- breaks
    }
  }

  
  ## calculate the observed distribution of occurrences
  histo  <- hist(x,breaks=br,plot=FALSE)
  counts <- histo$counts
  values <- histo$mid


  
  result <- test.fitting.poisson.histo(values=values, counts=counts, lambda=lambda, obs.mean=mean(x), histo, ...)
  
  return (result)    
}

## ##############################################################
## Data is provided as values and counts
test.fitting.poisson.histo <- function(values,  # a vector of observed values
				       counts, # a vector of counts, same length as values 
                                       histo=NULL, ## The breaks of the histogram
                                       obs.mean=NA, ## observed mean
				       lambda=NA,    # expected mean
				       subtitle="",  # subtitle for the plot
				       plot=TRUE,    # plot the result
				       xlim=NA,      # x limits fr the plot 
				       ylim=NA,      # y limits for the plot
				       log=""       # plot with a logarithmic scale
				       ) {

  values <- as.vector(as.matrix(values))
  counts <- as.vector(as.matrix(counts))
  
  ## Observed mean
  if (is.na(obs.mean)) {
    obs.mean <- sum(values*counts)/ sum(counts)
  }
  
  if (is.null(histo)) {
    breaks <- values
    x.values <- values
  } else {

    ## define all the integer values in the range of the user-specified breaks
    breaks <- histo$breaks
#    x.values <- round(histo$mids)
    x.values <- min(breaks):max(breaks)
  }
  
  
  ## lambda
  if (is.na(lambda)) {
    lambda <- obs.mean
  }

  ## calculate expected distribution of occurrences
  pois <- dpois(x.values,lambda)
  expected <- pois * sum(counts)
  ## plot(x.values, expected, panel.first=grid())


  ## regroup expected by class intervals if required
  if (!is.null(histo)) {
    exp.cumsum <- cumsum(expected)
    ## plot(x.values, exp.cumsum, panel.first=grid())
    
    class.interval <- breaks[2]-breaks[1]
    break.nb <- length(breaks) ## Number of breaks
    class.nb <- break.nb -1
    
    ## compute the min and max values for each class (natural numbers)
    class.min <- seq(from=breaks[1], by=class.interval, length.out=class.nb)
    class.max <- class.min + class.interval -1
    
    ## Expected number of observations by class
    exp.by.class <- exp.cumsum[x.values %in% class.max]
    exp.by.class[2:length(exp.by.class)] <- exp.by.class[2:length(exp.by.class)] - exp.by.class[1:(length(exp.by.class)-1)]
    expected <- exp.by.class
    
    
    ## Normalize expected to obtain the same sum as observed (avoid problems with rounding of class means)
    expected <- expected * sum(counts) / sum(expected)
  }

  print(data.frame(x.values, pois, expected))
  stop("HELLO")

  ## calculate expected distribution of occurrences
#  pois <- dpois(x.values,lambda)

  ## Add the tails of non-observed values to the expected
  if (min(x.values) > 0) {
    pois[1] = pois[1] + ppois(min(x.values)-1, lambda)
  }
  if (sum(pois) < 1) {
    pois[length(pois)] = pois[length(pois)]  + 1 - sum(pois)
  }

  expected <- pois * sum(counts)
#  stop(paste("HELLO", sum(expected)))
    

  ## plot(x.values, expected, type="h", panel.first=grid())
  
  ## Run the chi2 test, with specific options to regroup tails if
  ## required, in order to ensure the applicability condition (exp
  ## freq >= 5 for all classes)
  result <- chi2.test(obs=counts,exp=expected, min.exp=5)
  pval=result$chi2.Pval
  
  ## plot expected and observed distributions
  if (plot) {
    if (is.na(xlim)) {
      xlim        <- c(min(x.values),max(x.values))
    }
    if (is.na(ylim)) {
      ylim        <- c(0,max(max(expected),max(counts)))
    }
    xmin <- xlim[1]
    xmax <- xlim[2]
    ymin <- ylim[1]
    ymax <- ylim[2]
    
    
    colors <- c(exp="#000088",
                obs="#888888",
                exp2="#FF0088")
    
    ## Plot the observed distribution
    plot(values, counts, type="h", lwd=4, col=colors["obs"],
         main = paste("Poisson fitting ",subtitle),
  #       font.main = 4,
  #       font.axis = 4,
  #       font.lab = 4,
         xlab = "occurrences",
         ylab = "frequency",
         xlim = xlim,
         ylim = ylim,
         log=log
    )
    
    ## Plot the fitted Poisson
    lines(values, expected, type="l", lwd=2, col=colors["exp"])    

    if (max(values)/lambda < 2) {
      legend.location <- "topleft"
    } else {
      legend.location <- "topright"
    }
    legend(legend.location,
           c("expected","observed",
             paste(sep="", "lambda=", format(lambda, digits=4)),
             paste(sep="", "obs.mean=", format(m, digits=4)),
             paste(sep="", "pval=", format(pval, scientific=TRUE, digits=4))),
           col = colors,
           bty="n",
           lty = 1,
           lwd = c(2,4,0,0,0),
           cex = 1
           )
    
  }
  
  return(result)
}




