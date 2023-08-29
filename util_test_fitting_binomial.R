################################################################
## Fit a binomial distribution onto a set of observed values.
## + plot an histogram of the observed data + fitted distribution
## + test the goodness of fit
##
## Author: Jacques van Helden

library(stats)
source(file.path(dir.util, 'util_chi2_merged_tails.R'))
       
test.fitting.binomial <- function(x,  # vector of observed values (number of successes)
                                  trials,       # number of trials for each element of the enumeration (=max success value)
                                  proba=NA,     # probabiliy of success at each trial
                                  main="",  # main title for the plot
                                  plot=TRUE,    # plot the result
                                  xlim=NA,      # x limits fr the plot 
                                  ylim=NA,      # y limits for the plot
                                  log="",       # plot with a logarithmic scale
                                  return="auto", # return value ("auto" or "frame")
                                  plot.col=c(obs="gray",fit="blue")
                                  ) {

  #### calculate some parameters on the binomial
  m  <- mean(x)
  if (is.na(proba)) {
    proba <- m/trials
  }
  
  s  <- sd(x)
  s2  <- var(x)
  npq <- trials * proba * (1-proba)

  #### calculate the observed distribution of occurrences
  br  <- seq(-0.5,max(x)+0.5,1)
  histo  <- hist(x,breaks=br,plot=FALSE)
  observed <- histo$counts
  values <- histo$mids
  
  #### calculate expected distribution of occurrences
  expected <- dbinom(values,trials, proba) * length(x)

  ## Test fitting of the observed distribution with the binomial curve
  ## We merge the classes on the left and/or on the right tail of the expected
  ## distribution in order to ensure that sum(exp) = sum(obs)
  exp <- expected
  bino  <- length(x)*dbinom(0:trials,trials, proba)
  exp[1] <- sum(bino[1:(min(values)+1)])
  exp[length(exp)] <- sum(bino[(max(values)+1):trials])

  chi2.result <- chi2.test(observed, exp, min.occ=5, table=T)

  #### plot expected and observed distributions
  if (plot) {
    if (is.na(xlim)) {
      xlim        <- range(br)
    }
    if (is.na(ylim)) {
      ylim        <- c(0,max(c(expected,observed)))
    }
    xmin <- xlim[1]
    xmax <- xlim[2]
    ymin <- ylim[1]
    ymax <- ylim[2]
    plot(values, # to have class mid at occ
         observed,
         type = "h",
         lwd = 4,
         col = plot.col['obs'],
         main = main,
         font.main = 2,
         font.axis = 2,
         font.lab = 2,
         font.leg = 2,
         xlab = "occurrences",
         ylab = "frequency",
         xlim = xlim,
         ylim = ylim,
         log=log
         )

    lines(br,
          c(expected,expected[length(expected)]),
          type = "s",
          lwd = 2,
          col = plot.col['fit']
          )
    
    legend('topright',
           legend=c("observed",
             "fitted",
             paste("chi2=",format(chi2.result$chi2.obs,digits=3)),
             paste("df=",chi2.result$chi2.df),
             paste("Pval=",format(chi2.result$chi2.Pval, scientific=T,digits=2))
             ),
           col = c(plot.col,rep('white',3)),
	   bty="n",
           lty = 1,
           lwd = 2,
           cex = 1
           )
  }
  

  return(chi2.result)


}

#   if (return == "auto") {
#     return(chisq.test(observed,p=bino))
#   } else {
#     diff  <- expected - observed
#     diff2  <- diff^2
#     chisq.vect <- diff2/expected
#     chisq.obs <- sum(chisq.vect)
#     return(data.frame(
#                       values = values,
#                       observed = observed,
#                       expected = expected,
#                       diff = diff,
#                       diff.sq = diff2,
#                       chisq.vect = chisq.vect,
#                       chisq.obs
#                       ))
#   }
#}
