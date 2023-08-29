#' @title Chi2 test of homongeneity between vectors of Natural numbers
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Taking as input two vectors of Natural numbers, compute the
#' distributions (frequencies of each number), and  apply a chi-squared test
#' of homogeneity between the two distributions. The chi-squared test is computed 
#' from the scractch and compared to the result obtained with the R function
#' stats::chisq.test().
#'
#' @details
#' First version: ?. 
#' Last modification: 2015-10. 
#'
#' @param x1   First vector of values (must be Natural numbers)
#' @param x2   First vector of values (must be Natural numbers)
#' @param plot="none" If non-null, plot the two histograms in either semi-transparent 
#' histograms (plot="hist“) or polygon frequencies (plot="lines").
#' @param ... Additional parameters are passed to plot().
#' 
#' @examples
#'
#' ## In this example, we generate two random distributions of numbers 
#' ## following Poisson distributions, and test their fitting.
#' x1 <- rpois(n=1000, lambda=3.2)
#' x2 <- rpois(n=1000, lambda=3.2)
#' result <- chisq.homogeneity.test(x1,x2, plot="lines")
#' result <- chisq.homogeneity.test(x1,x2, plot="hist")
#' 
#' ## The next example performs the fitting between two Poisson distributions 
#' ## with different lambda values
#' x1 <- rpois(n=1000, lambda=3)
#' x2 <- rpois(n=1000, lambda=3.5)
#' result <- chisq.homogeneity.test(x1,x2, plot="hist")
#'
#' @export
chisq.homogeneity.test <- function(x1, ## Vector of values from the first sample
                                   x2,  ## Vector of values from the second sample
                                   plot="none", ## Plot the distributions
                                   ... ## Additional parameters are passed to plot
                                   ) {
  max.len <- max(c(x1,x2))
  breaks=(0:(max.len+1))-0.5
  
  n1 <- length(x1)
  hist1 <- hist(x1, breaks=breaks,plot=F)

  n2 <- length(x2)
  hist2 <- hist(x2, breaks=breaks,plot=F)

  chi2.table <- data.frame(x=hist1$mids, f1=hist1$counts,f2=hist2$counts)
  chi2.table$x.sum <- chi2.table$f1+chi2.table$f2
  chi2.sum <- apply(chi2.table, 2,sum)
  chi2.table$exp1 <- chi2.table$x.sum*n1/(n1+n2)
  chi2.table$exp2 <- chi2.table$x.sum*n2/(n1+n2)
  chi2.sums <- apply(chi2.table, 2,sum)

  attach(chi2.table)
  chi2.obs <- sum((f1-exp1)^2/exp1+(f2-exp2)^2/exp2)
  df <- nrow(chi2.table)-1
  p.val <- pchisq(q=chi2.obs,df=df,lower.tail=F)
  detach(chi2.table)


  if (plot == "hist") {
    plot(hist1,col=rgb(0,0,1, alpha=0.25),xlim=range(breaks), 
         xlab="Values", ylab="Frequencies", ...)
    plot(hist2,col=rgb(1,0.5,0, alpha=0.25),xlim=range(breaks), add=TRUE, ...)
    legend("topright", 
           legend=c(
             paste("chi2 =", signif(digits=3, chi2.obs)), 
             paste("pval =", signif(digits=3, p.val))))
  } else if (plot == "lines") {
    plot(chi2.table$x, chi2.table$f1, type="b", col="blue", 
         panel.first=abline(v=0:max.len, col="gray"),
         xlab="Values", ylab="Frequencies")
    lines(chi2.table$x, chi2.table$f2, type="b", col="darkgreen")
    legend("topright", 
           legend=c(
             paste("chi2 =", signif(digits=3, chi2.obs)), 
             paste("pval =", signif(digits=3, p.val))))
  }
  
  result <- list()
  result$n1 <- n1
  result$n2 <- n2
  result$chi2.table <- chi2.table
  result$chi2.obs <- chi2.obs
  result$df <- df
  result$p.val <- p.val
  
  ## Compare the manually computed result with the result from chisq.test function
  result$chisq.test.result <- chisq.test(chi2.table[,c("f1","f2")])
  
  return(result)
}

