## ##############################################################
## Normalization of each column of a microarray data set
##
## Developed by Jacques van Helden <Jacques.van-Helden@univ-amu.fr>
## First version: September 2002
## Last modification: December 2003
##
## By default, parameters (mean and sd) are estimated on the basis of
## quartiles. These estimators are robust to outliers.
##
## The result is a list containing different types of information
##     esimators   moment or quartile)
##     m.est       mean estimate for each column
##     s.est       standard dev estimate for each column
##     z           a table with the data converted to z scores
##     P.value     a table with the data converted to P values
##     E.value     a table with the data converted to E values
## Usage
##     	x.n <- normalize.chips(x,quartile.est=T)
## source (file.path(dir.util,'util_chip_standardization.R'))

## Deprecated name for the function
normalize.chips <- function (...) {
  standardize.chips(...)
}

standardize.chips <- function (x, # a data frame containing the expression profiles
                               quartile.est=T, # use mediant and quartiles for parameter estimation 
                               x.desc=NA, # Description of the genes
                               verbosity=1 # level of verbosity
                             ) {

  ## Report starting time
  if (verbosity >= 1) { print (paste(date(), "Standardization  started", sep= " - ")) }

  ## Data structure for storing the result
  result <- list()

  ## Dimensions of the data set
  n <- nrow(x)
  p <- ncol(x)
  
  ## Estimators of central tendency and dispersion
  if (quartile.est) {
    result$estimators <- 'quartiles'
    ## quartile-based estimations
    m.est <- apply(x,2,median,na.rm=T)
    iqr <- apply (x, 2, IQR, na.rm=T)
    s.est <- iqr/(qnorm(0.75) - qnorm(0.25))
  } else {

    ## Classical estimatorss (mean and standard deviation)
    result$estimators <- 'moments'
    m.est <- apply (x, 2, mean, na.rm=T)
    s.est <- apply (x, 2, sd, na.rm=T)
  }
    
  ## prepare space for standardized data
  z <- x
  P.value <- x
  for (i in 1:p) {
    ## Report starting time
    if (verbosity >= 2) { print (paste(date(), "Standardizing column", i, "of", p)) }

    ## calculate z-score
    z[,i] <- (z [,i] - m.est[i])/s.est[i]
    
    ## convert standardized log-ratios to P-values (with Normal assumption)
    P.value[,i] <- pnorm(abs(z[,i]),lower.tail=F)
  }

  ## convert P-values to E-values (taking n into account)
  if (verbosity >= 2) { print (paste(date(), "Calculating E-value")) }
  E.value <- P.value*n
  E.sig <- -log(E.value,base=10)  

  ## calculate max z-score and min E-value and P-value
  if (verbosity >= 2) { print (paste(date(), "Calculating maxima and minima per row")) }
  if (is.na(x.desc)) {
    extr <- data.frame(row.names=row.names(x), gene=row.names(x))
  } else {
    extr <- data.frame(row.names=row.names(x), gene=x.desc$NAME)
  }
  extr$z.min <- apply(z,1,min,na.rm=T)
  extr$z.max <- apply(z,1,max,na.rm=T)
  extr$z.abs.max <- apply(abs(z),1,max,na.rm=T)
  extr$P.value.min <- apply(P.value,1,min,na.rm=T)
  extr$E.value.min <- apply(E.value,1,min,na.rm=T)
  extr$E.sig.max <- apply(E.sig,1,max,na.rm=T)

  #### prepare and return the result
  if (verbosity >= 2) { print (paste(date(), "Building the result object")) }
  result$m.est <- m.est
  result$s.est <- s.est
  result$z <- z
  result$P.value <- P.value
  result$E.value <- E.value
  result$E.sig <- E.sig
  result$extr <- extr

  if (verbosity >= 1) { print (paste(date(), "Standardization  done", sep= " - ")) }
  
    return(result)
  
}

## ##############################################################
## demo for the utilization of standardize.chips()
standardize.chips.demo <- function () {
  ## generate a random data set (6000 rows, 10 columns)
  x <- rnorm(60000,mean=0.012,sd=0.43)
  x <- data.frame(matrix(x,nrow=6000,ncol=10))
  print ("data set")
  print (dim(x))
  
  ## we will test both quartile estimates and moment estimates
  for (q in c(T,F)) {
    print ("standardizing")

    ## standardize the data
    x.n <- standardize.chips(x,quartile.est=q)

    print (names(x.n))
    
    attach(x.n)
    print (c("estimators", estimators))
    print ("mean estimates")
    print (m.est)
    print ("sdandard deviation estimates")
    print (s.est)
    
    ## calculate mean of the z scores
    print ("column means after standardization")
    print (apply(z,2,mean))
    
    ## calculate the stndard deviation of the z-scores
    print ("column standard deviations after standardization")
    print(apply(z,2,sd))
    
    detach(x.n)
  }
}


################################################################
## Old name maintained for backward compatibility
select.act.repr <- function(...) {
  select.up.down(...)
}

################################################################
## Select lists of significantly up- and down-regulated genes for
## each column of a set of microarray data
## Input: the input is an object generated by the function
##    standardize.chips()
## Ouput
##    The output is a cluster table,
##    - the first column indicates the element (gene),
##    - the second indicates the cluster name (the condition, i.e. the column in the original expression table)
##    - the third column indicates the z-score for the selected gene in the selected condition
##    - the fourth column indicates the significance of the association between the gene and the condition)
select.up.down <- function (x.stand, ## Must be the result of standardize.chips
                            E.threshold = 0 ## Upper threshold on E-value per column
                            ) {
  result <- data.frame()
  for (c in 1:ncol(x.stand$z)) {
    if (verbosity >= 2) { print(paste("Selecing significant elements in column ", c)) }
    
    up <- !is.na(x.stand$E.sig[,c]) & (x.stand$E.value[,c] <= E.threshold) & (x.stand$z[,c] > 0)
    if (sum(up) > 0) {
      result <- rbind(result,
                      data.frame(id=row.names(x.stand$z)[up],
                                 cluster=paste(names(x.stand$z)[c], "up", sep="_"),
                                 z=x.stand$z[!is.na(up) & up,c],
                                 sig=x.stand$E.sig[!is.na(up) & up,c])
                      )
    }
    
    down <- !is.na(x.stand$E.sig[,c]) & (x.stand$E.value[,c] <= E.threshold) & (x.stand$z[,c] < 0)
    if (sum(down) > 0) {
      result <- rbind(result,
                      data.frame(id=row.names(x.stand$z)[down],
                                 cluster=paste(names(x.stand$z)[c], "down", sep="_"),
                                 z=x.stand$z[!is.na(down) & down,c],
                                 sig=x.stand$E.sig[!is.na(down) & down,c])
                      )
    }
  }
  return(result)
}


################################################################
## Apply the quantile normalization method as described in the affy
## manuel for the rma() method.
quantile.norm <- function (x, ## data frame
                           sample.to.plot=NULL, ## generate a plot with the pre- and post-normalization data for a selected sample number
                           prefix="data"
                           ) {
  

  ## Compute median per row for the column-sorted intensity values
  x.sorted.columns <- x
  for (col in 1:ncol(x)) {
    x.sorted.columns[,col] <- sort(x[,col], decreasing=T)
  }
  medians <- apply(x.sorted.columns, 1, median)

  ## Compute ranks for each sample
  x.ranks <- x
  for (col in 1:ncol(x)) {
    x.ranks[,col] <- nrow(x) - rank(x[,col])+1
  }
  

  ## Assign new values to all intensities, according to their rank
  x.quantile.norm <- x
  for (col in 1:ncol(x)) {
    x.quantile.norm[,col] <- medians[nrow(x) - rank(x[,col])+1]
  }

  
  ## Plot one sample if requested (option sample.to.plot)
  if (!is.null(sample.to.plot)) {
    x11(width=9,height=9)
    to.plot <- data.frame(x[,sample.to.plot],
                          x.quantile.norm[,sample.to.plot],
                          x.ranks[, sample.to.plot]
                          )
    colnames(to.plot) <- c(prefix,
                           paste(prefix, ".quant", sep=""),
                           paste(prefix, ".rank", sep=""))
    plot(to.plot, col='#BBBBBB')
  }
  
#  head(quantile.norm.check)
#  median.per.row <- apply(x, 1, median)
#  hist(median.per.probeset, breaks=seq(min(x),max(x), by=0.1), xlim=c(0,16))
  return (x.quantile.norm)
}
