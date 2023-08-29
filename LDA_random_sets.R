################################################################
##
## Random tests with random data.
##
## We generate a random data set and try to train a classifier (LDA)
## to recognize the two sets. The internal validation gives the
## impressio that the eror rate decreases as the number of variables
## increases (up to 100, then the hit rate increases
## again).
##
## However, when the leave-one-out test is performed, one observes
## that the error rate:
## - is much higher than with internal validation (not surprizing, tne
##   internal validation is biased since the object to classify was
##   used for training);
## - increases with the number of variables (over-fitting effect).

library(mda) # for confusion
library(MASS) # for LDA


################################################################
## Function to repet R times the random test and collect the error rates
rand.test <- function(test.data, # a n*p data frame
                        test.labels, # a vector of length n with the group labels
                        r=100, # repetitions of the random test
                        plot.histo=T, # plot the histogram
                        CV=T,
                        ... # additional parameters are transmitted to lda()
                        ) {

  ## Initialize the result object
  result <- list()

  ## Initialize the err.rates data frame
  err.rates <- rep(NA, r)

  ## Perform r random tests
  for (i in 1:r) {
    rand.groups <- sample(test.labels) ## Select random groups
    rand.test <- lda(test.data,rand.groups,CV=CV,...) ## LDA  
    
    ## Calculate confusion table
    if (CV) {
      conf <- confusion(rand.test$class,rand.groups) 
    } else {
      pred <- predict(rand.test,test.data)
      conf <- confusion(pred$class,rand.groups) 
    }
    err.rates[i] <- attr(conf,"error")
    
    print (paste("random test", i, err.rates[i]))
  }
  
  result$err.rates <- err.rates

  ## Plot the histogram of error rates
  if (plot.histo) {
    ## Calculate the class intervals for the histogram
    err.abs <- err.rates*n
    
    if (max(err.abs) > 0) {
      
      x.breaks <- 0:max(err.abs)/n
      
      ## Expected error rate with balanced predictions
      class.occ <- table(test.labels) # number of observations per class
      result$class.occ <- class.occ

      class.freq <- class.occ/n # class frequencies
      result$class.freq <- class.freq
      exp.hit.bal <- class.freq%*%class.freq
#      exp.hit.bal <- sum(class.freq^2)
      exp.err.bal <- 1 - exp.hit.bal
      result$exp.hit.bal <- exp.hit.bal
      result$exp.err.bal <- exp.err.bal
      
      ## Expected error rate with all elements predicted as the most frequent class
      exp.hit.maj <- max(class.freq)
      exp.err.maj <- 1-exp.hit.maj
      result$exp.hit.maj <- exp.hit.maj
      result$exp.err.maj <- exp.err.maj
      
      ## Calculate the limit of the X axis for the drawing
      x.max <- max(x.breaks, exp.err.bal,exp.err.maj)
      
      ## Draw histogram
      hist(err.rates,breaks=x.breaks,xlim=c(0,x.max))
      abline(v=exp.err.bal,col='#FF0000',lwd=2)
      abline(v=exp.err.maj,col='#BB0000',lwd=2)
      
      ## Mean error rate over the R randoms
      print(c(mean(err.rates),exp.err.bal,exp.err.maj))
    } else {
      print("no histogram: all errors are 0")
    }
  }


  result$err.rate.mean <- mean(err.rates)
  result$err.rate.sd <- sd(err.rates)
  return (result)
}

################################################################
## parameters for the test

## Number of repetitions
r <- 10

## Group labels
test.labels <- c(rep("M",30),
                 rep("B",70)
                 )

## Dimensions of the data set
p <- 8 # Default number of variables
n <- length(test.labels) # number of objects

## A vector with increasing numbers of variables, which will progressively lead to over-fitting
ps <- c(2,3,4,5,6,7,8,10,12,16,20,24,32,48,64,80,85,90,95,100,105,110,125,150,200,250,300,400,500) # a series of variable numbers

## perform the test
error.curves <- data.frame(p=ps,
                           err.LOO.mean=rep(NA,length(ps)),
                           err.LOO.sd=rep(NA,length(ps)),
                           err.internal.mean=rep(NA,length(ps)),
                           err.internal.sd=rep(NA,length(ps))
                  )
 
for (p in ps) {
  ## Generate a random data set (normal)
  test.data <- matrix(rnorm(n*p),nrow=n,ncol=p)

  ## LOO validation
  rand.loo <- rand.test(test.data, test.labels, r=r, CV=T)

  ## Update the table for drawing the curves
  error.curves[error.curves$p==p,"err.LOO.mean"] <- rand.loo$err.rate.mean
  error.curves[error.curves$p==p,"err.LOO.sd"] <- rand.loo$err.rate.sd
  
  ## internal validation
  rand.internal <- rand.test(test.data, test.labels, r=r, CV=F,plot.histo=F)

  ## Update the table for drawing the curves
  error.curves[error.curves$p==p,"err.internal.mean"] <- rand.internal$err.rate.mean
  error.curves[error.curves$p==p,"err.internal.sd"] <- rand.internal$err.rate.sd
}

################################################################
## Plot ERROR curves, LOO and internal
x11(width=7,height=7)
plot(error.curves$p,
     error.curves$err.LOO.mean,
     ylim=c(0,1),
     type="b",
     xlab='Number of variables (p)',
     ylab='Error rate',
     main='Number of variables and over-fitting\nrandom tests',
     col='#000088',
     panel.first=c(abline(h=0,col='#000000'),
                  abline(h=1,col='#000000'),
                  grid(col='#000000')),
     pch=16,
     lwd=2
     )
lines(error.curves$p,error.curves$err.internal.mean,type="b",col='#880000',lwd=2, pch=1)

## Expected error rate with balanced predictions
abline(h=rand.loo$exp.err.bal,col='#FF8800',lwd=2, lty='dotted')

## Expected error rate with all elements predicted as the most frequent class
abline(h=rand.loo$exp.err.maj,col='#FF0000',lwd=2, lty='dashed')

## Add a legend
legend("topleft",legend=c("LOO","internal", "expected (balanced classes)", "expected (majority class)"),
       col=c('#000088','#880000', '#FF8800', '#FF0000'),
       lty=c("solid", "solid", "dotted", "dashed"),
       lwd=2,
       bg='white',bty='o')


## Export the plot
setwd(dir.figures); export.plot(file.prefix="LDA_overfitting_randset_error_curve", export.formats=export.formats.plots, width=7,height=7)


################################################################
## Plot HIT curves, LOO and internal
x11(width=7,height=7)
plot(error.curves$p,
     1 - error.curves$err.LOO.mean,
     ylim=c(0,1),
     type="b",
     xlab='Number of variables (p)',
     ylab='Hit rate',
     main='Number of variables and over-fitting\nrandom tests',
     col='#000088',
     panel.first=c(abline(h=0,col='#000000'),
                  abline(h=1,col='#000000'),
                  grid(col='#000000')),
     lwd=2
     )
lines(error.curves$p, 1-error.curves$err.internal.mean,type="b",col='#880000',lwd=2)

## Expected error rate with balanced predictions
abline(h=rand.loo$exp.hit.bal,col='#FF8800',lwd=2, lty='dotted')

## Expected error rate with all elements predicted as the most frequent class
abline(h=rand.loo$exp.hit.maj,col='#FF0000',lwd=2, lty='dashed')

## Add a legend
legend("bottomleft",legend=c("LOO","internal", "expected (balanced classes)", "expected (majority class)"),
       col=c('#000088','#880000', '#FF8800', '#FF0000'),
       lty=c("solid", "solid", "dotted", "dashed"),
       lwd=2,
       bg='white',bty='o')

## Export the plot
setwd(dir.figures); export.plot(file.prefix="LDA_overfitting_randset_hit_curve", export.formats=export.formats.plots, width=7,height=7)
