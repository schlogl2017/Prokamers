## ##############################################################
## Fitting a binomial, a Poisson and a normal distribution on an
## observed distribution.
##
## Author: Jacques van Helden
##
## Running this script requires to first run the script config.R
##
## Single-command execution:
## source(file.path(dir.R.files, 'word_count_distrib_GATAA.R'))

## The book is in BW
in.colors <- F
if (in.colors) {
  plot.col <- c(obs='#0000AA',fit='#00AA00',grid='#DDDDDD')
} else { 
  plot.col <- c(obs='#AAAAAA',fit='black',grid='#DDDDDD')
}

## Store the observed distribution in a vector
distrib <- data.frame(occ=0:8,obs=c(223, 337, 247, 119, 54, 13, 6, 0, 1))
occ <- distrib$occ

################################################################
## Calculate parameters for fitting the binomial
L <- 800 ## sequence length
k <- 4 ## word length
pos.per.seq <- L - k + 1 ## Number of possible positions for the word in one sequence
N <- sum(distrib$obs) ## Number of sequences
occ.total <- sum(distrib$occ*distrib$obs) ## Total number of occurrences
occ.per.seq <- occ.total/N ## Average occurrences per sequence
occ.per.pos <- occ.per.seq/pos.per.seq ## Average occurrences per position

################################################################
## Calculate expected occurrences for 1000 sequences
distrib$binomial <- N*dbinom(occ, size=pos.per.seq, prob=occ.per.pos)
## Add the tail of the binomial to the last class (we merge the rightmost classes for the chi2 test)
distrib$binomial[9] <- sum(N*dbinom(8:pos.per.seq, size=pos.per.seq, prob=occ.per.pos))

X11(width=7,height=5)

## Plot the observed distribution
plot(distrib$occ,   # X axis values
     distrib$obs,   # Y axis values
     type='h',      # Plot type: histogram-like
     lwd=4,         # Line width
     col=plot.col['obs'], # Line color: blue
     main='Binomial fit - GATAA occurrences', # graph title
     xlab='Occurrences', # Label for the X axis
     ylab='Number of sequences', # Label for the Y axis
     ylim=c(0,400),  # Limits of the Y axis
     xlim=c(-3, 8),
     panel.first=c(abline(h=50*(0:7),col=plot.col['grid']),
       abline(v=0,col=plot.col['grid']))
    )
## Draw the binomial distribution over the observed distribution
## lines(distrib$occ,
##       distrib$binomial,
##       type='l',
##       lwd=2,
##       col=plot.col['fit']
##       )
rect(distrib$occ-0.5,distrib$binomial,
     distrib$occ+0.5,0,
     lwd=2,
     col=NA,
     border=plot.col['fit']
     )
legend('topright', col=plot.col[c('obs','fit')], lwd=c(4,2), lty=c('solid','solid'),legend=c('observed','fitted binomial'))

setwd(dir.figures); export.plot(file.prefix='word_count_fitting_GATAA_binomial', export.formats=export.formats.plots, width=7,height=5)

################################################################
## Fit the Poisson distribution. For this, we just need the expected
## number of occurrences per sequence
distrib$poisson <- N*dpois(occ, lambda=occ.per.seq)
## Add the tail of the Poisson to the last class (we merge the rightmost classes for the chi2 test)
distrib$poisson[9] <- sum(N*dpois(8:pos.per.seq, lambda=occ.per.seq))

## Draw the Poisson distribution over the observed distribution

## Plot the observed distribution
plot(distrib$occ,   # X axis values
     distrib$obs,   # Y axis values
     type='h',      # Plot type: histogram-like
     lwd=4,         # Line width
     col=plot.col['obs'], # Line color: blue
     main='Poisson fit - GATAA occurrences', # graph title
     xlab='Occurrences', # Label for the X axis
     ylab='Number of sequences', # Label for the Y axis
     ylim=c(0,400),  # Limits of the Y axis
     xlim=c(-3, 8),
     panel.first=c(abline(h=50*(0:7),col=plot.col['grid']),
       abline(v=0,col=plot.col['grid']))
     )
## Draw the binomial distribution over the observed distribution
## lines(distrib$occ,
##       distrib$poisson,
##       type='l',
##       lwd=2,
##       col=plot.col['fit']
##       )
rect(distrib$occ-0.5,distrib$poisson,
     distrib$occ+0.5,0,
     lwd=2,
     col=NA,
     border=plot.col['fit']
      )

arrows(occ.per.seq,380,occ.per.seq,350,col=plot.col['fit'],length=0.1, angle=15, code=2,lwd=3)
text(x=occ.per.seq, y=380,labels=paste('m=', format(occ.per.seq, digits=3),sep=""),col=plot.col['fit'],pos=3,font=2)
legend('topright', col=plot.col[c('obs','fit')], lwd=c(4,2), lty=c('solid','solid'),legend=c('observed','fitted Poisson'))
setwd(dir.figures); export.plot(file.prefix='word_count_fitting_GATAA_Poisson', export.formats=export.formats.plots, width=7,height=5)

## Fit the normal distribution. For this, we just need the expected
## number of occurrences per sequence

## Compute the mean and the standard deviation of the observed distribution
mean.est <- occ.per.seq
sample.var <- sum(distrib$obs*(occ-mean.est)^2)/N
var.est <- sample.var * N/(N-1)
sd.est <- sqrt(var.est)

## Another way to do it (a bit indirect: generate an enumeration)
## enumeration <- rep(distrib$occ,distrib$obs) ## expand the sample to obtain an enumeration
## var.est2 <- var(enumeration)
## sd.est2 <- sd(enumeration) ## Estimate the standard deviation from the sample
distrib$normal.dens <- N*dnorm(occ, mean=mean.est, sd.est)

## Draw the Normal distribution over the observed distribution
plot(distrib$occ,   # X axis values
     distrib$obs,   # Y axis values
     type='h',      # Plot type: histogram-like
     lwd=4,         # Line width
     col=plot.col['obs'], # Line color: blue
     main='Normal density fit - GATAA occurrences', # graph title
     xlab='Occurrences', # Label for the X axis
     ylab='Number of sequences', # Label for the Y axis
     ylim=c(0,400),  # Limits of the Y axis
     xlim=c(-3, 8),
     panel.first=c(abline(h=50*(0:6),col=plot.col['grid']),
       abline(v=0,col=plot.col['grid']))
     )

## The normal is defined in the negative range as well -> we will
## extend the plot of the theoretical curve to highlight this
## difference

## In addition, we show that the normal distribution is continuous,
## not discrete
cont.x <- seq(from=-3, to=8,by=0.02)
cont.norm <- N*dnorm(cont.x, mean=mean.est, sd.est)
lines(cont.x,
      cont.norm,
      type='l',
      lwd=2,
      col=plot.col['fit']
      )

## plot the discrete distribution obtained by normal approximation
ext.x <- -3:8
ext.norm.dens <- N*dnorm(ext.x, mean=mean.est, sd.est)
## lines(ext.x,
##       ext.norm.dens,
##       type='l',
##       lwd=2,
##       lty='dashed',
##       col=plot.col['fit']
##       )
rect(ext.x-0.5,ext.norm.dens,ext.x+0.5,0,
     lwd=2,
     lty='dashed',
     col=NA,
     border=plot.col['fit']
     )

arrows(mean.est,380,mean.est,325,col=plot.col['fit'],length=0.1, angle=15, code=2,lwd=3)
text(x=mean.est, y=380,labels=paste('m=', format(mean.est, digits=3),sep=""),col=plot.col['fit'],pos=3,font=2)
arrows(mean.est,375,mean.est+sd.est,375,col=plot.col['fit'],length=0.1, angle=15, code=3,lwd=2)
text(x=mean.est+sd.est*0.55, y=375,labels=paste('sd=', format(sd.est, digits=3),sep=""),col=plot.col['fit'],pos=1,font=2)
legend('topright', col=plot.col[c('obs','fit','fit')], lwd=c(4,2,2), lty=c('solid','solid','dashed'),legend=c('observed','normal (continuous)','normal (discretized)'))

setwd(dir.figures); export.plot(file.prefix='word_count_fitting_GATAA_normal_density', export.formats=export.formats.plots, width=7,height=5)


## ##############################################################
## Compute the normal area by class interval from x-0.5 to x+0.5, in
## order to treat properly the fact that counts are integer but the
## normal distribution is continuous.
##
## This should only make a small difference with the normal density
## computed above.
##
distrib$normal.area <- N*(pnorm(occ+0.5, mean=mean.est, sd.est)-pnorm(occ-0.5, mean=mean.est, sd.est))

## Draw the Normal distribution over the observed distribution
plot(distrib$occ,   # X axis values
     distrib$obs,   # Y axis values
     type='h',      # Plot type: histogram-like
     lwd=4,         # Line width
     col=plot.col['obs'], # Line color: blue
     main='Normal area fit - GATAA occurrences', # graph title
     xlab='Occurrences', # Label for the X axis
     ylab='Number of sequences', # Label for the Y axis
     ylim=c(0,400),  # Limits of the Y axis
     xlim=c(-3, 8),
     panel.first=c(abline(h=50*(0:6),col=plot.col['grid']),
       abline(v=0,col=plot.col['grid']))
     )

## The normal is defined in the negative range as well -> we will
## extend the plot of the theoretical curve to highlight this
## difference

## In addition, we show that the normal distribution is continuous,
## not discrete
cont.x <- seq(from=-3, to=8,by=0.02)
cont.norm <- N*dnorm(cont.x, mean=mean.est, sd.est)
lines(cont.x,
      cont.norm,
      type='l',
      lwd=2,
      col=plot.col['fit']
      )

## plot the discrete distribution obtained by normal approximation
ext.x <- -3:8
ext.norm.area <- N*(pnorm(ext.x+0.5, mean=mean.est, sd.est)-pnorm(ext.x-0.5, mean=mean.est, sd.est))
## lines(ext.x,
##       ext.norm.area,
##       type='l',
##       lwd=2,
##       lty='dashed',
##       col=plot.col['fit']
##       )
rect(ext.x-0.5,ext.norm.area,ext.x+0.5,0,
     lwd=2,
     lty='dashed',
     col=NA,
     border=plot.col['fit']
     )

arrows(mean.est,380,mean.est,325,col=plot.col['fit'],length=0.1, angle=15, code=2,lwd=3)
text(x=mean.est, y=380,labels=paste('m=', format(mean.est, digits=3),sep=""),col=plot.col['fit'],pos=3,font=2)
arrows(mean.est,375,mean.est+sd.est,375,col=plot.col['fit'],length=0.1, angle=15, code=3,lwd=2)
text(x=mean.est+sd.est*0.55, y=375,labels=paste('sd=', format(sd.est, digits=3),sep=""),col=plot.col['fit'],pos=1,font=2)
legend('topright', col=plot.col[c('obs','fit','fit')], lwd=c(4,2,2), lty=c('solid','solid','dashed'),legend=c('observed','normal (continuous)','normal (discretized)'))

setwd(dir.figures); export.plot(file.prefix='word_count_fitting_GATAA_normal_area', export.formats=export.formats.plots, width=7,height=5)




## ##############################################################
## Estimate the goodness of fit for the 3 fitted distributions
## (binomial, Poisson, and normal).

library(ctest)

## A rough application of the chisq test raises a problem, because
## some classes have exp < 5
chisq.test(distrib$obs,p=distrib$binomial/N)
chisq.test(distrib$obs,p=distrib$poisson/N)
chisq.test(distrib$obs,p=distrib$normal.dens/N)

## ##############################################################
## We will use a custom function which performs the chi2 test with
## some enhancements. 
source(file.path(dir.util, 'util_chi2_merged_tails.R'))

## Run the chi2 test with anb without merging the tails
min.exp <- 5
fitted.distrib <- 'binomial'
for (fitted.distrib in c('binomial', 'poisson', 'normal.dens')) {
  for (min.exp in c(0,5)) {
    print (paste("Chi2 test", fitted.distrib, "min.exp", min.exp))
    chi2.raw <- chi2.test(obs=distrib$obs, exp=distrib[, fitted.distrib], min.exp=NULL, table=T)
    setwd(dir.results); export.object(chi2.raw, file=paste('fitting_word_distrib_GATAA_chi2_', fitted.distrib,'_minexp', min.exp, sep=''), export.formats='print')
  }
}



stop('ending here')


## ##############################################################
## On the basis of the 3 fitted distributions, estimate the
## probability to observe at least 4, 10 and 25 occurrences,
## respectively, of the word GATAA in a 800bp sequence.

## Right tail of the binomial distribution, including the number of
## occurrences

## Prepare a table for storing the result
P.val <- data.frame(
		  occ=0:25,
#		  occ=c(3, 4,10,25),
		  binomial=NA,
		  Poisson=NA,
		  normal=NA)
print(P.val) ## Print the empty result table

for (row in 1:nrow(P.val)) {
  occ <- P.val[row, "occ"]
  P.val[row,"binomial"] <- pbinom(occ-1, size=pos.per.seq, prob=occ.per.pos, lower.tail=FALSE)
  P.val[row,"Poisson"] <- ppois(occ-1,  lambda=occ.per.seq, lower.tail=FALSE)
  P.val[row,"normal"] <- pnorm(occ,mean=occ.per.seq, sd.est, lower.tail=FALSE)
}
print(P.val)  ## Print the filled result table
