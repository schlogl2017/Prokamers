## ##############################################################
##
## Illustration of the sampling of the mean
##
## Author: Jacques van Helden
##
## Running this script requires to first run the script config.R
##
## Single-command execution:
## source(file.path(dir.R.files, 'sampling_mean_distribution.R'))


#### take a sample of size n from a uniform distribution 
#### and measure the mean
#### repeat r times
#### plot the distribution of the mean

plot.col <- c(histo='#888888',
              borders='#888888',
              points='black',
              grid='#CCCCCC',
              lines='#666666')

n.values <- c(1,2,4,16,64,128) # values of sample size for the histograms of mean distribution
#n.values <- c(1,2,4)
max.n <-max(n.values)
repetitions <- 10000         # random repetitions

################################################################
## Random sampling and mean computation

## Prepare a big table with all the sampled values
#values <- matrix(nrow=repetitions, ncol=max.n, data=runif(max.n*repetitions))
values <- matrix(nrow=repetitions, ncol=max.n, data=sample(1:6,max.n*repetitions,replace=T))

## Compute the sample mean and sample var for each repetition (row),
## and for each sample size (cumulative columns)
sample.means <- matrix(nrow=repetitions, ncol=max.n, data=NA)
sample.vars <- matrix(nrow=repetitions, ncol=max.n, data=NA)
for (r in 1:repetitions) {
  sample.means[r,] <- cumsum(values[r,])/(1:max.n)
  sample.vars[r,] <- cumsum(values[r,]^2)/(1:max.n) - sample.means[r,]^2
}


## Compute parameters of the sampling mean distribution
sampling.stats <- data.frame(sample.size=1:max.n,
                             mean.of.means=apply(sample.means,2,mean),
                             var.of.means=apply(sample.means,2,var))
sampling.stats$sd.of.means <- sqrt(sampling.stats$var.of.means)

## Compute the theoretical variance For a uniform discrete
## distribution with m possible values, the theoretical variance of
## the population is v=(m^2-1)/12.
##
## The variance of the sample mean is v^2/n
sampling.stats$var.of.means.theor <- (6^2-1)/(12*(1:max.n))
sampling.stats$sd.of.means.theor <- sqrt(sampling.stats$var.of.means.theor)

################################################################
#### plot histograms of sample mean distributions
X11(width=12,height=8)
par(font.lab=2)
par(font.axis=2)
par(mfrow=c(2,3))
par(cex=1)
par(cex.lab=1)
par(mai=c(0.9,0.9,0.2,0.2))

## Plot the distribution of sample means only for selected repetitions
for (n in n.values) {
  means <- sample.means[,n]
  h <- hist(means,
       freq=T,
       breaks=50, plot=F)
  plot(h,
       col=plot.col['histo'],
       border=plot.col['border'],
       main=NULL,
       xlim=c(1,6),
       ylim=c(0,1.4*max(h$counts)),
       xlab='sample mean (m)')
  text(1,1.2*max(h$counts),
       paste('n=',n,'\nmean(m) = ', format(sampling.stats[n,'mean.of.means'],dig=3),
             '\nsd(m) = ', format(sampling.stats[n,'sd.of.means'],dig=2)),font=2,pos=4)
}
setwd(dir.figures); export.plot(file.prefix=paste("sampling_mean_histograms_r", repetitions, sep=''), export.formats=export.formats.plots, width=12,height=8)


## Plot the distribution of sample variances only for selected repetitions
n.values.v <- c(2,3,4,8,16,128) # values of sample size for the histograms of variance distribution
for (n in n.values.v) {
  vars <- sample.vars[,n]
  h <- hist(vars,
       freq=T,
       breaks=50, plot=F)
  plot(h,
       col=plot.col['histo'],
       border=plot.col['border'],
       main=NULL,
       xlim=c(1,6),
       ylim=c(0,1.4*max(h$counts)),
       xlab='sample var (m)')
  text(1,1.2*max(h$counts),
       paste('n=',n,'\nvar(m) = ', format(sampling.stats[n,'var.of.vars'],dig=3),
             '\nsd(m) = ', format(sampling.stats[n,'sd.of.vars'],dig=2)),font=2,pos=4)
}
setwd(dir.figures); export.plot(file.prefix=paste("sampling_var_histograms_r", repetitions, sep=''), export.formats=export.formats.plots, width=12,height=8)

par(mfrow=c(1,1))


################################################################
#### plot the mean,var and sd as a function of n 
x11(width=7,height=5)
par(font.lab=2)
par(font.axis=2)

plot(1:max.n,
     sampling.stats$mean.of.means,
#     ylim=c(0.505,0.495),
     xlab='n (sample size)',
     ylab='mean of the means',
     main='mean of the sample means',
     col=plot.col['points'],
     type='p',
     pch=1,
     panel.first=c(grid(col=plot.col['grid'],lty='solid'),
       abline(h=3.5,col=plot.col['lines'],lwd=2)
       )
     )
legend('topright', legend=c('observed', 'theoretical'), col=plot.col[c('points','lines')], pch=c(1,NA), lwd=c(0,2))
setwd(dir.figures); export.plot(file.prefix=paste("sampling_mean_of_means_r",repetitions,sep=''), export.formats=export.formats.plots, width=7,height=5)

plot(1:max.n,
     sampling.stats$var.of.means,
     xlab='n (sample size)',
     ylab='variance of the means',
     main='variance of the sample means',
     col=plot.col['points'],
     type='p',
     panel.first=grid(col=plot.col['grid'],lty='solid'),
     pch=1
     )
lines(1:max.n, sampling.stats$var.of.means.theor,col=plot.col['lines'],lwd=2)
legend('topright', legend=c('observed', 'theoretical'), col=plot.col[c('points','lines')], pch=c(1,NA), lwd=c(0,2))
setwd(dir.figures); export.plot(file.prefix=paste("sampling_var_of_means_r",repetitions,sep=''), export.formats=export.formats.plots, width=12,height=8)

plot(1:max.n,
     sampling.stats$sd.of.means,
     xlab='n (sample size)',
     ylab='standard deviation of the means',
     main='standard deviation of the sample means',
     col=plot.col['points'],
     type='p',
     lwd=1,
     panel.first=grid(col=plot.col['grid'],lty='solid'),
     pch=1
     )
lines(1:max.n, sampling.stats$sd.of.means.theor,col=plot.col['lines'],lwd=2)
legend('topright', legend=c('observed', 'theoretical'), col=plot.col[c('points','lines')], pch=c(1,NA), lwd=c(0,2))
setwd(dir.figures); export.plot(file.prefix=paste("sampling_sd_of_means_r",repetitions,sep=''), export.formats=export.formats.plots, width=12,height=8)



## ##############################################################
## Fitting of normal curves on sample mean distributions taken from
## normal samples

X11(width=7, height=5)
plot.col['histo'] <- '#BBBBBB'
plot.col['lines'] <- '#000000'

## Generate a set of random numbers and store it is a matrix where
## each colomn represents one element of a given sample, and each row
## represents a repetition (from 1 to r).
n.values <- c(2, 3, 5, 10, 50)
max.n <- max(n.values) ## Number of columns in the table
r <- 10000 ## Number of repetitions
pop.mean <- 0
pop.sd <- 1
pop.var <- pop.sd^2
x <- matrix(rnorm(n=r*max.n,mean=pop.mean, sd=pop.sd), nrow=r, ncol=max.n)

## Draw an histogram of the total sample (size=max.n) over all repetitions
## and store the distribution in a variable called x.distrib 
x.distrib <- hist(x, breaks=100,plot=F)

## Plot a frequency polygon
plot(x.distrib$mids, x.distrib$density,
     main=paste("Distribution of sampled values (random normal, m=",pop.mean,"; sd=",pop.sd,sep=""),
     panel.first=grid(col=plot.col['grid'], lty='solid'),
     xlab="x",
     xlim=c(pop.mean-5*pop.sd,pop.mean+5*pop.sd),
     ylab="density",
     type='h',
     col=plot.col['histo'],
     lwd=4)

## Plot the theoretical distribution
lines(x.distrib$mids, dnorm(x.distrib$mids,mean=pop.mean, sd=pop.sd),type="l", col=plot.col['lines'],lwd=2)
legend('topright', legend=c('observed', 'theoretical'), col=plot.col[c('histo','lines')], lwd=c(4,2),bg='white',bty='o')

setwd(dir.figures); export.plot(file.prefix=paste("sampling_values_norm_fittinr", r,sep=''), export.formats=export.formats.plots, width=7,height=5)

## Calculate the distributions of sample means and fit the theoretical
## distribution
for (n in n.values) {

  ## Calculate the mean for each row of the sub-matrix (for each one of the r samples)
  ## You can call help(apply) to understand how the function apply() works
  sample.means <- apply(x[,1:n], 1, mean) ## The second argument (1) indicates that we apply the operation (mean) on each row
  
  ## Calculate the distribution of sample means
  sample.mean.distrib <- hist(sample.means, breaks=100, plot=F)
  ##  print(sample.mean.distrib)
  
  ## Plot the observed sample mean distribution
  plot(sample.mean.distrib$mids, sample.mean.distrib$density,
     main=paste('Distribution of sample mean, n=',n,sep=''),
     panel.first=grid(col=plot.col['grid'], lty='solid'),
     xlab="x",
     xlim=c(pop.mean-5*pop.sd,pop.mean+5*pop.sd),
     ylab="density",
     type='h',
     col=plot.col['histo'],
     lwd=4)

  ## Plot the theoretical distribution on top of the observed distribution
  lines(sample.mean.distrib$mids, dnorm(sample.mean.distrib$mids, mean=pop.mean, sd=pop.sd/sqrt(n)),type="l", col=plot.col['lines'],lwd=2)

  legend('topright', legend=c('observed', 'theoretical'), col=plot.col[c('histo','lines')], lwd=c(4,2),bg='white',bty='o')

  setwd(dir.figures); export.plot(file.prefix=paste('sampling_mean_norm_fitting_n',n,'_r',r,sep=''), export.formats=export.formats.plots, width=7,height=5)
}


## Calculate the distributions of variances and fit the theoretical distribution
n.values <- c(2,3,4,5,10,15,20,30,50)
for (n in n.values) {

  ## Calculate the var for each row of the sub-matrix (for each one of the r samples)
  ## You can call help(apply) to understand how the function apply() works
  sample.vars <- apply(x[,1:n], 1, var) ## The second argument (1) indicates that we apply the operation (var) on each row
  
  ## Calculate the distribution of sample vars
  sample.var.distrib <- hist(sample.vars, breaks=100, plot=F)


  ## Plot the observed sample var distribution
  plot(sample.var.distrib$mids, sample.var.distrib$density,
     main=paste('Distribution of sample var, n=',n,sep=''),
     panel.first=grid(col=plot.col['grid'], lty='solid'),
     xlab="x",
     xlim=c(0,max(sample.vars)),
     ylab="density",
     type='h',
     col=plot.col['histo'],
     lwd=4)

  ## Compute class intervals
  ci <- sample.var.distrib$breaks[2:length(sample.var.distrib$breaks)] -  sample.var.distrib$breaks[1:(length(sample.var.distrib$breaks)-1)]
  
  ## Compute the theoretical density
  var.density.theor <- (n/pop.var) * dchisq(sample.var.distrib$mids * n / pop.var, df=(n-1))
  sum(var.density.theor*ci) ## This should give ~ 1

  ## Plot the theoretical distribution on top of the observed distribution
  lines(sample.var.distrib$mids, var.density.theor,type="l", col=plot.col['lines'],lwd=2)

  legend('topright', legend=c('observed', 'theoretical'), col=plot.col[c('histo','lines')], lwd=c(4,2),bg='white',bty='o')

  setwd(dir.figures); export.plot(file.prefix=paste('sampling_var_chisq_fitting_n',n,'_r',r,sep=''), export.formats=export.formats.plots, width=7,height=5)
}
