################################################################
##
## Illustrate the principles of estimation by infering yeast orf 
## lengths from the 3rd chromosome genes only
##
##
## Author: Jacques van Helden
##
## Running this script requires to first run the script config.R
##
## Single-command execution:
## source(file.path(dir.R.files, 'central_tendency_yeast_orf_lengths.R'))

#### read the 3rd chromosome data
sample.file <- file.path(dir.data, "orf_lengths", "yeast_chr3_orf_length_enum.txt")
sample  <- scan(sample.file)
sample.size  <- length(sample)  #size of the sample
sample.mean  <- mean(sample)    # mean of the sample
sample.var <- sum(sample^2)/sample.size - sample.mean^2  # variance of the sample
sample.sd  <- sqrt(sample.var)          # standard deviation of the sample
sample.stand.err <- sample.sd/sqrt(sample.size) # standard error on the mean

## Estimate population parameters from the sample
est.pop.mean <- sample.mean
est.pop.var <- sample.var * sample.size/(sample.size-1)
est.pop.sd <- sqrt(est.pop.var)
est.stand.err <- sample.sd / sqrt(sample.size - 1)

## Confidence interval around the mean
conf.alpha <- 0.05
conf.t.theor <- qt(p=1-conf.alpha/2,df=(sample.size-1)) 
conf.d <- conf.t.theor * sample.stand.err
conf.min <- sample.mean - conf.d
conf.max <- sample.mean + conf.d
  
#### read the whole genome data
pop.file <- file.path(dir.data, "orf_lengths", "yeast_orf_length_enum.txt")
pop  <- scan(pop.file)
pop.size  <- length(pop)  #size of the sample
pop.mean  <- mean(pop)    # mean of the sample
pop.var <- sum(pop^2)/pop.size - pop.mean^2  # variance of the sample
pop.sd  <- sqrt(pop.var)          # standard deviation of the sample

#### print the resuls
result <- data.frame(
                     sample.size = sample.size,
                     sample.mean=sample.mean,
                     sample.median = median(sample),
                     sample.var = sample.var,
                     sample.sd = sample.sd,
                     sample.stand.err=sample.stand.err,

                     est.pop.mean = est.pop.mean,
                     est.pop.sd = est.pop.sd,
                     est.pop.var = est.pop.var,

                     conf.alpha = conf.alpha,
                     conf.d = conf.d,
                     conf.min = conf.min,
                     conf.max = conf.max,
                     
                     pop.size = pop.size,
                     pop.mean=pop.mean,
                     pop.median = median(pop),
                     pop.sd = pop.sd,
                     pop.var=pop.var
                     )

print(result)
