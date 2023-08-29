## ##############################################################
## generate a random data set


## ##############################################################
## Generate a multi-group data set with binomial distributions
##
random.binom <- function(
                         group.size=c(60,40), # group sizes
                         vars=5, # number of variables
                         proba = c(0.0002,p2=0.00005), # probabilties of success for the different groups
                         trials=10000 # number of trials for the binomial distribution
                         ) {

  g <- length(group.size)

  if (g < 2) {
    stop ("There must be at least two groups")
  }
  if (length(proba) != g) {
    stop ("The probability vector must have the same length as group.sizes")
  }
  
  x <- matrix(rbinom(vars*group.size[1],trials,proba[1]),nrow=group.size[1],ncol=vars)
  group <- rep(1,group.size[1])
  
  for (i in 2:g) {
    x <- rbind(x, matrix(rbinom(vars*group.size[i],trials,proba[i]),nrow=group.size[i],ncol=vars))
    group <- c(group,rep(i,group.size[i]))
  }

  x <- data.frame(x)
  x$group <- group
  return(x)
}


