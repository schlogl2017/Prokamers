################################################################
#
# Poisson-based distances
#
# Reference:
# van Helden (2004). Metrics for comparing regulatory sequences
# on the basis of pattern counts.  Bioinformatics. 2004 Feb
# 12;20(3):399-406


################################################################
## Calculate a similarity matrix 
## on the basis of Poisson probabilities
##
## For the time being, I assume independence of the p
## Poisson-distributed variables.
##
sim.poisson <- function (x, # a multivariate-poisson distributed dataset
                         product=T, ## return the geometric mean instead of arithmetic mean
#                         diag=T, # return the diagonal
                         upper=T # only return the upper triangle of the matrix
                         ) {
  x <- as.matrix(x)      # coerce the data into a matrix
  n <- dim(x)[1]         # number of objects
  p <- dim(x)[2]         # number of variables
  m <- apply(x, 2, mean) # vector of variable means
  d <- matrix(data=0,nrow=n,ncol=n) # similarity matrix
  for (i in 1:n) {
    for (j in i:n) {
      C <- pmin(x[i,],x[j,])    # common sites 
      p.C <- (1-ppois(C-1,m))^2	# probability to observe >= C common sites
      if (product) {
        sim <- 1-exp(sum(log(p.C))/p)
#        sim <- sum(log(p.C))/p
      } else {
        sim <- sum(1-p.C)/p
      }
      d[i,j] <- sim
      d[j,i] <- sim
    }
  }

  d <- as.data.frame(d)
  row.names(d) <- row.names(x)
  names(d) <- row.names(x)
#  d <- as.dist(d,diag=diag,upper=upper) ### has to b suppressed because it sets all the diagonal elements to 0
  attr(d,"method") <- "sim.poisson.product"
  if (product) {
    attr(d,"model") <- "product"
  } else {
    attr(d,"model") <- "additive"
  }
  attr(d,"Labels") <- row.names(x)
  attr(d,"Upper") <- upper
#  attr(d,"Diag") <- diag
  return(d)
}


################################################################
## Calculate a distance matrix 
## on the basis of Poisson probabilities
##
## For the time being, I assume independence of the p
## Poisson-distributed variables.
##
dist.poisson <- function (x, # a multivariate-poisson distributed dataset
#			  diag=T, # return the diagonal
			  upper=T,  # only return the upper triangle of the matrix
			  over=T, # calculate the distance in terms of over-representation [P >=x = 1 - F(x-1)]
                          product=F # take prod(1-P) instead of the sum(P). Beware, this returns a similarity !
			  ) {
  x <- as.matrix(x)      # coerce the data into a matrix
  n <- dim(x)[1]         # number of objects
  p <- dim(x)[2]         # number of variables
  m <- apply(x, 2, mean) # vector of variable means
  d <- matrix(data=0,nrow=n,ncol=n) # distance matrix
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (over) {
        ## calculate the distance in terms of over-representation [P >=x = 1 - F(x-1)]
        p.i <- ppois(x[i,]-1,m) # 1 - probability to observe >= x[i,]
        p.j <- ppois(x[j,]-1,m) # 1 - probability to observe >= x[j,]
      } else {
        ## calculate the distance as the difference between the left-tails
        p.i <- ppois(x[i,],m) # probability to observe <= x[i,]
        p.j <- ppois(x[j,],m) # probability to observe <= x[j,]
      }
      if (product) {
        d.ij <- prod(1-abs(p.i-p.j))^(1/p)
      } else {
        d.ij <- sum(abs(p.i-p.j))/p
      }
      d[i,j] <- d.ij
      d[j,i] <- d.ij
    }
  }
  d <- as.data.frame(d)
  row.names(d) <- row.names(x)
  names(d) <- row.names(x)
#  d <- as.dist(d,diag=diag,upper=upper) ### has to be suppressed because it sets diagonal to 0 
  attr(d,"method") <- "dist.poisson"
  if (over) {attr(d,"method") <- paste (attr(d,"method"), "over", sep=".")}
  if (product) {
    attr(d,"model") <- "product"
  } else {
    attr(d,"model") <- "additive"
  }
  attr(d,"Labels") <- row.names(x)
  attr(d,"Upper") <- upper
#  attr(d,"Diag") <- diag
  return(d)
}

################################################################
##
## Mixed metric calculated as the weighted difference between
## poisson-based similarity and distance.
##
sim.poisson.mixed <- function (x, ## a multi-variate table, one row per object, one column per variable
                               beta = 0, ## constant added to the difference
                               alpha = 1, ## weighting of the similarity substraction
                               over=T, # transmitted to dist.poisson
                               product=T, # transmitted to sim.poisson and dist.poisson
#			       diag=F, # only return the diagonal of the matrix
                               ...
                               ) {
  sim.p <- sim.poisson(x,product=product,...)
  dist.p <- dist.poisson(x,over=over,product=product,...)
#  if (product) {
#    s <- as.data.frame(sim.p * (1-dist.p))
#    s <- max(s) - s ### convert similarity into distance
#  } else {
  s <- as.data.frame(sim.p - alpha * dist.p + beta)
#  }

  ################################################################
  ## document the resulting object

  ## method
  attr(s,"method") <- "sim.poisson.mixed"

  ## product
  if (product) {
      attr(s,"product") <- T
      attr(s,"method") <- paste (attr(s,"method"), "product", sep=".")
  }

#  ## logarithms
#  if (log) {
#      attr(s,"log") <- T
#      attr(s,"method") <- paste (attr(s,"method"), "log", sep=".")
#  }

  ## over-representation
  if (over) {
      attr(s,"over") <- T
      attr(s,"method") <- paste (attr(s,"method"), "over", sep=".")
  }

  attr(s,"Labels") <- row.names(x)
#  attr(s,"Diag") <- diag
  return(s)
}


################################################################
##
## Conceptual illustration of the Poisson-based distance
## (This script draws the Figure 1 of the paper)
plot.poisson.example <- function (x = 3, 
                                  y = 5,
				  m = 3.5,
				  xmax=10,
				  over=T) {
  plot(0:(xmax+1)-0.5,dpois(0:(xmax+1),m),type="s",lwd=2,col=1,xlab='occ',ylab='proba')
  if (over) {
    lines(x:(y-1),dpois(x:(y-1),m),type="h",lwd=5,col="#008800",lty=1)
    lines(y:xmax,dpois(y:xmax,m),type="h",lwd=3,col=6,lty=2)
    dist <- abs(ppois(y-1,m)-ppois(x-1,m))
    legend(6,0.2,c("Poisson","A or more","B or more", paste("d=", dist, sep="")),col=c(1,"#008800",6,0),lwd=c(2,5,3,0),lty=c(1,1,2),cex=0.8)
  } else {
    lines(0:x,dpois(0:x,m),type="h",lwd=3,col=4,lty=2)
    lines((x+1):y,dpois((x+1):y,m),type="h",lwd=5,col=2,lty=1)
    dist <- abs(ppois(x,m)-ppois(y,m))
    legend(6,0.2,c("Poisson","Common occurrences","Distinct occurrences",paste("d=",dist,sep="")),col=c(1,4,2,0),lwd=c(2,3,5,0),lty=c(1,2,1),cex=0.8)
  }

  return(dist)
}

