################################################################
## Definition of Poisson-based multivariate distance and 
## similarity metrics. 
## 
## Comparison between different multivariate similarity and
## distance metrics. 
##
################################################################

## source(paste(dir.util,'util_distances.R', sep='/'))

#### load multivariate analysis library
library(stats)
source(paste(dir.util,'poisson_distances.R', sep='/'))

#### define which metrics are distances/similarities
metrics.dist <- c("binary.distance", 
	          "canberra.distance", 
	          "poisson.distance.distinct", 
	          "poisson.distance.over",
	          "manhattan.distance", 
	          "euclidian.distance")

metrics.mixed <- c(
                  "poisson.mixed.over",
                  "poisson.mixed.distinct",
                  "poisson.mixed.over.product",
                  "poisson.mixed.distinct.product",
#                 "park.similarity.v1",
                 "park.similarity"
                  )

metrics.sim <- c(
                 "correlation.coefficient",
                 "poisson.similarity",
	         "poisson.similarity.product"
                 )


################################################################
##
## Similarity measurement as defined by Park et al 
## (second submission to bioinformatics, 2002)
## This is only calculated for one type of sites, 
## the combination of putative and known sites can be done 
## by calling twice this function. 
##
sim.park <- function(x, # a multivariate-poisson distributed dataset
#                     diag=T, # return the diagonal
                     upper=T, # only return the upper triangle of the matrix
                     alpha = 0.1,
                     beta = 0.5
                     ) {
  x <- as.matrix(x)      # coerce the data into a matrix
  n <- dim(x)[1]         # number of objects
  p <- dim(x)[2]         # number of variables
  m <- apply(x, 2, mean) # vector of variable means
  d <- matrix(data=NA,nrow=n,ncol=n) # similarity matrix
  for (i in 1:n) {
    for (j in i:n) {
      C <- pmin(x[i,],x[j,])            # common sites 
      T <- x[i,]+x[j,]			# total sites
      D <- (2/T + alpha)*C		# similarity
      P <- beta * T *(C == 0)		# penalty
      S <- sum(na.omit((D-P)/sqrt(m)))	# score
      d[j,i] <- S
      d[i,j] <- S
    }
  }
  d <- as.data.frame(d)
  row.names(d) <- row.names(x)
  names(d) <- row.names(x)
#  d <- as.dist(d,diag=diag,upper=upper) ### has to b suppressed because it sets all the diagonal elements to 0
  attr(d,"Labels") <- row.names(x)
  attr(d,"Upper") <- upper
#  attr(d,"Diag") <- diag
  return(d)
}


################################################################
##
## Similarity measurement as defined by Park et al 
## (second submission to bioinformatics, 2002)
## This is only calculated for one type of sites, 
## the combination of putative and known sites can be done 
## by calling twice this function. 
##
sim.park.first.version <- function(x, # a multivariate-poisson distributed dataset
#                                   diag=T, # return the diagonal
                                   upper=T, # only return the upper triangle of the matrix
                                   alpha = 0.1
                                   ) {
  x <- as.matrix(x)      # coerce the data into a matrix
  n <- dim(x)[1]         # number of objects
  p <- dim(x)[2]         # number of variables
  m <- apply(x, 2, mean) # vector of variable means
  d <- matrix(data=NA,nrow=n,ncol=n) # similarity matrix
  for (i in 1:n) {
    for (j in i:n) {
      g <- pmin(x[i,],x[j,])  # common regulatory elements
      C <- g / sqrt(m)	# weighted common regulatory elements
      T <- x[i,]+x[j,]	# total of regulatory elements
      S <- sum(na.omit(C*(2/T + alpha)))
      d[j,i] <- S
      d[i,j] <- S
    }
  }
  d <- as.data.frame(d)
  row.names(d) <- row.names(x)
  names(d) <- row.names(x)
#  d <- as.dist(d,diag=diag,upper=upper) ### has to b suppressed because it sets all the diagonal elements to 0
  attr(d,"Labels") <- row.names(x)
  attr(d,"Upper") <- upper
#  attr(d,"Diag") <- diag
  return(d)
}

################################################################
#
# Comparison of different distance and similarity metrics
# Returns a data frame with one column per metric, and one row 
## per pair of objects.
# Usage:
#    comp <- compare.distances(x)
# 
compare.distances <- function(x, # a data frame containing pattern counts 
			      alpha=1 # weighting of the Poisson distance in the mixed metrics
			      ) {
  #### similarity metrics
  poisson.similarity <-  as.vector(as.matrix(sim.poisson(x,product=F)))   ## Poisson similarity with sum of proba
  poisson.similarity.product <-  as.vector(as.matrix(sim.poisson(x,product=T)))   ## Poisson similarity with product of probabilities
  park.similarity.v1 <-  as.vector(as.matrix(sim.park.first.version(x)))   ## Park distance, first submission
  park.similarity <-  as.vector(as.matrix(sim.park(x)))   ## Park distance, second submission
  correlation.coefficient <- as.vector(cor(t(x)))   ## correlation 

  #### distance metrics
  poisson.distance.distinct <- as.vector(as.matrix(dist.poisson(x,upper=T,over=F)))   ## Poisson distance in terms of proba difference
  poisson.distance.over <- as.vector(as.matrix(dist.poisson(x,upper=T,over=T)))   ## Poisson distance in terms of over-representation
  binary.distance <-  as.vector(as.matrix(dist(x,method="binary",diag=T,upper=T))) ## binary distance
  canberra.distance <- as.vector(as.matrix(dist(x,method="canberra",diag=T,upper=T))) ## canberra distance
  manhattan.distance <- as.vector(as.matrix(dist(x,method="manhattan",diag=T,upper=T)))  ## manhattan distance
  euclidian.distance <- as.vector(as.matrix(dist(x,method="euclidian",diag=T,upper=T))) ## euclidian distance

  #### mixed metrics, combining distance and similarity
  poisson.mixed.distinct <-  poisson.distance.distinct - alpha * poisson.similarity
  poisson.mixed.over <-  poisson.distance.over - alpha * poisson.similarity
  poisson.mixed.distinct.product <-  poisson.distance.distinct - alpha * poisson.similarity.product
  poisson.mixed.over.product <-  poisson.distance.over - alpha * poisson.similarity.product

  ## compare distances
  compa <- data.frame(
                      poisson.similarity=poisson.similarity,
                      poisson.similarity.product=poisson.similarity.product,
                      park.similarity.v1 =park.similarity.v1 ,
                      park.similarity=park.similarity,
                      correlation.coefficient=correlation.coefficient,

                      poisson.distance.distinct=poisson.distance.distinct,
                      poisson.distance.over=poisson.distance.over,
                      binary.distance=binary.distance,
                      canberra.distance=canberra.distance,
                      manhattan.distance=manhattan.distance,
                      euclidian.distance=euclidian.distance,
                      poisson.mixed.distinct=poisson.mixed.distinct,
                      poisson.mixed.over=poisson.mixed.over,
                      poisson.mixed.distinct.product=poisson.mixed.distinct.product,
                      poisson.mixed.over.product,poisson.mixed.over.product
                      )
  
  return(compa)

}


## ##############################################################
##
## Hierarchical Clustering on the basis of different metrics
## and different clustering methods
##
pattern.count.clust <- function (x,
                                 prefix="clustering",
                                 dist.metric="poisson", ## poisson, manhattan, euclidian, ...
                                 clust.method="ward", ## see help(hclust) for a list of supported methods
                                 export=F,
                                 cex=0.7,
				 export.formats=c("ps"),
				 width=25,
				 height=10
                                 ) {
  par(cex=cex)
  title <- paste(prefix, dist.metric, clust.method, sep=" - ")
  dim(x)
  if (dist.metric=="poisson") {
    d <- sim.poisson.mixed(x)
  } else {
    d <- dist(x,method=dist.metric)
  }
  hc <- hclust(dist(d),method=clust.method)
  plot(hc,main=title)
  if (export) export.plot(file.prefix=paste("hclust", prefix, dist.metric, clust.method, sep="_"), export.formats=export.formats,width=width,height=height)
  par(cex=1)
  return (hc)
}



################################################################
#### calculate distances between pairs of elements
#### several metrics are supported
profile.pair.distances <- function(profiles, ## A data frame containing a multivariate table
                                   pairs, ## A two-column data frame, one row per pair, two elements per row
                                   metrics=c("euclidian","manhattan") ## Metrics for the calculation of distance
                                   ) {
  
  for (metric in metrics) {
    verbose (paste("Calculating pair distances, metric=", metric),1)
    
    ## initialise distance
    pairs[,metric] <- rep(-1,pair.nb)
    
    ## pairs
    for (i in 1:pair.nb) {
      e1 <- as.vector(pairs[i,1])
      e2 <- as.vector(pairs[i,2])
      verbose (paste(i,e1,e2),2)
      if (metric == "dist.poisson.distinct") {
        pairs[i,metric] <- dist.poisson(profiles[c(e1,e2),],over=F)[1,2]
      } else if (metric == "dist.poisson.over") {
        pairs[i,metric] <- dist.poisson(profiles[c(e1,e2),],over=T)[1,2]
      } else if (metric == "sim.poisson.add") {
        pairs[i,metric] <- 1 - sim.poisson(profiles[c(e1,e2),],product=F)[1,2]
      } else if (metric == "sim.poisson.product") {
        pairs[i,metric] <- 1 - sim.poisson(profiles[c(e1,e2),],product=T)[1,2]
      } else if (metric == "hypergeometric") {
        p1 <- profiles[e1,]
        p2 <- profiles[e2,]
        common.counts <- sum(p1 & p2)
        min.count <- min(counts[e1],counts[e2])
        max.count <- max(counts[e1],counts[e2])
        pairs[i,metric] <- phyper(q=common.counts-1,
                                       m=min.count,
                                       n=p-min.count,
                                       k=max.count,
                                       lower.tail=F,
                                       log.p=F)
      } else if (metric == "hypergeometric.log") {
        p1 <- profiles[e1,]
        p2 <- profiles[e2,]
        common.counts <- sum(p1 & p2)
        min.count <- min(counts[e1],counts[e2])
        max.count <- max(counts[e1],counts[e2])
        pairs[i,metric] <- phyper(q=common.counts-1,
                                       m=min.count,
                                       n=p-min.count,
                                       k=max.count,
                                       lower.tail=F,
                                       log.p=T)
      } else if (metric == "binary.rel") {
        p1 <- profiles[e1,]
        p2 <- profiles[e2,]
        max.count <- max(counts[e1],counts[e2])
        common.counts <- sum(p1 & p2)
        pairs[i,metric] <- 1 - common.counts/max.count
      } else if (metric == "common") {
        p1 <- profiles[e1,]
        p2 <- profiles[e2,]
        pairs[i,"common"] <- sum(p1 & p2)
      } else {
        pairs[i,metric] <- dist(profiles[c(e1,e2),], method=metric)
      }
    }
  }
  return (pairs)
}

