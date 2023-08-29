## ##############################################################
## a VERY inefficient algorithm for clustering with mobile centers,
## mainly used for illustrative purposes (show the result at each
## step). The algorithm allows to specify a distance metrics.

mobile.centers <- function(x, ## a data frame
			   k=NA, ## number of mobile centers
			   centers=NA, ## initial centers
			   iter.max=10, ## Maximum number of iterations
			   dist.method="euclidian", ## Distance metrics
			   plot=T, ## Plot the result (in the two first dimensions of the data space)
			   ... ## Additional parameters are passed to the plot function
			   ) {
    
  result <- list()

  n <- dim(x)[1]
  p <- dim(x)[2]

  ## seed centers
  if (is.na(centers)) {
    if (is.na(k)) {
      stop ("either k or centers should be specified")
    } else {
      centers <-  x[sample(n,k),]
    }
  } else {
    k <- dim(centers)[1]
  }

  distances <- as.data.frame(matrix(nrow=dim(x)[1],ncol=k))
  cluster <- vector(length=dim(x)[1],mode="numeric")
  previous.cluster <- cluster

  for (iter in 1:iter.max) {
    
    ## calculate distances
    for (i in 1:n) {
      for (c in 1:k) {
        distances[i,c] <- dist(rbind(x[i,], centers[c,]),method=dist.method)
      }
    }
    
    ## assign class
    cluster <- apply (distances, 1, which.min)

    ## stop iterations if clusters do not evolve anymore
    if (sum(cluster == previous.cluster) == n) {
      break
    }
    previous.cluster <- cluster
    
    ## calculate new centers
    for (c in 1:k) {
      centers[c,] <- apply(x[cluster==c,],2,mean)
    }
    

  }

  ## Plot the result, in the first dimensions of the data space
  plot(x[,c(1,2)], 
       col=cluster+1, 
       main=paste('iter.max =', iter.max, '; iterations = ',iter),
       pch=19,
       xlab=NA,
       ylab=NA,
       ...
       )
  points(centers, col = '#0000DD', pch = 8)
  points(seed.centers, col='#000000',pch=19)

  ## calculate the result 
  result$iterations <- iter
  result$k <- k
  result$n <- n
  result$p <- p
  result$iter.max <- iter.max
  result$cluster <- cluster
  result$centers <- centers
  result$size <- as.vector(table(cluster))
  return(result)
}
