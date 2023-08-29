## ##############################################################
## Comparison of distance metrics with random data

n <- 133
p <- 11

x <- matrix(rnorm(n*p), nrow=n,ncol=p)


## ##############################################################
## Calculate dissimilarity matrix

## Euclidian distance
e <- dist(x, method="euclidian") 

## dissimilarity based on Pearson's correlation coefficient
c <- as.dist(1-cor(t(x)))

## average dot product
adp <- (as.dist(as.matrix(x) %*% t(as.matrix(x))))/p
adp <- max(adp) -adp ## Transform the similarity into a dissimilarity

## ##############################################################
## Compare the metrics

## Euclidian versus correlation
plot(e,c,
     main='metrics comparisons',
     xlab='Euclidian distance',
     ylab='1 -correlation'
     )
setwd(dir.figures); export.plot(file.prefix="rnorm_eucl_vs_corr", export.formats=export.formats.plots, width=8,height=8)

## Euclidian versus average dot product
plot(e,adp,
     main='metrics comparisons',
     xlab='Euclidian distance',
     ylab='max(avg dot product) -(avg dot product)'
     )
setwd(dir.figures); export.plot(file.prefix="rnorm_eucl_vs_dotp", export.formats=export.formats.plots, width=8,height=8)

## correlation versus average dot product
plot(c,adp,
     main='metrics comparisons',
     xlab='1 - (correlation)',
     ylab='max(avg dot product) -(avg dot product)'
     )
setwd(dir.figures); export.plot(file.prefix="rnorm_cor_vs_dotp", export.formats=export.formats.plots, width=8,height=8)

## ##############################################################
## Calculate a tree on the basis of the matrix

################################################################
## effect of the linkage method
par(cex.main=2)
par(cex=0.65)
par(mfrow=c(4,1))
for (m in c("single","average","complete","ward")) {
  tree <- hclust(e,method=m)
  plot(tree,main=paste("Random numbers (normal) ;",m,"linkage", "; Euclidian distance"), xlab="")
}
par(mfrow=c(1,1))
setwd(dir.figures); export.plot(file.prefix="rnorm_hierarchical_linkage_method", export.formats=export.formats.plots, width=16,height=12)

################################################################
## effect of the dissimilarity metric
m <- "complete"

for (m in c('single','average','complete', 'ward')) {
  par(mfrow=c(3,1))

  ## Euclidian distance
  tree <- hclust(e,method=m)
  plot(tree,main=paste("Random numbers (normal) ;",m,"linkage", "; Euclidian distance"), xlab="")
 
  ## Average dot product
  tree <- hclust(adp,method=m)
  plot(tree,main=paste("Random numbers (normal) ;",m,"linkage", "; Average dot product"), xlab="")

  ## Coefficient of correlation
  tree <- hclust(c,method=m)
  plot(tree,main=paste("Random numbers (normal) ;",m,"linkage", "; Correlation"), xlab="")
 
  par(mfrow=c(1,1)) 
  setwd(dir.figures); export.plot(file.prefix=paste("rnorm", m ,
						    "hierarchical_distance_metrics",
						    sep='_'),
				  export.formats=export.formats.plots, width=16,height=12)
}

