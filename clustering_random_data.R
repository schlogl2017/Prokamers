################################################################
##
## Apply clustering algorithm on random data to show that the
## agglomeration rule has an impact, even when the data is not
## structured.
##
################################################################

## Generate a data frame with random numbers
x <- matrix(nrow=50, ncol=20,rnorm(n=50*20, m=0, sd=1))
x <- data.frame(x)

## calculate a distance matrix
d <- dist(x)

X11(width=10,height=8)
par(mfrow=c(2,2))

## Apply hierarchical clustering
plot(t <- hclust(d, method='single'), main="Random data - Single lignage",xlab="", cex=0.7)

## Apply hierarchical clustering
plot(t <- hclust(d, method='average'), main="Random data - Average lignage",xlab="", cex=0.7)

## Apply hierarchical clustering
plot(t <- hclust(d, method='complete'), main="Random data - Complete lignage",xlab="", cex=0.7)

## Apply hierarchical clustering
plot(t <- hclust(d, method='ward'), main="Random data - Ward lignage",xlab="", cex=0.7)

par(mfrow=c(1,1))

setwd(dir.figures); export.plot(file.prefix='hclust_linkage_effect_random_data', export.formats=export.formats.plots, width=10,height=8)

