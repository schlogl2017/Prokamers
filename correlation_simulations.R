################################################################
##
## Compare the properties of Pearson, spearman correlation + mutual
## information.
##
## TO DO/ compute mutual information and see if it assigns a higher score to the circular relationship

################################################################
##
## Study case #1: when the relationship is linear, Pearson and
## Spearman detect it equally well.
##
## Generate a data set with linear relationship: y = x + some.noise
n <- 10000
noise.sd <- 0.2

## Generate a data set that roughly mimics the distribution of affymetrix log-transformed data (range between 2 and 16, asymmetry)
x <- sqrt((abs(rnorm(mean=10,sd=2, n=n))))
x <- 16*(max(x) -x) / (max(x)  - min(x))
hist(x, breaks=100, main='Fake microarray data', xlab='Signal', ylab='Number of probesets')

some.noise <- rnorm(mean=0,sd=noise.sd, n=n)
y <- x + some.noise

## plot x,y
x11(width=7, height=7)
plot(x,y, main='linear relationship')

## Compute Pearson's correlation
print(cor.pearson <- cor(x,y, method="pearson"))
print(cor.spearman <- cor(x,y, method="spearman"))

################################################################
## Study case 2: when the datasets are related by a non-linear
## relationship, Spearman detects it well,

gamma <- 2
y <- x^gamma + some.noise

## plot x,y
x11(width=7, height=7)
plot(x,y, main=paste('exponential relationship, gamma=', gamma))

## Compute Pearson's correlation
print(cor.pearson <- cor(x,y, method="pearson"))
print(cor.spearman <- cor(x,y, method="spearman"))


################################################################
## generate a data set that is not monotonous
r <- (max(x) - min(x))/2
y <- sqrt(r^2 - (x-r)^2) + some.noise

## plot x,y
x11(width=7, height=7)
plot(x,y, main='circular relationship')

## Compute Pearson's correlation
print(cor.pearson <- cor(x,y, method="pearson"))
print(cor.spearman <- cor(x,y, method="spearman"))
