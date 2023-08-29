################################################################
#
# Exercises on hypothesis testing
#

library(ctest)


################################################################
# Test of conformity of the mean
# Taking successively each column of a DNA chip result, test
# if the means of log-ratios equal 0.  

source("load_phosphate.R")
data <- na.omit(pho[,-1])
names(data) <- paste("chip",1:8,sep="")
data.title <- "phosphate expression data"

n<- dim(data)[1]
p<- dim(data)[2]

################################################################
# Execute the Student test manually.
#
sample<- data[,1]

####calculated the observed t-value
m<- mean(sample)
sd <- sd(sample)
n<- length(sample)
mu <- 0
t.obs <- abs(m-mu)/(sd/sqrt(n))

#### calculate the theoretical t-value
alpha <- 0.05
t.theor <- abs(qt(alpha/2,df=n-1))
accept <- t.obs < t.theor 

################################################################
# Apply  Student test to each column successively
#
result <- data.frame(mean=apply(data,2,mean),sd=apply(data,2,sd))
result$t.obs <- rep(NA,p)
result$df <- rep(NA,p)
result$proba <- rep(NA,p)
result$accept <- rep(NA,p)

for (i in 1:p) {
  stud <- t.test(data[,i],mu=0)
  result[i,"t.obs"] <- stud$statistic
  result[i,"accept"] <- abs(t.theor) > abs(stud$statistic)
  result[i,"df"] <- stud$parameter
  result[i,"proba"] <- stud$p.value
}
