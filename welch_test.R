## ##############################################################
## Apply the Welch version of the t-test to compare the mean of two samples
welch.test <- function (x1, ## Enumeration of values from the first sample
                        x2 ## Enumeration of values from the second sample
                        ) {

  ## Estimate parameters from the first sample
  mean1 <- mean(x1) ## mean of the sample
  n1 <- length(x1)
  ssd1 <-  sum((x1-mean1)^2) ## Sum of squared deviations for sample 1
  var.sample1 <-ssd1/n1 ## var of the sample 1 (with n)
  sd.sample1 <- sqrt(var.sample1)
  var.pop1 <- var(x1) ## estimator of the var of the population 1 (with n-1)
  sd.pop1 <- sd(x1) ## estimator of the sd of the population 1 (with n-1)
  st.err1 <- sd.pop1/sqrt(n1) ## Estimation of the dispersion of the mean

  ## Estimate parameters from the second sample
  mean2 <- mean(x2) ## mean of the sample
  n2 <- length(x2)
  ssd2 <-  sum((x2-mean2)^2) ## Sum of squared deviations for sample 1
  var.sample2 <-ssd2/n2 ## var of the sample 2 (with n)
  sd.sample2 <- sqrt(var.sample2)
  var.pop2 <- var(x2) ## estimator of the var of the population 1 (with n-1)
  sd.pop2 <- sd(x2) ## estimator of the sd of the population 2 (with n-1)
  st.err2 <- sd.pop2/sqrt(n2) ## Estimation of the dispersion of the mean

  ## Compare the populations
  diff <- mean1 - mean2
  diff.sd <- sqrt(var.sample1/(n1-1)+var.sample2/(n2-1))
  t.obs <- diff/diff.sd
  df <-(var.sample1/(n1-1) + var.sample2/(n2-1))^2/(((var.sample1/(n1-1))^2)/(n1-1)+((var.sample2/(n2-1))^2)/(n2-1))
  p.val <- pt(q=abs(t.obs),df=df,lower.tail=F)
  result <- list()

  ## Report the result
  result$parameters <- data.frame(n=c(n1,n2),
                                  mean=c(mean1,mean2),
                                  var.sample=c(var.sample1,var.sample2),
                                  sd.sample=c(sd.sample1,sd.sample2),
                                  st.err=c(st.err1,st.err2)
                                  )
  names(result$parameters) <- c("x1","x2")

#   result$n1 <- n1
#   result$mean1 <- mean1
#   result$var.sample1 <- var.sample1
#   result$sd.sample1 <- sd.sample1
#   result$st.err1 <- st.err1

#   result$n2 <- n2
#   result$mean2 <- mean2
#   result$var.sample2 <- var.sample2
#   result$sd.sample2 <- sd.sample2
#   result$st.err2 <- st.err2

  result$diff <- diff
  result$diff.sd <- diff.sd
  result$t.obs <- t.obs
  result$df <- df
  result$p.val <- p.val
  return(result)
}

