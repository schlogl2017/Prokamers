################################################################
## An underlying assumption of the chi2 is that the expected number of
## occurrences is >= 5 in each class.
##
## This function takes as input a vector of observed frequencies (must
## be absolute frequencies, i.e. occurrences), and a vector of
## expected frequencies, and merges the left and right tails,
## respectively, in order to ensure that the terminal classes have a
## frequency >= 5.
##
## Loading:
## source(file.path(dir.util, 'util_chi2_merged_tails.R'))


merge.tails <- function(obs, ## A vector of integer numbers (observed occurrences)
                        exp, ## A vector of real numbers (expected occurrences)
                        min.exp = 5 ## Min number of observations in each tail class
                        ) {
  result <- list()
  result$obs.ori <- obs
  result$exp.ori <- exp

  ## Check input data
  if (sum(exp) < min.exp) {
    stop('the sum of all expected values must be at least 5')
  }
  
  ## Merge the left tail
  exp.leftsum <- cumsum(exp)
  obs.leftsum <- cumsum(obs)
  first.left <-   min(which(exp.leftsum > min.exp))
  result$first.left <- first.left

  exp[first.left] <- exp.leftsum[first.left]
  exp <- exp[first.left:length(exp)]

  obs[first.left] <- obs.leftsum[first.left]
  obs <- obs[first.left:length(obs)]

  ## Merge the right tail
  exp.rightsum <- rev(cumsum(rev(exp)))
  obs.rightsum <- rev(cumsum(rev(obs)))
  last.right <-   max(which(exp.rightsum > min.exp))
  result$last.right <- last.right

  exp[last.right] <- exp.rightsum[last.right]
  exp <- exp[1:last.right]

  obs[last.right] <- obs.rightsum[last.right]
  obs <- obs[1:last.right]

  ## Adjust expexcted vector to ensure that the sum is same as the
  ## observed values vector (avoir rounding problems)
  exp <- exp * sum(obs) / sum(exp)
  
  ## return the result
  result$obs.merged <- obs
  result$exp.merged <- exp
  return(result)
}

## ##############################################################
## This implementation of the chi2 is implemented for didactic
## purposes.
##
## It is somewhat redundant with the chisq.test() funcion in stats
## library, but presents the following interest:
##
## 1) With table=T, returns a table with all the computation details
## (this can be useful to detect the classes that most influence the
## chi2.obs statistics)
##
## 2) The P-value is computed with the pchisq() function, which is not
## limited to 1e-16, in contrast with the function chisq.test()
##
## 3) The option min.exp automatically merges the left and/or right
## tail, respectively, to meet the assumption that all the exp should
## be >=5.
##
## Author: Jacques van Helden
chi2.test <- function(obs, ## A vector with observed occurences
                      exp, ## A vector with expected occurrences
                      min.exp=5, ## If not null, merge the tails
                      table=T
                      ) {

  ## Check input data
  if (abs(sum(obs) - sum(exp)) > 1) {
    stop(paste ("The sum of observed (",sum(obs), ") and expected (",sum(exp), ") cannot differ.", sep=""))
  }


  result <- list()
  if (!is.null(min.exp)) {
    result$obs.ori <- obs
    result$exp.ori <- exp
    merged <- merge.tails(obs,exp,min.exp=min.exp)
    obs <- merged$obs.merged
    exp <- merged$exp.merged
  }
  result$obs <- obs
  result$exp <- exp
  
  p <- exp/sum(exp)
  result$p <- p
  result$chisq.test <- chisq.test(obs,p=p)

  ## Manually compute the chi2
  diff  <- exp - obs
  diff2  <- diff^2
  chi2.vect <- diff2/exp
  chi2.obs <- sum(chi2.vect)
  
  
  ## add  the chi2 details  in a table  
  if (table) {
    result$table <- data.frame(
                               obs = obs,
                               exp = exp,
                               p = p,
                               diff  = diff,
                               diff.sq  = diff2,
                               chi2.vect = chi2.vect
                               )
  }

  result$chi2.obs <- chi2.obs
  result$chi2.df <- length(exp)-1
  result$chi2.Pval <- pchisq(q=chi2.obs, df=result$chi2.df,lower.tail=F)
  return(result)
}

## Quick test
# N <- 500
# p <- 0.02
# R <- 1000
# x <- rbinom(R,size=N,p=p)
# max.x <- max(x)
# b <- (0:(max.x+2))-0.5
# h <- hist(x,breaks=b,col='gray')
# obs <- h$counts
# exp <- R*dbinom(h$mids,size=N,prob=p)
# exp[length(exp)] <- sum(R*dbinom(max(h$mids):N,size=N,prob=p))

# chi2.test(obs,exp)


