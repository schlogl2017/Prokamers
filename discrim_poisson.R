## ##############################################################
## Bayesian supervised classification based on an assumption of
## Poisson-distributed variables.
##
## source(paste(dir.util,'discrim_poisson.R',sep='/'))

library(mda) # for confusion

discrim.poisson <- function(x,group.labels, CV=F, ...) {
  if (CV) {
    result <- discrim.poisson.loo(x,group.labels, ...)
  } else {
    result <- discrim.poisson.test(x,group.labels,x,group.labels)
  }
  return(result)
}


## ##############################################################
## Leave-one-out with discrim.poisson
discrim.poisson.loo <- function(x,group.labels,...) {
  result <- list()
  groups <- names(table(group.labels))
  class <- rep(NA, length(group.labels))

  ## vector for storing posterior probabilities
  post <- matrix(ncol=length(groups), nrow=nrow(x))
  post <- as.data.frame(post)
  row.names(post) <- row.names(x)

  ## vector for storing inverse probabilities
  inv <- matrix(ncol=length(groups), nrow=nrow(x))
  inv <- as.data.frame(inv)
  row.names(inv) <- row.names(x)


  message.verbose("Cross-validation by Leave-one-out ",2)
  for (i in 1:length(group.labels)) {
    if (!is.na(group.labels[i])) {
      message.verbose(paste("Leave-one-out ", i),3)
      disc <- discrim.poisson.test(x[-i,],group.labels[-i], test=x[i,],test.labels=group.labels[i])
      class[i] <- disc$class
      inv[i,] <- disc$P.inverse
      post[i,] <- disc$P.posterior
    }
  }
  names(post) <- names(disc$P.posterior)
  names(inv) <- names(disc$P.inverse)

  result$lev <- groups
  result$P.inverse <- inv
  result$P.posterior <- post
  result$class <- class
  result$conf <- confusion(class, group.labels)
  return(result)
}


## ##############################################################
## Discriminant analysis for ordinal data, based on a Poisson
## distribution
##
## TO DO: use sum of logs to increase precision
discrim.poisson.test <- function (x, # training set. This should be a multivariate Poisson distributed dataset
                                  group.labels, # a vector containing one label for each element of the training set
                                  test=NULL, # (optional) testing set. This should be a multivariate Poisson distributed dataset
                                  test.labels = NULL,  # (optional) a vector containing one label for each element of the testing set
                                  prior = NA # prior group membership probabilities
                                  ) {

  train <- as.data.frame(x[!is.na(group.labels),])      # coerce the data into a matrix
  train.labels <- group.labels[!is.na(group.labels)]

  message.verbose("Poisson-based Bayesian classification",3)

  ## data dimensions
  train.n <- dim(train)[1]         # number of training objects
  p <- dim(train)[2]         # number of variables
  groups <- levels(as.factor(train.labels))

  ## counts and priors
  counts <- table(as.matrix(train.labels))
  if (is.na(prior)) {
    prior <- counts/sum(counts)
  } else {
    ## check dimension
    if (length(prior) != length(groups)) {
      stop ("The vector of priors must have the same length as the number of groups")
    }
    ## check prior names
    if (sort(names(prior)) != sort(groups)) {
      stop ("prior must be a named vector, with one entry per group")
    }
  }

  ## ##############################################################
  ## calculate group means
  means <- matrix(nrow=length(groups),ncol=p)
  row.names(means) <- groups
  for (g in groups) {
    means[g,] <- apply(train[train.labels==g,],2,mean)
  }
  means <- as.data.frame(means)
  names(means) <- names(train)


  ## ##############################################################
  ## initialize the result
  result <- list()
  result$prior <- prior
  result$counts <- counts
  result$means <- means
  result$lev <- groups
  result$N.train <- train.n

  ## ##############################################################
  ## Analyse the testing set
  if (!is.null(test)) {
    message.verbose("Classifying test data",3)
    test.n <- dim(test)[1]
    result$N.test <- test.n

    ## check number of variables
    test.p <- ncol(test)
    if (test.p != p) {
      stop ("Testing set must have the same number of variabes as training set")
    }

    ## ##############################################################
    ## calculate inverse probability P(x|g)
    ## initialize probability tables
    message.verbose("Calculating inverse probability",4)
    P.inverse <- matrix(nrow=test.n,ncol=length(groups))
    P.inverse <- as.data.frame(P.inverse)
    row.names(P.inverse) = row.names(test)
    names(P.inverse) = groups
    col.probas <- matrix(nrow=test.n,ncol=p)
    for (g in groups) {
      for (j in 1:p) {
        m <- means[g,j]
        col.probas[,j] <- dpois(as.matrix(test)[,j],m)
      }
      P.inverse[,g] <- apply(col.probas,1,prod)
    }

    ## ##############################################################
    ## calculate posterior probability P(g|x)
    ## P(g|x) = P(x|g)*P(g)/P(x)
    ## initialize probability tables
    message.verbose("Calculating posterior probability",4)
    P.posterior <- P.inverse
    P.x <- crossprod(t(P.inverse),as.matrix(prior))
    names(P.posterior) <- groups
    for (g in groups) {
      P.posterior[,g] <- P.posterior[,g]*prior[g]/P.x
    }

    result$P.inverse <- P.inverse
    result$P.x <- P.x
    result$P.posterior <- P.posterior

    ## ##############################################################
    ## predict class membership
    message.verbose("Predicting class",4)

    ### WARNING: I still have to improve the treatment of the cases where all the P.inverse=0
    if(sum(P.inverse)==0) {
      result$class <- NA
    } else {
      result$class <- groups[apply (P.posterior, 1, which.max)]
    }

    ## ##############################################################
    ## calculate confusion table
    if ((!is.null(test.labels)) && (length(test.labels) > 1)) {
      message.verbose("Confusion table",4)
      result$conf <- confusion (result$class, test.labels)
    }
  }

  message.verbose("Poisson-based Bayesian classification done",3)
  ## return the result
  return(result)
}


#  ## ##############################################################
#  ## Leave-one-out evaluation of Poisson-based bayesian discriminant
#  ## analysis
#  discrim.poisson.cv.old <- function (x, # a multivariate-poisson distributed dataset
#                                      train.labels # a vector containing one label for each element of x
#                                      ) {
#    x <- as.data.frame(x)      # coerce the data into a matrix
#    n <- dim(x)[1]         # number of objects
#    p <- dim(x)[2]         # number of variables
#    groups <- names(table(train.labels))
#
#    counts <- table(train.labels)
#    prior <- counts/sum(counts)
#    means <- matrix(nrow=length(groups),ncol=p)
#    row.names(means) <- groups
#    for (g in groups) {
#      means[g,] <- apply(x[train.labels==g,],2,mean)
#    }
#    means <- as.data.frame(means)
#    names(means) <- names(x)
#
#
#    P.inverse <- matrix(nrow=0,ncol=length(groups))
#    P.posterior <- matrix(nrow=0,ncol=length(groups))
#
#    for (i in 1:n) {
#      evaluation <- discrim.poisson.test(train=x[-i,],
#                                         train.labels = train.labels[-i],
#                                         test=x[i,])
#      P.inverse <- rbind(P.inverse,evaluation$P.inverse)
#      P.posterior <- rbind(P.posterior,evaluation$P.posterior)
#    }
#
#    P.inverse <- as.data.frame(P.inverse)
#    row.names(P.inverse) = row.names(x)
#    names(P.inverse) = groups
#
#    P.posterior <- as.data.frame(P.posterior)
#    row.names(P.posterior) = row.names(x)
#    names(P.posterior) = groups
#
#    class <- groups[apply (P.posterior, 1, which.max)]
#    conf <- confusion (class, train.labels)
#
#    result <- list()
#    result$prior <- prior
#    result$counts <- counts
#    result$means <- means
#    result$lev <- groups
#    result$N <- n
#    result$P.inverse <- P.inverse
#    result$P.posterior <- P.posterior
#    result$class <- class
#    result$confusion <- conf
#    return(result)
#  }
