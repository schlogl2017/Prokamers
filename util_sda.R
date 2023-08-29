##############################################################
##
## Stepwise discriminant analysis procedures
##
################################################################
## source(file.path(dir.util, "util_sda.R"))

## load the multivariate analysis libraries
library(stats)
library(mda)
library(MASS) ## Contains lda, qda

## load utilities
source(file.path(dir.util, "util.R"))
source(file.path(dir.util, "util_plots.R"))

################################################################
##
## A common call for different methods of PDA
## (predictive discriminant analysis)
##
################################################################
pda <- function (data,                # a data frame
                 group.labels,        # a vector containing the training labels
                 test.labels = NULL,  # a vector containing labels for an independent testing set
                 pda.method="lda",        #  "lda" (linear) or "qda" (quadratic)
                 loo=T,               # leave-one-out evaluation
		 conf=T,               # calculate the confusion matrix
                 prior=NA
                 ) {

  message.verbose("predictive discriminant analysis", 3)
  message.verbose(paste("prior probabilities", prior), 3)
  message.verbose(paste("LOO", loo), 3)
  message.verbose(paste("confusion", conf), 3)


  if (length(test.labels==0)) {
    test.labels <- NULL;
  }

  data <- as.matrix(data)
  ##  data.training <- as.matrix(data)
  data.training <- as.matrix(data[!is.na(group.labels),])
  labels.training <- group.labels[!is.na(group.labels)]

  ## check that there is no variable with all values identical, since
  ## this is not supported by lda() and qda()
  sd.per.var <- apply(data.training,2,sd)
  if (min(sd.per.var) == 0) {
    message.verbose("pda() problem: some variables are constant for all the training objects. Skipped, an error value arbitrarily set to 1.2. ", 0)
    result <- list()
    result$confusion <- list()
    attr(result$confusion, "error") <- 1.2
    result$sd.per.var <- sd.per.var
    return (result)
  } else if (pda.method == "qda") {
    if (is.na(prior)) {
      result <- qda(data.training,labels.training,na.action=na.omit,CV=loo)
    } else {
      result <- qda(data.training,labels.training,na.action=na.omit,CV=loo,prior=prior)
    }
  } else if (pda.method == "lda") {
    if (is.na(prior)) {
      result <- lda(data.training,labels.training,na.action=na.omit,CV=loo)
    } else {
      result <- lda(data.training,labels.training,na.action=na.omit,CV=loo,prior=prior)
    }
  } else if (pda.method == "poisson") {
    if (is.na(prior)) {
      result <- discrim.poisson(data,group.labels,CV=loo)
    } else {
      result <- discrim.poisson(data,group.labels,CV=loo, prior=prior)
    }
  } else {
    print("ERROR: unknown PDA method")
    return(NULL)
  }

  ## add some information
  result$training.names <- names(labels.training)
  result$training.label <- as.factor(labels.training)
  result$pda.method <- pda.method
  result$loo <- loo

  ## calculate confusion table
  if (pda.method != "poisson") {
    if (conf) {

      ## Confusion table with the training set
      if (loo) {
        ## Leave-one-out validation
        result$validation <- "LOO"
        message.verbose("PDA - LOO validation", 4)
        result$confusion <- confusion(result$class,result$training.label)
      } else {
        ## Internal validation
        result$validation <- "internal"
        result$class <- predict(result,as.matrix(data.training))$class
        verbose ("PDA - internal validation",4)
        result$confusion <- confusion(result$class,result$training.label)
      }

      ##      } else {
      if (!is.null(test.labels)) {
        ## Confusion table with the trainng set
        result$validation <- "testing set"
        verbose ("PDA - calculating confusion with testing set", 0)
        if (loo) {
          verbose ("PDA - treating testing set - LOO", 0)
          test.result <- pda(data,group.labels,test.labels, loo=F, conf=T,prior=prior)
          result$confusion.test <- test.result$confusion
        } else {
          verbose ("PDA - treating testing set - no LOO", 0)
          result$test.pred <- predict(result, data)
          result$test.labels <- test.labels
          result$confusion.test <- confusion(result$test.pred$class,test.labels)
        }
      }
    }
  }
  return(result)
}

################################################################
## Compute the hit rate by comparing vectors of
## known and predicted classes.
compute.hit.stats <- function (known.class,
                               predicted.class) {
  result <- list()

  ## Build a confusion table of known versus predicted class
  result$confusion <- table(known.class, predicted.class)
  result$confusion.tab <- as.data.frame(result$confusion[,]) ## the [,] is necessary to keep the structure of the table
  result$known.sum <- apply(result$confusion.tab, 1, sum)
  result$pred.sum <- apply(result$confusion.tab, 2, sum)

  ## Add margins to the confusion table
  result$confusion.tab[, "known.sum"] <- result$known.sum
  result$confusion.tab["pred.sum",] <- apply(result$confusion.tab, 2, sum)

  ## Number of predictions
  nb.pred <- sum(!is.na(predicted.class)) ## Total number of predictions
  nb.nas <- sum(is.na(predicted.class)) ## NA values

  ## Hit rate
  hits <- known.class == predicted.class
  nb.hits <- sum(na.omit(hits))
  hit.rate <- round(nb.hits / nb.pred, digits=3)

  ## Error rate
  errors <- known.class != predicted.class
  nb.errors <- sum(na.omit(errors))
  error.rate <- round(nb.errors/ nb.pred, digits=3)

  ## Compute the number of hits
  ## (we need to omit NA values because LDA fails to assign a group to some objects).
  result$rates <- c("nb.pred" = nb.pred,
                    "nb.hits" = nb.hits,
                    "nb.errors" = nb.errors,
                    "nb.nas" = nb.nas,
                    "hit.rate" = hit.rate,
                    "error.rate" = error.rate)
  return(result)
}


################################################################
## Predictive discriminant analysis with stepwise forward variable
## selection by progressively adding variables (in a specified order)
## and evaluating the error rate at each step
stepwise.pda <- function (data,              # a data frame
                          group.labels,      # a character vector of length N, containing group labels
                          test.labels = NULL,  # a vector containing labels for an independent testing set
                          data.title=NULL,     # a prefix for the gaph titles and export file names
                          var.order = NULL,    # variable ordering
                          max.p = NA,	     # max number of variables to incorporate
			  df=3,              # degrees of freedom for qda.check.variables()
                          pda.method="lda",  # discriminant analysis method
                          stepwise.method="forward", # variable selection method
                          ## supported: "forward", "single var"
                          ##                          loo=T,	     # leave-one-out validation
                          display.plot=T,    # plot the error rate as a function of the number of variables
			  last.min=T,         # when several combination of variables have the same error rates, retain the last minimu m (with the highes number of variables)
                          pc.transform=F,
                          ...
                          ) {


  ## coerce the data to a data frame
  ## and select objects with non-NA labels
  result <- list()
  nona <- !is.na(group.labels)
  data <- data.frame(data[nona,])
  group.labels <- group.labels[nona]
  n <- nrow(data)
  p <- ncol(data)

  ## max number of variables
  if (is.na(max.p)) {
    max.p <- p
  }
  if (pda.method == "qda") {
    min.n <- min(table(group.labels)) # nb of elements in the smallest group
    max.p <- min(max.p,min.n-df)
                                        #    max.p <- min(max.p,p)
  }
  result$max.p <- max.p

  ## initialize the error table
  groups <- levels(as.factor(group.labels))
  g <- length(groups)
  result$errors <- vector()
  if (!is.null(test.labels)) {
    result$errors.test <- vector()
    result$confusion.table <- matrix(nrow=0,ncol=4+g^2+1)
  } else {
    result$confusion.table <- matrix(nrow=0,ncol=4+g^2)
  }

  ## single-variable based ordering
  if ((is.null(var.order)) &&
      (stepwise.method=="single var")){
    o <- variable.ordering(data,group.labels,single.var=T, ...)
    var.order <- o$order.single.var
  }

  ## verbose
  verbose (paste(date(), "Starting stepwise PDA",pda.method,max.p,sep="        "),1)

  ## calculate error rate for increasing nubers of variables
  selected.vars <- vector(mode="numeric")
  error.min <- 2 # initialization
  for (step in 1:max.p) {
    message.verbose(paste(date(), "Stepwise analysis, starting step", step), 0)
    if (stepwise.method == "forward") {
      ## reevaluate variable order by forcing already selected vars in the model
      o <- variable.ordering(data,group.labels,single.var=T,forced.vars=selected.vars, ...)
      var.order <- append(selected.vars, o$order.single.var)
      message.verbose(paste("step", step, "variable ordering", sep="        "),3)
      message.verbose(var.order,3)
    }
    selected.vars <- var.order[1:step]
    message.verbose(paste("step", step, "selected vars", sep= "        "), 3)
    message.verbose(selected.vars,3)

    ##    verbose (paste("","stepwise PDA",pda.method, max.p, step, names(data)[var.order[step]], sep="        "),3)
    data.reduced <- data[,selected.vars]
## THIS DOES NOT WORK ANYMORE NOM DE DIEU    data.reduced <- data.frame(data[,selected.vars])

    ## DEBUGGING TEMP
#    print(paste("selected vars", selected.vars))
#    print(dim(data.reduced))
#    print(length(group.labels))
#    print(length(test.labels))
#    print(pda.method)
#    print(data.reduced)

    next.pda <- pda(data.reduced,group.labels,test.labels, pda.method=pda.method,...)

#    stop("HELLO")

    result$errors[step] <- attr(next.pda$confusion,"error")
    conf.vect <- as.vector(next.pda$confusion)
    conf.row <- c(step,
                  result$errors[step],
                  selected.vars[step],
                  names(data)[selected.vars[step]],
                  conf.vect
                  )

    ## include the result of the confusion on separate test set
    if (!is.null(test.labels)) {
      result$errors.test[step] <- attr(next.pda$confusion.test,"error")
##      print (test.labels)
##      stop("BINGO")
##      conf.vect.test <- as.vector(next.pda$confusion.test)
##      conf.row <- c(conf.row,
##      result$errors.test[step],
##                       conf.vect.test
##                       )
#      result$confusion.table <- rbind(result$confusion.table,
#                                      c(step,
#                                        result$errors[step],
#                                        result$errors.test[step],
#                                        selected.vars[step],
#                                        names(data)[selected.vars[step]],
#                                        conf.vect
#                                        ))
    }
    ##    } else {
      result$confusion.table <- rbind(result$confusion.table,
                                      c(step,
                                        result$errors[step],
                                        selected.vars[step],
                                        names(data)[selected.vars[step]],
                                        conf.vect
                                        ))
#    }




    if ((last.min && (result$errors[step] <= error.min))  ### if the error rates are equal, take the lowest number of variables
	|| (!last.min && (result$errors[step] <= error.min))) { ### if the error rates are equal, take the highest number of variables, for robustness
      error.min <- result$errors[step]
      result$best.pda <- next.pda
      result$best.p <- step
      result$best.vars <- selected.vars
      result$best.var.names <- names(data)[selected.vars]
      result$best.error <- result$errors[step]
    }
    message.verbose(paste("step", step, "errors", sep="        "),3)
    message.verbose(result$errors,3)
  }

  ## set the names of the error table
  result$confusion.table <- as.data.frame(result$confusion.table)
  result$confusion.table[3] <- var.order[1:max.p]
  names <- vector()
  names <- c(names,"step")
  names <- c(names,"error")
  names <- c(names,"added.var")
  names <- c(names,"added.var.name")
  for (k in groups) {
    for (p in groups) {
      names <- append (names, paste(k, p,sep="."))
    }
  }
  names(result$confusion.table) <- names
  result$pc.transform <- pc.transform
#  result$errors <- errors

  ## Variable ordering: the sequence of numbers indicate the order of
  ## the of variables (each variable is specified by its column number
  ## in the data table)
  result$var.order <- var.order

  ## Variable rank: indicates the rank of each variable in the sorted list
  result$var.rank <- order(var.order)

  ## plot the error curve
  main <- paste(stepwise.method,  "stepwise ", toupper(pda.method), " - Error rate",sep = " ")
  if (!is.null(data.title)) {
    main <- paste(data.title, main, sep=" - ")
  }
  if (display.plot) {
    plot(1:length(result$errors),
         result$errors,
         col="#FF0000",
         main=main,
         xlab="number of components",
         ylab="Error rate",
         ylim=c(0,max(na.omit(result$errors))),
         type="b",,
         lwd=2,
         pch=7,
         panel.first = c(
           abline(h=0:1,col=1),
           abline(v=0:1,col=1),
           abline(h=seq(from=0, to=1, by=0.02), col="#EEEEEE"),
           abline(h=seq(from=0, to=1, by=0.1), col="darkgrey"),
           grid(lty="solid", col="darkgrey")
           )
         )
    ##	 text(1:max.p,result$errors,labels=names(data)[var.order],pos=3)
  }

  return(result)
}


################################################################
## Leave-one-out test outside the variable selection
## Warning: this is statistically unbiased, but computationally very
## inefficient: the stepwise discriminant analysis is performed N
## times succcessively.
##
## Select variables on the basis of the average rank.
##
loo.outside.stepwise.pda <- function(x, # A data frame with the objects to discriminate
                                     grouping, # a vector specifying the class for each observation
                                     export.interm.results=F, # export each intermediate result
                                     max.p=ncol(x), # Maximal number of variables to select
                                     ... # Arguments are passed to stepwise.pda
                                     ) {

  ## filter out objects with undefined label
  no.na <- !is.na(grouping)
  x <- x[no.na,]
  grouping <- grouping[no.na]

  ## dimensions

  ## Number of objects (rows)
  n <- nrow(x)
  message.verbose(paste("Number of labelled objects=", n))

  ## Number of variables (columns)
  p <- ncol(x)
  message.verbose(paste("Number of variables=", p))

  ## Number of groups
  groups <- unique(grouping)
  g <- length(groups)

  ## Initialize variables
  i <- 1 ## Iterator
  result <- list() ## the result object returned at the end o the procedure

  prediction <- matrix(ncol=g,nrow=n)
  prediction <- data.frame(prediction)

  ## Table to store the variable ordering, for each iteration of the leave-one-out
  var.order.table <- data.frame(matrix(nrow=n,ncol=p),row.names=row.names(x))
  names(var.order.table) <- names(x)

  ## Table to store the variable ranking
  var.rank.table <- data.frame(matrix(nrow=n,ncol=p),row.names=row.names(x))
  names(var.rank.table) <- names(x)

  ## error table
  ## one row per element (the left out element), and one column per selected variable
  errors <- data.frame(matrix(nrow=n,ncol=max.p),row.names=row.names(x))
  names(errors) <- 1:max.p

  ## confusion table
  conf.table <- data.frame()

  ## selected variables
  ## One row per object, one column per variable
  ## Boolean matrix indicating, for each object (row) which variables
  ## (columns) were selected or not when this object was left out
  selected.vars <- data.frame(matrix(0, nrow=n,ncol=p))
  names(selected.vars) <- names(x)
  row.names(selected.vars) <- row.names(x)

  pred.class <-rep(NA, length(grouping))
  names(pred.class) <- names(grouping)

  ## Iterate over the objects
  i <- 1
  for (i in 1:n) {
    print (paste("loo.outside.stepwise.pda, object", i,"of",n))

    ## Build a training set with all entries except the ith one
    train.labels <- grouping
    train.labels[i] <- NA

    ## Build a testing set with the ith entry
    test.labels <- rep(NA,n)
    test.labels[i] <- grouping[i]

    ## select variables without the ith object and test the prediction with this object
    spda <- stepwise.pda(x,train.labels,test.labels=test.labels, max.p=max.p,...)
    var.order.table[i,] <- spda$var.order
    var.rank.table[i,] <- spda$var.rank

    errors[i,] <- spda$errors
#    print(errors)

    conf.table <- rbind(conf.table, data.frame(i,
                                               leftout=row.names(x)[i],
                                               group=grouping[i],
                                               best.p=spda$best.p,
                                               best.err=spda$best.err,
                                               conf=t(data.frame(as.vector(as.matrix(spda$best.pda$confusion))))
                                        #                                         best.vars=t(data.frame(spda$best.vars))
                                         ))
#    print(conf.table)

    leftout <- row.names(x[i,])
    message.verbose(paste("Variable selection, discarding object", i, leftout))


#    mpda <- pda(x[-i,],grouping[-i])

    ## select variables
    selected <- spda$best.vars
    selected.vars[i,selected] <- 1

    ################################################################
    ## Predictive LDA : calculate the posterior probability for the left out object,
    ## with the predictor trained on all the other objects, and the selected variables
    spda.pred <- pda(x[-i,selected],grouping[-i],loo=F,pda.method="lda")  ## Build the predictor
    pred <- predict(spda.pred,data.frame(x[,selected]))

    ## build a predictive function
    prediction[i,] <- pred$posterior[i,]
    if (i==1) {
      names(prediction) <- names(data.frame(pred$posterior))
      row.names(prediction) <- row.names(data.frame(pred$posterior))
    }
    pred.class[i] <- as.vector(pred$class)[i]
##    print (prediction)
    if (export.interm.results) {
      tmp.file.name <- paste("tmp","LOO_outside",sep="_")
      export.object(cbind(prediction,pred.class,grouping),tmp.file.name,export.formats="table")
    }

  }

  message.verbose("LOO iterations done", 1)

  prediction$pred.class <- pred.class
  prediction$grouping <- grouping
                                        #    print (pred$posterior)



  result$prediction <- prediction
#  result$grouping <- grouping
#  result$pred.class <- pred.class
  result$confusion <- confusion(pred.class,grouping)

  result$var.order.table <- var.order.table

  result$var.rank.table <- var.rank.table
  result$var.rank.mean <- apply(var.rank.table, 2, mean)
  result$var.rank.sd <- apply(var.rank.table, 2, sd)

  result$errors <- errors
  result$conf.table <- conf.table
  result$selected.vars <- selected.vars
  result$selected.vars.freq <- apply(selected.vars, 2,sum)
  result$selected.vars.per.obj <- apply(selected.vars, 1,sum)

  message.verbose("Finished LOO outside stepwise variable selection", 1)

  return(result)
}



################################################################
##
## K-fold cross validation with stepwise discriminant analysis
##
## Select variables on the basis of the average rank.
##
k.fold.cv.stepwise.pda <- function(x, # A data frame with the objects to discriminate
                                     grouping, # a vector specifying the class for each observation
                                     k=10, # Number of groups for the K-fold cross-validation
                                     export.interm.results=F, # export each intermediate result
                                     max.p=ncol(x), # Maximal number of variables to select
                                     ... # Arguments are passed to stepwise.pda
                                     ) {


  ## filter out objects with undefined label
  no.na <- !is.na(grouping)
  x <- x[no.na,]
  grouping <- grouping[no.na]

  ## dimensions

  ## Number of objects (rows)
  n <- nrow(x)
  message.verbose(paste("Number of labelled objects=", n))

  ## Number of variables (columns)
  p <- ncol(x)
  message.verbose(paste("Number of variables=", p))

  ## Number of groups
  groups <- unique(grouping)

  ## Initialize variables
  result <- list() ## the result object returned at the end of the procedure

  prediction <- matrix(ncol=g,nrow=n)
  prediction <- data.frame(prediction)

  ## Table to store the variable ordering, for each iteration of the leave-one-out
  var.order.table <- data.frame(matrix(nrow=n,ncol=p),row.names=row.names(x))
  names(var.order.table) <- names(x)

  ## Table to store the variable ranking
  var.rank.table <- data.frame(matrix(nrow=n,ncol=p),row.names=row.names(x))
  names(var.rank.table) <- names(x)

  ## error table
  ## one row per element (the left out element), and one column per selected variable
  errors <- data.frame(matrix(nrow=n,ncol=max.p),row.names=row.names(x))
  names(errors) <- 1:max.p

  ## confusion table
  conf.table <- data.frame()

  ## selected variables
  ## One row per object, one column per variable
  ## Boolean matrix indicating, for each object (row) which variables
  ## (columns) were selected or not when this object was left out
  selected.vars <- data.frame(matrix(0, nrow=n,ncol=p))
  names(selected.vars) <- names(x)
  row.names(selected.vars) <- row.names(x)

  pred.class <-rep(NA, length(grouping))
  names(pred.class) <- names(grouping)

  ## Assign each training object to one of the K validation groups
  validation.group.breaks <- c(1, round((1:k)*n/k)) ## Number of objects per group
  random.order <- sample(1:n)
  validation.group <- sample(1:10, size=n,rep=T)

  ## Iterate over the k groups
  g <- 1 ## Iterator
  for (g in 1:k) {

    ## test objects
    i <- random.order[validation.group.breaks[g]:validation.group.breaks[g+1]]

    message.verbose(paste("k.fold.cv.stepwise.pda group", g, "; ", length(i), "objects of",n), 1)

    ## Build a training set with all entries except the ith one
    train.labels <- grouping
    train.labels[i] <- NA
    table(train.labels)

    ## Build a testing set with the ith entry
    test.labels <- rep(NA,n)
    test.labels[i] <- grouping[i]
    table(test.labels)

    ## select variables without the ith object and test the prediction with this object
    spda <- stepwise.pda(x,train.labels,test.labels=test.labels, max.p=max.p,...)
    var.order.table[i,] <- spda$var.order
    var.rank.table[i,] <- spda$var.rank

    errors[i,] <- spda$errors
#    print(errors)

    conf.table <- rbind(conf.table,
                        data.frame(i,
                                   leftout=row.names(x)[i],
                                   group=grouping[i],
                                   best.p=rep(spda$best.p, length(i)),
                                   best.err=rep(spda$best.err, length(i)),
                                   conf=t(matrix(rep(t(data.frame(as.vector(as.matrix(spda$best.pda$confusion)))), length(i)),ncol=length(i)))
                                        #                                         best.vars=t(data.frame(spda$best.vars))
                                   )
                        )

                                        #    print(conf.table)

    leftout <- row.names(x[i,])
    message.verbose(paste("k=", k, "; Variable selection, discarding ", length(i), " objects", sep=""))


#    mpda <- pda(x[-i,],grouping[-i])

    ## select variables
    selected <- spda$best.vars
    selected.vars[i,selected] <- 1

    ################################################################
    ## Predictive LDA : calculate the posterior probability for the left out object,
    ## with the predictor trained on all the other objects, and the selected variables
    spda.pred <- pda(x[-i,selected],grouping[-i],loo=F,pda.method="lda")  ## Build the predictor
    pred <- predict(spda.pred,data.frame(x[,selected]))

    ## build a predictive function
    prediction[i,] <- pred$posterior[i,]
    if (i==1) {
      names(prediction) <- names(data.frame(pred$posterior))
      row.names(prediction) <- row.names(data.frame(pred$posterior))
    }
    pred.class[i] <- as.vector(pred$class)[i]
##    print (prediction)
    if (export.interm.results) {
      tmp.file.name <- paste("tmp",k, "fold_CV",sep="_")
      export.object(cbind(prediction,pred.class,grouping),tmp.file.name,export.formats="table")
    }

  }

  message.verbose("CV iterations done", 1)

  prediction$pred.class <- pred.class
  prediction$grouping <- grouping
                                        #    print (pred$posterior)



  result$prediction <- prediction
#  result$grouping <- grouping
#  result$pred.class <- pred.class
  result$confusion <- confusion(pred.class,grouping)

  result$var.order.table <- var.order.table

  result$var.rank.table <- var.rank.table
  result$var.rank.mean <- apply(var.rank.table, 2, mean)
  result$var.rank.sd <- apply(var.rank.table, 2, sd)

  result$errors <- errors
  result$conf.table <- conf.table
  result$selected.vars <- selected.vars
  result$selected.vars.freq <- apply(selected.vars, 2,sum)
  result$selected.vars.per.obj <- apply(selected.vars, 1,sum)

  message.verbose("Finished K-fold CV with stepwise variable selection", 1)

  return(result)
}


## ##############################################################
##
## Iterative predictive linear analysis
## First cycle:
## - perform a predictive linear analysis
## At each new cycle :
## - detect the misclassified units from the previous cycle
## - remove these misclassified units from the training set
## - perform a predictive linear analysis
##
iterative.stepwise.pda <- function (data,
				    group.labels,
				    max.cycles=2, # maximum number of cycles for the iteration
				    performed.cycles=0, # count of the cycles already performed
                                    export.formats.plots=NA, # export results at each cycle
                                    export.formats.obj=NA, # export results at each cycle
				    data.type = "data", # prefix for saving files
				    ...
				    ) {
  ## Data title for stepwise.sda()
  data.title <- paste(data.type, "iterative", sep="_")

  ## verbose
  verbose ("Iterative stepwise PDA")
  verbose (paste("Performed cycles : ", performed.cycles, "; Remaining cycles : ", max.cycles))

  ## check number of cycles
  if (max.cycles < 1) {
    stop ("cycles should be >=1")
  }

  s <- stepwise.pda (data,group.labels, data.title=data.title, ...)
  performed.cycles <- performed.cycles + 1
  s$miscl <- misclassified.units(data,group.labels,s$best.pda)
  group.labels.filtered <- group.labels
  group.labels.filtered[row.names(s$miscl)] <- NA
  s$group.labels.filtered <- group.labels.filtered
  s$max.cycles <- max.cycles
  s$performed.cycles <- performed.cycles


  ## Export the resulting object and graphic
  file.name <- paste(data.type,"iterative_sda_cycle",performed.cycles, sep="_")
  if (!is.na(export.formats.obj)) {
    setwd(dir.results); export.object(s, file=file.name, export.formats=export.formats.obj)
  }
  if (!is.na(export.formats.plots)) {
    setwd(dir.figures); export.plot(file.prefix=file.name, export.formats=export.formats.plots)
  }
#   file.name <- paste(data.type,"iterative_sda_cycle",performed.cycles, sep="_")
#   export.object(s, file=file.name, export.formats=export.formats.obj)
#   export.plot(file.prefix=file.name, export.formats=export.formats.plots)


  if (sum(is.na(group.labels.filtered) & !is.na(group.labels)) == 0) {
#  if (group.labels.filtered == group.labels) {
    message.verbose("iteration converged because group labels are not modified anymore", 0)
    return (s)
  } else if (max.cycles == 1) {
    return (s)
  } else {
    ## next cycle
    return (iterative.stepwise.pda(data,group.labels.filtered,
                                   max.cycles=max.cycles-1,
                                   performed.cycles=performed.cycles,
                                   export.formats.plots=export.formats.plots,
                                   export.formats.obj=export.formats.obj,
                                   data.type=data.type,
                                   ...
                                   ))
  }

}



## ##############################################################
## Export profiles of prior classes, posterior classes, and
## misclassified units
##
export.sda.plots <- function(data,
                     	     group.labels,
			     export.dir=".",
			     file.prefix="",
			     data.title="",
                             princomp=F,
			     export.formats.plots = c("postscript","jpg"),
			     export.formats.obj = NULL,
                             lda=T,
                             qda=T,
			     ...
			     ) {
  wd <- getwd()
  setwd(export.dir)



  if (princomp) {
    ## principal components with prcomp
    pc <- prcomp.with.group.labels(data,group.labels,data.title=data.title, display.plot=T)
    file.name <- paste(file.prefix,'PC_prcomp',sep='_')
    export.plot(file.prefix=file.name,export.formats=export.formats.plots)
    export.object(pc,file.prefix=file.name,export.formats=export.formats.obj)

    ## principal components with princomp
    pc2 <- princomp.with.group.labels(data,group.labels,data.title=data.title,display.plot=T)
    file.name <- paste(file.prefix,'PC_princomp',sep='_')
    export.plot(file.prefix=file.name,export.formats=export.formats.plots)
    export.object(pc2,file.prefix=file.name,export.formats=export.formats.obj)

    data <- pc2$scores
    file.prefix <- paste(file.prefix,"PC",sep="_")
    data.title <- paste(data.title," (PC)")
  }

  data.title.ori <- data.title
  ## linear or quadratic discriminant analysis
  methods <- vector()
  if (lda) {methods <- append(methods, "lda") }
  if (qda) {methods <- append(methods, "qda") }

  ##  for (method in c("lda", "qda")) {
  for (pda.method in methods) {
    if (pda.method == "qda") {
      data.reduced <- qda.check.variables(data, group.labels,...)

    } else {
      data.reduced <- data
    }

    ## LOO of internal validation
    for (loo in c(T, F)) {
      ## plot profiles by predicted class
      message.verbose("Profiles by predicted class", 2)
      plot.profiles.by.predicted.class(data.reduced,group.labels,data.title=data.title,...)
      file.name <- paste(file.prefix,'profile_predicted_class', pda.method, loo, sep='_')
      export.plot(file.prefix=file.name,
                  export.formats=export.formats.plots,
                  width=8,
                  height=10,
                  horizontal=F)

      ## plot misclassified units
      message.verbose("Misclassified units", 2)
      misclass <- plot.misclassified.units(data.reduced,group.labels,data.title=data.title,display.plot=T,...)
      file.name <- paste(file.prefix,'misclassified_profiles',pda.method,loo,sep='_')
      export.plot(file.prefix=file.name,
                  export.formats=export.formats.plots,
                  width=8,
                  height=10,
                  horizontal=F)
      export.object(misclass,file.prefix=file.name,export.formats=export.formats.obj)

      ## plot the calibration results (predicted class as character, predefined class as color)
      message.verbose("PDA calibration", 2)
      calib <- plot.pda.calib(data.reduced,group.labels,data.title=data.title,...)
      file.name <- paste(file.prefix,'pda_calib', pda.method, loo, sep='_')
      export.plot(file.prefix=file.name,
                  export.formats=export.formats.plots,
                  width=8,
                  height=10,
                  horizontal=F)
      export.object(calib,file.prefix=file.name,export.formats=export.formats.obj)
    }
    rm (data.reduced)
    data.title <- data.title.ori
  }

  setwd(wd)
}

################################################################
##
## Check the number of variables for quadratic discriminant analysis.
## If required, reduce the data to the max.p most discriminant
## variables. The selection is based on a single-variable qda, with
## LOO evaluation.
## usage:
##    data.checked <- check.variables(data,group.labels)
##
qda.check.variables <- function (data,
				 group.labels,
				 df=3,
				 max.p=NA,
				 data.type = 'data', ## Prefix for titles
                                 ...) {

  p <- dim(data)[2]
  ## check the number of variables
  if (is.na(max.p)) {
    max.p <- p
  }
                                        #  if (pda.method == "qda") {
  min.n <- min(table(group.labels)) # nb of elements in the smallest group
  max.p <- min(p,min.n-df)
                                        #  }

  ## data reduction
  if (max.p < p) {
    message.verbose(paste("Too many variables (",p, ").",
                  "Selecting", max.p, "variables"), 1)
    ## single-var based selection of the max.p best variables
    vo <- variable.ordering(data,group.labels, pda.method="qda", ...)
    selected.vars <- vo$order.single.var[1:max.p]
    data.reduced <- data[,selected.vars]
    data.type <- paste(data.type, 1, max.p, sep="-")
    data.title <- paste(data.title, "reduced to", max.p, "variables")
    message.verbose("selected variables", 1)
    message.verbose(selected.vars,1)
    return (data.reduced)
  } else {
    return(data)
  }
}



################################################################
## Calculate the hit rate for each variable separately.
variable.ordering <- function (data,               # explanatory variables
                               group.labels,	   # grouping variable
			       forced.vars=vector(), # force specified variables to be in the model
                               ## the ordering is performed on the other varibles
#                               pda.method="lda",	   # the method for assessing variable order
#                               loo=T,		   # use leave-one-out for assessment of error rate
                               ## when FALSE, error rate is assessed by internal validation
                               all.var=F,          # calculate confusion table for al variables together
                               single.var=T,	   # calculate error rate for each variable independently
                               ## (and return in first position the vars with smallest error)
                               deleted.var=F,	   # calculate error rate when each variable is deleted
                               ## (and return in first position the vars with highest error)
                               ... ## Additional arguments are passed to pda()
                               ) {

  result <- list()

  n <- dim(data)[1]; result$n <- n
  p <- dim(data)[2]; result$p <- p
  group.labels.training <- group.labels[!is.na(group.labels)]
  n.training <- length(group.labels.training); result$n.training <- n.training
  result$forced.vars <- forced.vars
  vars.to.test <- setdiff (1:p, forced.vars); result$tested.vars <- vars.to.test
#  result$pda.method <- pda.method
#  result$loo <- loo

  ## calculate the hit rate with all variables
  if (all.var) {
    pda.all.var<- pda(data,group.labels,conf=T, ...)
    result$confusion.all.var <- pda.all.var$confusion
  }

  ## calculate the hit rate with all variables
  if (length(forced.vars >0)) {
    pda.forced.var<- pda(data[,forced.vars],group.labels,conf=T, ...)
    result$confusion.forced.var <- pda.forced.var$confusion
    result$err.forced.var <- attr(pda.forced.var$confusion,"error")
  }

  ## calculate error rates with single variable analysis
  if (single.var) {
    err.single.var <- rep(length.out=p,9999)
    for (i in vars.to.test) {
#      pda.single.var<- pda(data[,union(forced.vars, i)],group.labels,conf=T, ...)
      pda.single.var<- pda(data[,union(forced.vars, i)],group.labels,conf=T)
      err.single.var[i] <- attr(pda.single.var$confusion, "error")
      message.verbose(paste("Variable ordering ", i, "/", length(vars.to.test), "; error=", err.single.var[i], sep=""), 5)
    }
    result$err.single.var <- err.single.var  # single variable hit rate
    message.verbose(paste("Variable ordering - errors per variable: ",paste(result$err.single.var), sep=""), 3)
    result$order.single.var <- setdiff(order(err.single.var),forced.vars)
    message.verbose(paste("Variable ordering - variable order: ",paste(result$order.single.var), sep=""), 3)
  }

  ## calculate error rates by individual variable deletion
  if (deleted.var) {
    err.deleted.var <- rep(length.out=p,-9999)
    for (i in vars.to.test) {
      pda.deleted.var<- pda(data[,-i],group.labels, ...)
      err.deleted.var[i] <- attr(pda.deleted.var$confusion, "error")
      message.verbose(paste("Variable ordering (deletion) ", i, "/", length(vars.to.test), "; error=", err.deleted.var[i], sep=""), 5)
    }
    result$err.deleted.var <- err.deleted.var  # deleted variable hit rate
    message.verbose(paste("Variable ordering - errors per deleted variable: ",paste(result$err.deleted.var), sep=""), 3)
    result$order.deleted.var <- setdiff(rev(order(err.deleted.var)),forced.vars)
    message.verbose(paste("Variable ordering - deleted variable order: ",paste(result$order.deleted.var), sep=""), 3)
  }

  ## return the result
  return(result)
}

################################################################
## Perform a permutation test on stepwise predictive discriminant
## analysis. The data set is used as such, but group labels are
## randomly permuted. The permutation test is repeated n times, and
## the program returns an object containing
## - a data frame with the error profile for each repetition
## - a mean of the error profiles
## - a standard deviation of the error profiles
##
stepwise.pda.permut.test <- function (data,              # a data frame
                                      group.labels,      # a character vector of length N, containing group labels
                                      max.p=NA,          # max number of variables
                                      test.number=10, # Number of reppetitions of the permutation test
                                      ...) {
  message.verbose(paste("stepwise PDA, ",test.number," permutation tests", sep=""), 1)
  groups <- names(table(group.labels))
  g <- length(table(group.labels))
  if (is.na(max.p)) {
    max.p <- ncol(data)
  }
  result <- list()
  result$confusion.table <- matrix(ncol=g^2+1,nrow=0)
  result$errors <- data.frame(matrix(nrow=max.p,ncol=test.number+2))
  names(result$errors) <- c(paste("permut",1:test.number, sep="."), "err.mean", "err.sd")
  for (i in 1:test.number) {
    message.verbose(paste("stepwise PDA, permutation test ", i, "/", test.number,sep=""), 1)
    group.labels.permuted <- rep(NA, length(group.labels))
    group.labels.permuted[!is.na(group.labels)] <- sample(group.labels[!is.na(group.labels)])
    one.permut.test <- stepwise.pda(data,group.labels.permuted,max.p=max.p,...)
    result$confusion.table <- rbind(result$confusion.table,
                                    c(attr(one.permut.test$best.pda$confusion,"error"),
                                      as.vector(as.matrix(one.permut.test$best.pda$confusion))))
    result$errors[1:length(one.permut.test$errors),i] <- one.permut.test$errors
    result$best.p <- one.permut.test$best.p
    result$best.error <- one.permut.test$best.error
    result$best.pda <- one.permut.test$best.pda
 }

  message.verbose(paste("stepwise PDA, permutation tests: confusion table", sep=""), 1)
  result$confusion.table <- data.frame(result$confusion.table)
  conf.names <- vector()
  conf.names <- c(conf.names, "error")
  for (k in groups) {
    for (p in groups) {
      conf.names <- append (conf.names, paste(k, p,sep="."))
    }
  }
  names(result$confusion.table) <- conf.names

  message.verbose(paste("stepwise PDA, permutation tests: calculating errors", sep=""), 1)
  if (test.number > 1) {
    result$errors[,"err.mean"] <- apply(result$errors[,1:test.number],1,mean,na.rm=T)
    result$errors[,"err.sd"] <- apply(result$errors[,1:test.number],1,sd,na.rm=T)
    result$errors[,"err.sterr"] <-     result$errors[,"err.sd"]/test.number^0.5
  } else {
    result$errors[,"err.mean"] <- result$errors
  }
  message.verbose(paste("stepwise PDA, permutation tests: done", sep=""), 1)
  return(result)
}

################################################################
##
## Compare the error rate curve obtained by stepwise
## predictive discriminant analysis, with different methods
## - LDA and QDA
## - apply or not PCA before analysis
## - same analysis on random data
##
################################################################
compare.stepwise.pda.methods <- function(data,
					 group.labels,
#					 test.labels=NULL,
                                         data.title="",
					 display.plots=T, # display intermediate plots
					 export.result=F, # export intermediate results
					 export.formats.plots=c("postscript","jpg"),
					 export.formats.obj=c("print"),
                                         lda=T, ### test lindear discriminant analysis
                                         qda=T, ### test quadratic discriminant analysis
					 permut.test=T, ### add a test with permuted group labels
                                         test.number=10, ## Number of permutation tests
					 rand.test=F, ### add a test with multivariate normal random data
					 without.pc=T,   ### apply tests without PC transformation
					 with.pc=T,   ### apply tests after PC transformation
					 max.p=NA,
					 data.type = 'data', ## Prefix for graphic titles, ...
					 ...
					 ) {
  if ((!lda) && (!qda)) {
    stop( "Either lda or qda must be TRUE")
  }

  data <- data.frame(data)
  n <- nrow(data)
  p <- ncol(data)
  if (is.na(max.p)) {
    max.p <- p
  }
  data.type.ori <- data.type
  data.title.ori <- data.title

  ## columns in the result
  ncol <- 0

  ## stepwise PDA
  message.verbose("stepwise PDA", 1)

  if (without.pc) {
    message.verbose("stepwise PDA with untransformed data",1)
    if (lda) {
      (stepwise.lda <- stepwise.pda(data,group.labels,pda.method="lda",max.p=max.p, display.plot=display.plots,data.title=data.title,...))
      ncol <- ncol + 1
      if (export.result) {
        file.name <- paste(data.type, "LDA_LOO", sep="_")
        export.plot(file.prefix=file.name, export.formats=export.formats.plots)
        export.object(stepwise.lda,file=file.name,export.formats=export.formats.obj)
      }
    }

    if (qda) {
      (stepwise.qda <- stepwise.pda(data,group.labels,pda.method="qda",max.p=max.p, display.plot=display.plots,data.title=data.title,...))
      ncol <- ncol + 1
      if (export.result) {
        file.name <- paste(data.type, "QDA_LOO", sep="_")
        export.plot(file.prefix=file.name, export.formats=export.formats.plots)
        export.object(stepwise.qda,file=file.name,export.formats=export.formats.obj)
      }
    }

  }

  if (with.pc) {

    ## stepwise PDA with principal components
    message.verbose("stepwise PDA with principal components",1)
    pc.data <- princomp.with.group.labels(data,group.labels,display.plot=F)
    data.type <- paste(data.type.ori, "PC",sep="_")
    data.title <- paste (data.title.ori, " (PC)", sep="")

    if  (lda) {
      (pc.stepwise.lda <- stepwise.pda(pc.data$scores,group.labels,pda.method="lda",max.p=max.p, display.plot=display.plots,data.title=data.title, pc.transform=T,...))
      ncol <- ncol + 1
      if (export.result) {
        file.name <- paste(data.type, "LDA_LOO", sep="_")
        export.plot(file.prefix=file.name, export.formats=export.formats.plots)
        export.object(pc.stepwise.lda,file=file.name,export.formats=export.formats.obj)
      }
    }

    if (qda) {
      (pc.stepwise.qda <- stepwise.pda(pc.data$scores,group.labels,pda.method="qda",max.p=max.p, display.plot=display.plots,data.title=data.title,pc.transform=T,...))
      ncol <- ncol + 1
      if (export.result) {
        file.name <- paste(data.type, "QDA_LOO", sep="_")
        export.plot(file.prefix=file.name, export.formats=export.formats.plots)
        export.object(pc.stepwise.qda,file=file.name,export.formats=export.formats.obj)
      }
    }
  }

  ## ##############################################################
  ## permutation test : stepwise PDA with the same data set but
  ## permuted labels
  if (permut.test) {
    message.verbose("stepwise PDA with permuted labels",1)
    data.type <- paste(data.type.ori, "permut", sep="_")
    data.title <- paste ("permuted", data.title.ori)
    permut.labels <- sample(group.labels)

    if (without.pc) {
      if (lda) {
#        message.verbose("stepwise PDA with permuted labels, LDA, no PC",1)
        (permut.stepwise.lda <- stepwise.pda.permut.test(data,group.labels,test.number=test.number,pda.method="lda",max.p=max.p, display.plot=display.plots,data.title=data.title,...))
#        (permut.stepwise.lda <- stepwise.pda(data,permut.labels,pda.method="lda",max.p=max.p, display.plot=display.plots,data.title=data.title,...))
	ncol <- ncol + 1
        if (export.result) {
          file.name <- paste(data.type, "LDA_LOO", sep="_")
          export.plot(file.prefix=file.name, export.formats=export.formats.plots)
          export.object(permut.stepwise.lda,file=file.name,export.formats=export.formats.obj)
        }
      }

      if (qda) {
        (permut.stepwise.qda <- stepwise.pda.permut.test(data,group.labels,test.number=test.number,pda.method="qda",max.p=max.p, display.plot=display.plots,data.title=data.title,...))
#        (permut.stepwise.qda <- stepwise.pda(data,permut.labels,pda.method="qda",max.p=max.p, display.plot=display.plots,data.title=data.title,...))
        ncol <- ncol + 1
        if (export.result) {
          file.name <- paste(data.type, "QDA_LOO", sep="_")
          export.plot(file.prefix=file.name, export.formats=export.formats.plots)
          export.object(permut.stepwise.qda,file=file.name,export.formats=export.formats.obj)
        }
        }
    }
#    }

    if (with.pc) {

      ## stepwise PDA with PCA with permuted labels
      message.verbose("stepwise PDA with PCA with permuted labels", 1)
      data.type <- paste(data.type.ori, "permut", "PC", sep="_")
      data.title <- paste ("permuted", data.title.ori, "(PC)")

      if (lda) {
        (permut.pc.stepwise.lda <- stepwise.pda.permut.test(pc.data$scores,group.labels,test.number=test.number,pda.method="lda",max.p=max.p, display.plot=display.plots,data.title=data.title,pc.transform=T,...))
#        (permut.pc.stepwise.lda <- stepwise.pda(pc.data$scores,permut.labels,pda.method="lda",max.p=max.p, display.plot=display.plots,data.title=data.title,pc.transform=T,...))
        ncol <- ncol + 1
        if (export.result) {
          file.name <- paste(data.type, "LDA_LOO", sep="_")
          export.plot(file.prefix=file.name, export.formats=export.formats.plots)
          export.object(permut.pc.stepwise.lda,file=file.name,export.formats=export.formats.obj)
        }
      }

      if (qda) {
        (permut.pc.stepwise.qda <- stepwise.pda.permut.test(pc.data$scores,group.labels,test.number=test.number,pda.method="qda",max.p=max.p, display.plot=display.plots,data.title=data.title,pc.transform=T,...))
#        (permut.pc.stepwise.qda <- stepwise.pda(pc.data$scores,permut.labels,pda.method="qda",max.p=max.p, display.plot=display.plots,data.title=data.title,pc.transform=T,...))
        ncol <- ncol + 1
        if (export.result) {
          file.name <- paste(data.type, "QDA_LOO", sep="_")
          export.plot(file.prefix=file.name, export.formats=export.formats.plots)
          export.object(permut.pc.stepwise.qda,file=file.name,export.formats=export.formats.obj)
        }
      }
    }
  }

  ## stepwise PDA with random data
  if (rand.test) {
    message.verbose("stepwise PDA with random data",1)
    rand <- randomize(data)
    data.type <- paste(data.type.ori, "rand", sep="_")
    data.title <- paste ("random", data.title.ori)


    if (without.pc) {
      if (lda) {
        (rand.stepwise.lda <- stepwise.pda(rand,group.labels,pda.method="lda",max.p=max.p, display.plot=display.plots,data.title=data.title,...))
	ncol <- ncol + 1
        if (export.result) {
          file.name <- paste(data.type, "LDA_LOO", sep="_")
          export.plot(file.prefix=file.name, export.formats=export.formats.plots)
          export.object(rand.stepwise.lda,file=file.name,export.formats=export.formats.obj)
       }
      }

      if (qda) {
       (rand.stepwise.qda <- stepwise.pda(rand,group.labels,pda.method="qda",max.p=max.p, display.plot=display.plots,data.title=data.title,...))
        ncol <- ncol + 1
        if (export.result) {
          file.name <- paste(data.type, "QDA_LOO", sep="_")
          export.plot(file.prefix=file.name, export.formats=export.formats.plots)
          export.object(rand.stepwise.qda,file=file.name,export.formats=export.formats.obj)
        }
      }
    }

    if (with.pc) {

      ## stepwise PDA with PCA of random data
      message.verbose("stepwise PDA with PCA of random data", 1)
      pc.rand <- princomp.with.group.labels(rand,group.labels,display.plot=F)
      data.type <- paste(data.type.ori, "rand", "PC", sep="_")
      data.title <- paste ("random", data.title.ori, "(PC)")

      if (lda) {
        (rand.pc.stepwise.lda <- stepwise.pda(pc.rand$scores,group.labels,pda.method="lda",max.p=max.p, display.plot=display.plots,data.title=data.title,pc.transform=T,...))
        ncol <- ncol + 1
        if (export.result) {
          file.name <- paste(data.type, "LDA_LOO", sep="_")
          export.plot(file.prefix=file.name, export.formats=export.formats.plots)
          export.object(rand.pc.stepwise.lda,file=file.name,export.formats=export.formats.obj)
        }
      }

      if (qda) {
        (rand.pc.stepwise.qda <- stepwise.pda(pc.rand$scores,group.labels,pda.method="qda",max.p=max.p, display.plot=display.plots,data.title=data.title,pc.transform=T,...))
        ncol <- ncol + 1
        if (export.result) {
          file.name <- paste(data.type, "QDA_LOO", sep="_")
          export.plot(file.prefix=file.name, export.formats=export.formats.plots)
          export.object(rand.pc.stepwise.qda,file=file.name,export.formats=export.formats.obj)
        }
      }
    }
  }


   ## ##############################################################
   ## Synthetize the results
   ##
   legend.labels <- vector()
   errors <- matrix(nrow=p,ncol=ncol)
   pca.methods <- vector()
   result <- list()   ## create the result object

   col <- 1
   if (lda) {
     if (without.pc) {
       errors[1:length(stepwise.lda$errors),col] <- stepwise.lda$errors; col <- col+1
       legend.labels <- append (legend.labels, "LDA")
       pca.methods <- append (pca.methods, "stepwise.lda")
       result$stepwise.lda <- stepwise.lda
     }
     if (with.pc) {
       errors[1:length(pc.stepwise.lda$errors),col] <- pc.stepwise.lda$errors; col <- col+1
       legend.labels <- append (legend.labels, "LDA with PC")
       pca.methods <- append (pca.methods, "pc.stepwise.lda")
       result$pc.stepwise.lda <- pc.stepwise.lda
     }
     if (permut.test) {
       if (without.pc) {
         errors[1:length(permut.stepwise.lda$errors$err.mean),col] <- permut.stepwise.lda$errors$err.mean; col <- col+1
#         errors[1:length(permut.stepwise.lda$errors),col] <- permut.stepwise.lda$errors; col <- col+1
         legend.labels <- append (legend.labels, "LDA, permuted")
         pca.methods <- append (pca.methods, "permut.stepwise.lda")
         result$permut.stepwise.lda <- permut.stepwise.lda
       }
       if (with.pc) {
         errors[1:length(permut.pc.stepwise.lda$errors$err.mean),col] <- permut.pc.stepwise.lda$errors$err.mean; col <- col+1
#         errors[1:length(permut.pc.stepwise.lda$errors),col] <- permut.pc.stepwise.lda$errors; col <- col+1
         legend.labels <- append (legend.labels, "LDA with PC, permuted")
         pca.methods <- append (pca.methods, "permut.pc.stepwise.lda")
         result$permut.pc.stepwise.lda <- permut.pc.stepwise.lda
       }
     }
     if (rand.test) {
       if (without.pc) {
         errors[1:length(rand.stepwise.lda$errors),col] <- rand.stepwise.lda$errors; col <- col+1
         legend.labels <- append (legend.labels, "LDA, random data")
         pca.methods <- append (pca.methods, "rand.stepwise.lda")
         result$rand.stepwise.lda <- rand.stepwise.lda
       }
       if (with.pc) {
         errors[1:length(rand.pc.stepwise.lda$errors),col] <- rand.pc.stepwise.lda$errors; col <- col+1
         legend.labels <- append (legend.labels, "LDA with PC, random data")
         pca.methods <- append (pca.methods, "rand.pc.stepwise.lda")
         result$rand.pc.stepwise.lda <- rand.pc.stepwise.lda
       }
     }
   }
  if (qda) {
    if (without.pc) {
      errors[1:length(stepwise.qda$errors),col] <- stepwise.qda$errors; col <- col+1
      legend.labels <- append (legend.labels, "QDA")
      pca.methods <- append (pca.methods, "stepwise.qda")
      result$stepwise.qda <- stepwise.qda
    }
    if (with.pc) {
      errors[1:length(pc.stepwise.qda$errors),col] <- pc.stepwise.qda$errors; col <- col+1
      legend.labels <- append (legend.labels, "QDA with PC")
      pca.methods <- append (pca.methods, "pc.stepwise.qda")
      result$pc.stepwise.qda <- pc.stepwise.qda
    }
    if (permut.test) {
      if (without.pc) {
        errors[1:length(permut.stepwise.qda$errors$err.mean),col] <- permut.stepwise.qda$errors$err.mean; col <- col+1
#        errors[1:length(permut.stepwise.qda$errors),col] <- permut.stepwise.qda$errors; col <- col+1
        legend.labels <- append (legend.labels, "QDA, permuted")
        pca.methods <- append (pca.methods, "permut.stepwise.qda")
        result$permut.stepwise.qda <- permut.stepwise.qda
      }
      if (with.pc) {
        errors[1:length(permut.pc.stepwise.qda$errors$err.mean),col] <- permut.pc.stepwise.qda$errors$err.mean; col <- col+1
#        errors[1:length(permut.pc.stepwise.qda$errors),col] <- permut.pc.stepwise.qda$errors; col <- col+1
        legend.labels <- append (legend.labels, "QDA with PC, permuted")
        pca.methods <- append (pca.methods, "permut.pc.stepwise.qda")
        result$permut.pc.stepwise.qda <- permut.pc.stepwise.qda
      }
    }
    if (rand.test) {
      if (without.pc) {
        errors[1:length(rand.stepwise.qda$errors),col] <- rand.stepwise.qda$errors; col <- col+1
        legend.labels <- append (legend.labels, "QDA, random data")
        pca.methods <- append (pca.methods, "rand.stepwise.qda")
        result$rand.stepwise.qda <- rand.stepwise.qda
      }
      if (with.pc) {
        errors[1:length(rand.pc.stepwise.qda$errors),col] <- rand.pc.stepwise.qda$errors; col <- col+1
        legend.labels <- append (legend.labels, "QDA with PC, random data")
        pca.methods <- append (pca.methods, "rand.pc.stepwise.qda")
        result$rand.pc.stepwise.qda <- rand.pc.stepwise.qda
      }
    }
  }

  message.verbose("pca.methods",1)
  message.verbose(paste (pca.methods, sep=";"),1)
  errors <- as.data.frame(errors)
  names(errors) <- pca.methods


  ## restore parameters
  data.type <- data.type.ori
  data.title <- data.title.ori

  plot.error.profiles(errors,main=paste(data.title,"- Stepwise PDA - Error rates"),legend.labels=legend.labels)

  if (test.number > 1) {
    if (lda) {
      for (i in 1:test.number) {
        lines(row.names(permut.stepwise.lda$errors), permut.stepwise.lda$errors[,i],col="#FFCCCC")
      }
    }
    if (qda) {
      for (i in 1:test.number) {
        lines(row.names(permut.stepwise.qda$errors), permut.stepwise.qda$errors[,i],col="#CCCCCC")
      }
    }
  }

  #################################################################
  ## Generate a summary confusion table with the best pdas
  verbose ("Summarizing results", 2)
  groups <- levels(as.factor(group.labels))
  g <- length(groups)
  confusion.table <- matrix(nrow=0,ncol=g^2 + 2)

  for (m in pca.methods) {
    method <- get(m)
    confusion.table <- rbind(confusion.table,
                             c(
                               method$best.p,
                               method$best.error,
                               as.vector(method$best.pda$confusion)
                               )
                             )
  }
  confusion.table <- as.data.frame(confusion.table)
  names <- vector()
  names <- c(names, "p")
  names <- c(names, "error")
  for (k in groups) {
    for (p in groups) {
      names <- append (names, paste(k, p,sep="."))
    }
  }

  names(confusion.table) <- names
  row.names(confusion.table) <- pca.methods


  ## compare the results obtained with the different approaches
  best.error <- 1
  for (m in pca.methods) {
    method <- get(m)
    if (method$best.error < best.error) {
      best.error <- method$best.error
      best.vars <- method$best.vars
      best.method <- m
      best.pda <- method
    }
  }

  result$max.p <- max.p
  result$errors <- errors
  result$legend.labels <- legend.labels
  result$confusion.table <- confusion.table
  result$best.errors <- best.error
  result$best.vars <- best.vars
  result$best.method <- best.method
  result$best.pda <- best.pda
  return(result)

}

################################################################
#
# Plot error rates returned by compare.stepwise.pda.methods(...)
#
plot.error.profiles <- function (errors=NA,
                                 main="error profiles",
                                 legend.labels=NA,
				 legend.x = NA,
				 legend.y = NA,
                                 plot.legend=T,
                                 col=NA,
				 x.values = NA,
				 ...
                                 ) {

  ## legend labels
  if (is.na(legend.labels)) legend.labels <- names(errors)

  ## curve colors
  if (is.na(col)) {col <- palette}
  if (length(col) < ncol(errors)) {
    col <- rep(col, length.out=ncol(errors))
  }


  verbose ("Plotting error rates", 2)
  par(font=2)
  par(font.lab=2)

  if (is.na(x.values)) {
    x.values <- 1:length(errors[,1])
  }

  if (plot.legend) {
     xlim <- c(1,max(x.values)*1.25)
  } else {
     xlim <- c(1,max(x.values))
  }
  plot(x.values,
       errors[,1],
       ylim=c(0,max(errors,na.rm=T)),
       main=main,
       xlab="number of variables",
       ylab="Error rate",
       type="b",
       pch=1,
       col=col[1],
       panel.first = c(grid(col=1),abline(h=0,col=1),abline(v=0,col=1)),
       lwd=2,
       xlim <- xlim,
       ...
       )
  if (dim(errors)[2] > 1) {
    for (m in 2:dim(errors)[2]) {
      x.values.selected <- x.values[1:length(errors[,m])]
      lines(x.values.selected,errors[,m],ylim=c(0,max(errors,na.rm=T)),type="b",pch=m,col=col[m],lwd=2)
    }
  }
  par(font.lab=1)
  par(font=1)
  if (plot.legend) {
    if (is.na(legend.x)) { legend.x <- max(x.values)*1.02}
    if (is.na(legend.y)) { legend.y <- max(errors,na.rm=T)*0.8 }
    legend(legend.x,
           legend.y,
	   legend.labels,
	   col=col[1:dim(errors)[2]],
	   lwd=2,
#           cex=0.8,
	   bty="n",
	   pch=1:dim(errors)[2]
           )
  }
}

################################################################
## Flexible Discriminant Analysis
################################################################

## ## TRAINING TO TREAT
## ## apparently, fda requires at least 3 families
## tfda <- fda(family ~ .,data=training)
## plot(tfda)
## confusion(tfda)
## fda.pred <- predict(tfda,training)
## tfda.validation <- data.frame(
##                               family=training$family,
##                               predict=fda.pred,
##                               row.names=row.names(training)
##                               )
## confusion(fda.pred,training$family)

################################################################
## Linear discriminant analysis
################################################################

################################################################
## Return misclassified units for a give predictive
## discriminant analysis
##
misclassified.units <- function (data, # data fame containing the multivariate table
				 group.labels, # group labels
				 mpda # predictive discriminant analysis object
                                      ) {
  compa <- data.frame(training.label <- mpda$training.label,
	              class <- mpda$class,
		      row.names=mpda$training.names)
  names(compa) <- c("training.label", "predicted.class")

  #### detect misclassified units
  failures <- compa[compa$predicted.class != compa$training.label,]
  return (failures)
}

################################################################
## Perform predictive linear discriminant analysis,
## and evaluate the result both by internal analysis, and by LOO.
## returns the confusion tables and the profiles of failures.
## Optionally, plots the failure profiles.
##
plot.misclassified.units <- function (data, # data fame containing the multivariate table
                          	      group.labels, # group labels
				      mpda, # predictive discriminant analysis object
				      data.title = "Misclassifications",
				      xlab="variable",
				      ylab="value",
				      display.plot=T # plot the failure profiles
                                      ) {
  failures <- misclassified.units(data,group.labels,mpda)
#    compa <- data.frame(training.label <- mpda$training.label,
#  	              class <- mpda$class,
#  		      row.names=mpda$training.names)
#    names(compa) <- c("training.label", "predicted.class")

#    #### detect misclassified units
#    failures <- compa[compa$predicted.class != compa$training.label,]

  ## Plot profiles of misclassified units
  if (display.plot) {
    par(mai=c(0.7,0.7,0.3,0.1))
    groups <- levels(as.factor(failures$predicted.class))
    par(mfrow=n2mfrow(sum(table(failures$predicted.class) > 0)))
    group.palette <- get.palette(data,group.labels)
    for (g in groups) {
      if(table(failures$predicted.class) > 0) {
        group.failures <- row.names(failures)[!is.na(failures$predicted.class) & failures$predicted.class == g]
        if (length(group.failures > 0)) {
          plot.profiles(data,group.failures,
                        main=paste(data.title, "Classified as",g),
                        colpal=group.palette[g],
                        xlab=xlab,
                        ylab=ylab
                        )
        }
      }
    }
    par(mfrow=c(1,1))
    par(mai=c(1,1,1,1))
  }

  #### return the result
  result <- list()
  result$pda <- mpda
  result$confusion <- mpda$confusion
  result$failures <- failures

  return (result)
}

################################################################
## Plot profiles of objects as a function of their predicted class
##
plot.profiles.by.predicted.class <- function (data,
                                              group.labels,
                                              data.title="",
#                                              method="lda",
#                                              loo = T,
                                              ylim=NA,
                                              ...
                                              ) {

  if (method =="qda") {
    data <- qda.check.variables(data,group.labels)
  }
  mpda <- pda(data,group.labels,...)

  group.palette <- get.palette(data,group.labels)
  groups <- levels(as.factor(group.labels))
  if (is.na(ylim)) {
    ylim <- c(min(data[!is.na(group.labels),],na.rm=T),max(data[!is.na(group.labels),],na.rm=T))
  }
  main <- paste(data.title, method, sep=" - ")
  if (loo) {
    main <- paste(main, "LOO", sep=" - ")
  } else {
    main <- paste(main, "Internal", sep=" - ")
  }

  par(mai=c(0.3,0.5,0.3,0.1))
  par(font.main=2)
  par(cex.main=1.5)
  par(mfrow=n2mfrow(length(groups)))
  for (g in groups) {
    ## detect NA in the predictions (qda sometimes allocates NA)
    problems <-  mpda$training.names[is.na(mpda$class)]
    if (length(problems) > 0) {
      verbose ("discarded units",1)
      verbose (problems)
    }

    g.pred <- mpda$training.names[!is.na(mpda$class) & mpda$class == g]
    plot.profiles(data,g.pred,
                  main=paste(main, "-", "Predicted as", g),
                  colpal=group.palette[g])
  }
  par(mfrow=c(1,1))
  par(mai=c(1,1,1,1))
}


################################################################
#
# Plots the training set, with a color representing the training
# class, and a character representing the predicted class.
#
################################################################

plot.pda.calib <- function (data,
			    group.labels,
                            pda.method="lda",
#                            loo=T,
			    data.title = "",
			    main=NA,
			    mai=c(1.5,1.5,1,1),
			    cex=1.5,
			    cex.main=1,
			    cex.lab=1,
			    font=2,
			    font.main=2,
			    font.lab=2,
                            ...
			    ) {## plot the calibration result

  par(mai=mai)
  par(cex=cex)
  par(cex.main=cex.main)
  par(cex.lab=cex.lab)
  par(font=font)
  par(font.main=font.main)
  par(font.lab=font.lab)

  if (pda.method =="qda") {
    data <- qda.check.variables(data,group.labels)
  }

  mpda <- pda(data,group.labels,...)
  coords <- calc.coordinates(data,group.labels)
  coords.training <- coords[!is.na(group.labels),]
  colors.training <- get.colors.training(data,group.labels)

  if (is.na(main)) {
    main <- paste(data.title, "-",  pda.method)
    if (loo) {
      main <- paste(main, "-",  "LOO")
    } else {
      main <- paste(main, "-",  "Internal")
    }
  }
  plot(coords.training,
       col=colors.training,
       main=main,
       pch=substr(as.vector(mpda$class),1,1),
       xlab=names(coords)[1],
       ylab=names(coords)[2]
       )
  par(mai=c(1,1,1,1))
  par(cex=1)
  par(cex.main=1)
  par(cex.lab=1)
  par(font=1)
  par(font.main=1)
  par(font.lab=1)

}



################################################################
#
# IN CONSTRUCTION
#
plot.lda.prediction <- function (data,
				 group.labels,
				 tlda=NA, #linear discriminant function
                                 priors=NA, # if not provided, calculated automatically
				 data.title="") {
  ## define reasonable prior for whole-genome predictions
  if (is.na(priors)) {
    priors <- append(table(group.labels),length(is.na(group.labels)))
    priors <- priors/sum(priors)
  }
#  names(priors)[length(priors)] <- "NA"

  tlda$prior <- priors
  lda.pred<- data.frame(predict(tlda,data))
  lda.pred.pho <- row.names(lda.pred[lda.pred$class=="PHO",])
  lda.pred.met <- row.names(lda.pred[lda.pred$class=="MET",])

  ## predict in the whole gene set
  ##  if (print.files) {postscript(paste(fig.dir,file.prefix,'_LDA_prediction.ps',sep=""))}
  plot(lda.pred$x.LD1,lda.pred$x.LD2,
       main=paste("Whole data set predictions - ", data.title),
       col="#999999",
       xlab="Linear discriminant function 1",ylab="Linear discriminant function 2"
       )
  if (length(lda.pred.pho) > 0) {
    points(lda.pred[lda.pred.pho,"x.LD1"],lda.pred[lda.pred.pho,"x.LD2"],col=palette.families["PHO"],pch=19)
    text(lda.pred[lda.pred.pho,"x.LD1"],lda.pred[lda.pred.pho,"x.LD2"],group.labels=lda.pred.pho,pos=1,font=2,cex=0.8,col=palette.families["PHO"])
  }
  if (length(lda.pred.met) > 0) {
    points(lda.pred[lda.pred.met,"x.LD1"],lda.pred[lda.pred.met,"x.LD2"],col=palette.families["MET"],pch=19)
    text(lda.pred[lda.pred.met,"x.LD1"],lda.pred[lda.pred.met,"x.LD2"],group.labels=lda.pred.met,pos=1,font=2,cex=0.8,col=palette.families["MET"])
  }

  ##  if (print.files) {dev.off()}

  ## count predictions by class
  dim(lda.pred[lda.pred$class=="PHO",])
  dim(lda.pred[lda.pred$class=="MET",])
  dim(lda.pred[lda.pred$class=="CTL",])
}

################################################################
## Quadratic discriminant analysis
################################################################
my.qda <- function (data,
		    group.labels,
		    data.title="") {
  tqda <- qda(group ~ .,data=training,na.action=na.omit)
  tqda.loo <- qda(group ~ .,training,na.action=na.omit,CV=TRUE)

  qda.compa <- data.frame(
                          group=training$group,
                          loo=tqda.loo$class,
                          class=predict(tqda,training)$class,
                          predict(tqda,training)$posterior,
                          row.names=row.names(training)
                          )

  ## confusion matrices
  confusion(predict(tqda)$class,training$group)
  confusion(tqda.loo$class,training$group)

  failures.loo <- qda.compa[(qda.compa$group != qda.compa$loo),]
  failures.loo
  failures <- qda.compa[(qda.compa$group != qda.compa$class),]
  failures
  training[row.names(failures),]
  ##biplot(tqda$PHO,tqda$MET, main=paste("linear discriminant analysis - ", data.title))

  ## plot the calibration result
  ##  if (print.files) {postscript(paste(fig.dir,file.prefix,'_QDA_internal_analysis_plot.ps',sep=""))}
  pca.training <- prcomp(data[!is.na(group.labels),])
  par(mai=c(0.7,0.7,0.7,0.1))
  par(font=2)
  par(cex=0.7)
  par(cex.main=2)
  par(cex.axis=1)
  plot(pca.training$x[,1:2],
       col=colors.training,
       main=paste("QDA - internal analysis - ", data.title),
       pch=substr(as.vector(qda.compa$class),1,1),
       xlab="PC1",ylab="PC2"
       )
  ##  if (print.files) {dev.off()}

  ##  if (print.files) {postscript(paste(fig.dir,file.prefix,'_QDA_leave_one_out_plot.ps',sep=""))}
  plot(pca.training$x[,1:2],,
       col=colors.training,
       main=paste("QDA - LOO - ", data.title),
       pch=substr(as.vector(qda.compa$loo),1,1),
       xlab="PC1",ylab="PC2"
       )
  par(cex=1)
  par(cex.main=1)
  par(cex.axis=1)
  par(font=1)
  par(mai=c(1,1,1,1))
  ##  if (print.files) {dev.off()}

  ## define reasonable prior for whole-genome predictions
  priors <- c(0.01,0.01,0.98)
  names(priors) <- c("PHO","MET","CTL")
  tqda$prior <- priors

  ## predict in the whole gene set
  qda.pred <- data.frame(predict(tqda,data))
  qda.pred.pho <- row.names(qda.pred[qda.pred$class=="PHO",])
  qda.pred.met <- row.names(qda.pred[qda.pred$class=="MET",])
  qda.pred.ctl <- row.names(qda.pred[qda.pred$class=="CTL",])

  length(qda.pred.pho)
  length(qda.pred.met)

  ## plot the predictions on the PCA plane
  ##  if (print.files) {postscript(paste(fig.dir,file.prefix,'_QDA_prediction.ps',sep=""))}
  pca.data <- prcomp(data)
  par(mai=c(1,1,0.7,0.1))
  plot(pca.data$x[,1:2],
       col=palette["grayB"],
       main=paste("QDA - Whole data set predictions - ", data.title),
       xlab="PC1",ylab="PC2"
       )

  if (length(qda.pred.pho) > 0) {
    points(pca.data$x[qda.pred.pho,1:2],col=palette.families["PHO"],pch=19)
    text(pca.data$x[qda.pred.pho,1:2],group.labels=qda.pred.pho,pos=1,font=2,cex=0.8,col=palette.families["PHO"])
  }
  if (length(qda.pred.met) > 0) {
    points(pca.data$x[qda.pred.met,1:2],col=palette.families["MET"],pch=19)
    text(pca.data$x[qda.pred.met,1:2],group.labels=qda.pred.met,pos=1,font=2,cex=0.8,col=palette.families["MET"])
  }
  par(mai=c(1,1,1,1))
  ##  if (print.files) {dev.off()}

################################################################
## Plot profiles of predicted PHO and MET genes
##  if (print.files) {postscript(paste(fig.dir,file.prefix,'_QDA_pred_profiles.ps',sep=""))}
  par(mai=c(0.7,0.7,0.7,0.1))
  par(mfrow=c(2,2))
  if (length(qda.pred.pho > 0)) {
    plot.profiles(data,qda.pred.pho,
                  main="QDA - Predicted as PHO",
                  colpal=colors.all[qda.pred.pho])
  }
  if (length(qda.pred.met > 0)) {
    plot.profiles(data,qda.pred.met,
                  main="QDA - Predicted as MET",
                  colpal=colors.all[qda.pred.met])
  }
  if (length(qda.pred.ctl > 0)) {
    plot.profiles(data,qda.pred.ctl,
                  main="QDA - Predicted as CTL",
                  colpal=colors.all[qda.pred.ctl])
    ##	colors.training[qda.pred.ctl])
  }
  par(mfrow=c(1,1))
  ##  if (print.files) {dev.off()}

  ## count predictions by class
  length(qda.pred.pho)
  length(qda.pred.met)
  length(qda.pred.ctl)

  par(mai=c(1,1,1,1))
}

## WORKS UNDER WINDOWS BUT NOT UNDER LINUX
##biplot(tlda, main=paste("linear discriminant analysis - ", data.title))


## ###############################################################
## Plot distributions of discriminant function,
## per training group
## usage: plot.lda.function.distrib(pred)
plot.lda.function.distrib <- function (pred,            # prediction object
				       group.labels,    # a vector with the training group labels
			               pooled = F,      # a separate plot with the pooled distribution
				       by.group = T     # plot distributions by group
			               ) {
  groups <- levels(pred$class)
  group.palette <- palette[1:length(groups)]

  xlim <- c(min(pred$x),max(pred$x))
  breaks <- pretty(xlim,n=50)

  if ((pooled) && (by.group)) {
    par(mfrow=c(2,1))
  }

  if (pooled) {
    h <- hist(pred$x,
              breaks=breaks,
	      plot=F)
    plot(h$mids,
         h$counts,
         lwd=3,
         breaks=breaks,
         type="l",
         xlab = "LD",
         ylab = "counts",
         main="linear discriminant function, pooled distribution",
         col="#888888")
  }
  if (by.group) {
    for (g in 1:length(groups)) {
      group <- groups[g]
      h <- hist(pred$x[group.labels==group],
                breaks=breaks,
                plot=F)
      if (g == 1) {
        plot(h$mids,
             h$counts,
	     lwd=3,
	     breaks=breaks,
	     type="l",
	     xlab = "LD",
	     ylab = "counts",
	     main="linear discriminant function, distribution by training group",
	     col=group.palette[g])
      } else {
        lines(h$mids,h$counts,lwd=3,type="l", col=g)
      }
    }
    legend((0.6*xlim[2]-0.4*xlim[1]),
           max(h$counts),
	   legend=groups,
	   col=group.palette,
	   lwd=3)

  }

  if ((pooled) && (by.group)) {
    par(mfrow=c(1,1))
  }
}



## ##############################################################
## Methods developed by Najib Naamane
## ##############################################################

################################################################
## K.fold.iterative.stepwise.pda
## Developed by Najib naamane
kfold.iterative.stepwise.pda <- function (x, ## Data frame witht eh data +  group labels
                                          ngroup = n, ## Number of groups for the K-fold validation
                                          np, ## Number of positives
                                          nn, ## Number of negatives
                                          max.cycles,
                                          method="lda") {

  ##variables
  call <- match.call()
  l <- dim(x)[1]
  w<-dim(x)[2]
  results<-vector()
  conftables<-list()
  validationTables<-list()
  results$rounds<-list()
  results$validationTable<-vector()
  results$m.sets<-list()
  results$t.sets<-list()
  results$folds<-list()
  results$variables<-matrix(data= 0,nrow= w-2,ncol=max.cycles*ngroup)
  results$variables<-as.data.frame(results$variables)
  row.names(results$variables)<-names(x)[2:(w-1)]
  names <- vector()
  oldone = 0
  col = 0

  ##split the ids of the data set into n stratified sets

  if (ngroup < 2) {
    stop("ngroup should be greater than or equal to 2")
  }
  if (ngroup > l) {
    stop("ngroup should be less than or equal to the number of observations")
  }
  if (ngroup == l) {
    groups <- 1:l
    leave.out <- 1
  }
  if (ngroup < l) {

    leave.outp <- trunc(np/ngroup)
    o <- sample(1:np)
    groupsp <- vector("list", ngroup)
    for (j in 1:(ngroup - 1)) {
      jj <- (1 + (j - 1) * leave.outp)
      groupsp[[j]] <- o[jj:(jj + leave.outp - 1)]
    }
    groupsp[[ngroup]] <- o[(1 + (ngroup - 1) * leave.outp):np]


    leave.outn <- trunc(nn/ngroup)
    o <- sample((np+1):l)
    groupsn <- vector("list", ngroup)
    for (j in 1:(ngroup - 1)) {
      jj <- (1 + (j - 1) * leave.outn)
      groupsn[[j]] <- o[jj:(jj + leave.outn - 1)]
   }
    groupsn[[ngroup]] <- o[(1 + (ngroup - 1) * leave.outn):nn]
  }

  ## actions for each fold
  for (j in 1:ngroup) {
    verbose ("fold analysis started")

    ## construct the training set

    a<- x[setdiff(c(1:np),groupsp[[j]]),]
    b<- x[setdiff(c((np+1):l),groupsn[[j]]),]
    t.set<-rbind(a,b)

    ##t.set<-rbind( x[setdiff(c(1:np),groupsp[[j]]),] , x[setdiff(c((np+1):l),groupsn[[j]]),])
    results$t.sets[[j]]<-t.set
    ## construct the test set
    c<- x[groupsp[[j]],]
    d<- x[groupsn[[j]],]
    m.set<-rbind(c,d)
    ## m.set<-rbind( x[groupsp[[j]],] , x[groupsn[[j]],] )
    results$m.sets[[j]]<-m.set

    ## perform an iterative PDA

    u <- najib.iterative.stepwise.pda(t.set,max.cycles = max.cycles,performed.cycles=0,method=method)
    results$folds[[j]]<-u

    ##chek if the number of rounds of the last fold is greater the precedent.

    if ((u$performed.cycles > oldone) && (j!=1)) {
      for(i in (oldone+1):u$performed.cycles) {
        results$rounds[[i]]<-results$rounds[[i-1]]
      }
    }

    ##report evaluation metrics for each round

    for(k in 1:u$performed.cycles){

      ##conftables[[k]]<-prediction(t.set[,u$performed.pdas[[k]]$best.var.names],t.set[,32], m.set[,2:31], m.set[,32],confusion=TRUE,method=method,name= "results.txt")$confusion.table
      varnames<-u$performed.pdas[[k]]$best.var.names
      ma<-as.matrix(u$performed.pdas[[k]]$data[,varnames])
      ve<-u$performed.pdas[[k]]$data[,w]
      if (method == "qda") {
        qda<-qda(ma,ve)
        predict<-predict(qda,as.matrix(m.set[,varnames]))
      }else if(method == "lda") {
        lda<-lda(ma,ve)
        predict<-predict(lda,as.matrix(m.set[,varnames]))
      } else {
        print("ERROR: unknown method")
        return(NULL)
      }

      conftables[[k]]<-confusion(predict$class,m.set[,w])

      ##conftables[[k]]<-prediction(a,b,c,d,confusion=TRUE,method=method,name= "results.txt")$confusion.table

      validationTables[[k]]<-validation(as.integer(conftables[[k]][4]),as.integer(conftables[[k]][1]),as.integer(conftables[[k]][2]),as.integer(conftables[[k]][3]))

      if(j==1){
        results$rounds[[k]]<-vector()
      }
      results$rounds[[k]]<-cbind(results$rounds[[k]],validationTables[[k]][,1])
      ##construct variables table
      names <- cbind(names, paste("R",k,"F",j,sep="."))
      col = col+1
      results$variables[intersect(row.names(results$variables),varnames),col]<-1
    }

    ##chek if the number of rounds of the last fold is smaller than the precedent.
    if ((u$performed.cycles < oldone) && (j!=1)) {
      for (i in (u$performed.cycles+1):length(results$rounds)){
        results$rounds[[i]]<-cbind(results$rounds[[i]],validationTables[[k]][,1])
      }
    }
    oldone<-length(results$rounds)
  }

  ##calculate occurences of the variables
  names(results$variables)<-names
  occur<-vector()
  for(i in 1:(w-2)){
    occur<-append(occur,sum(results$variables[i,]))
  }
  results$variables<-cbind(results$variables,occur)

  ##calculate the means for the n fold results
  for(k in 1:length(results$rounds)){
    sum<-vector()
    mean<-vector()
    sd<-vector()
    for(i in 1:13){
      sum<-append(sum,sum(results$rounds[[k]][i,]))
      mean<-append(mean,mean(results$rounds[[k]][i,]))
      sd<-append(sd,sd(results$rounds[[k]][i,]))
    }
    results$rounds[[k]]<-cbind(results$rounds[[k]],sum,mean,sd)
  }

  ##plots
  results$toPlot<-vector()
  lr<-length(results$rounds)
  for(i in 1:lr){
    results$toPlot<-cbind(results$toPlot,results$rounds[[i]][c("Sn", "Sp","PPV","Acc.g"),(ngroup+2)])
  }
  results$toPlot<-round(results$toPlot,digits=2)
  barplot(t(results$toPlot),
          beside=TRUE,
          xlab="Evaluation metrics",
          ylab="Percentage",
          main=paste("Average performance of",
            method,"with different preprocessed data sets"),
          legend.text= c(paste("round",1:lr)),
          xlim = c(0,(1.5*lr*3)+5),ylim = c(0,1))

                                        #toPlot<-vector()
                                        #for(i in 1:length(x5F10R1$rounds)){
                                        #  toPlot<-cbind(toPlot, x5F10R1$rounds[[i]][5:7,7])
                                        #}
                                        #barplot(t(toPlot),beside=TRUE)

  return(results)
}


## #########################################################
## Iterative predictive linear analysis
## First cycle:
## - perform a predictive linear analysis
## At each new cycle :
## - detect the misclassified units from the previous cycle
## - remove these misclassified units from the training set
## - perform a predictive linear analysis
############################################################

najib.iterative.stepwise.pda <- function (dataset,
                                       max.cycles=2, # maximum number of cycles for the iteration
                                       performed.cycles=0, # count of the cycles already performed
                                       export.formats.plots=NA, # export results at each cycle
                                       export.formats.obj=NA, # export results at each cycle
                                       data.type = "data", # prefix for saving files
                                       method="lda",
                                       performed.pdas= NA,
                                       ...
                                       ) {


  ## verbose
  verbose ("Iterative stepwise PDA")
  verbose (paste("Performed cycles : ", performed.cycles, "; Remaining cycles : ", max.cycles))

                                        # check number of cycles
  if (max.cycles < 1) {
    stop ("cycles should be >=1")
  }
                                        # chek
  if (performed.cycles==0){
    performed.pdas<-vector("list")
  }

  w<-dim(dataset)[2]
  data<-dataset[,2:(w-1)]
  group.labels<-dataset[,w]

  s <- stepwise.pda (data,group.labels, pda.method= method, ...)


  performed.cycles <- performed.cycles + 1
  s$miscl <- misclassified.units(data,group.labels,s$best.pda)
  temp<- group.labels
  temp[as.numeric(row.names(s$miscl))] <- NA
  s$group.labels.filtered<-temp[which(temp[]!="NA")]
  s$max.cycles <- max.cycles
  s$performed.cycles <- performed.cycles
  s$group.labels<- group.labels
  s$data<-dataset
  s$filtered.data<-dataset[which(temp[]!="NA"),]
  performed.pdas[[performed.cycles]]<-s
  s$performed.pdas<-performed.pdas



                                        #file.name <- paste(data.type,"iterative_sda_cycle",performed.cycles, sep="_")
                                        #export.object(s, file=file.name, export.formats=export.formats.obj)
                                        #export.plot(file.prefix=file.name, export.formats=export.formats.plots)

  if (length(s$group.labels.filtered)==length(group.labels)) {
    message.verbose("iteration converged because group labels are not modified anymore", 0)
    return (s)
  } else if (max.cycles == 1) {
    return (s)
  } else {
    ## next cycle
    return (najib.iterative.stepwise.pda(s$filtered.data,
                                      max.cycles=max.cycles-1,
                                      performed.cycles=performed.cycles,
                                      export.formats.plots=export.formats.plots,
                                      export.formats.obj=export.formats.obj,
                                      method= method,
                                      performed.pdas=performed.pdas,
                                      ...
                                      ))
  }

}


#############
## validation
validation <- function(TP,TN,FP,FN) {

  Sn <- TP/(TP+FN) ## Sensitivity
  Sp <- TN/(TN+FP) ## Specificity (in the sense of rejection power)
  PPV <- TP/(TP+FP) # Positive Predictive Value
  NPV <- TN/(TN+FN) # Negative Predictive Value
  Acc <- (PPV+Sn)/2 ## Accuracy, as arithmetic mean
  Acc.g <- sqrt(PPV*Sn) ## Accuracy, as geometric mean
  Approximate.Coefficient <- 0.5*( TP/(TP+FN) + TP/(TP+FP) + TN/(TN+FP)
                                  + TN/(TN+FN)) - 1## Accuracy following Guigo
  Hit.rate <- (TP+TN)/(TP+TN+FP+FN)
  Error.rate <- (FP+FN)/(TP+TN+FP+FN)

  result <- t(data.frame(
                         TP = round(TP,digits=0), ## True positive
                         TN = round(TN,digits=0), ## True negative
                         FN = round(FN,digits=0), # False negative
                         FP = round(FP,digits=0), # False positive
                         Sn=round(Sn,digits=2),
                         Sp=round(Sp,digits=2),
                         PPV=round(PPV,digits=2),
                         NPV=round(NPV,digits=2),
                         Acc=round(Acc,digits=2),
                         Acc.g=round(Acc.g,digits=2),
                         Approximate.Coefficient = round(Approximate.Coefficient,digits=2),
                         Hit.rate = round(Hit.rate,digits=3),
                         Error.rate = round(Error.rate, digits=3)
                         ))


  return(result)

}

