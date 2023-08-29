################################################################
#
# Utilities for multivariate analysis
#
################################################################

#### load the multivariate analysis libraries
#library(mda)
library(stats)
#library(mda)
library(MASS)

#### load utilities
source(paste(dir.util,'util.R', sep='/'))
source(paste(dir.util,'util_descr_stats.R', sep='/'))
source(paste(dir.util,'util_chip_analysis.R', sep='/'))
source(paste(dir.util,'util_sda.R', sep='/'))

################################################################
#
# Perform all analyses and export the plots to files in different formats.
#
################################################################

export.all.plots <- function(data,
                             group.labels,
                             export.dir=".",
                             data.title="",
                             file.prefix="",
			     export.formats.plots=c("postscript","jpg"),
			     export.formats.obj=c("print"),
                             ## for stepwise PDA
                             max.p=NA,
                             lda=T,
                             qda=T
                             ) {

  ## check the data
  data <- check.data(data)
  file.prefix.ori <- file.prefix
  data.title.ori <- data.title

  ## determine the max nnumber of variables
  if (is.na(max.p)) {
    max.p <- dim(data)[2]
  }

  ## perform all analyses with the data
  for (princomp in c(F,T)) {
    export.mva.plots(data,group.labels,export.dir=dir.figures,
                     file.prefix=file.prefix,
                     data.title=data.title,
                     princomp=princomp,
                     export.formats.plot=export.formats.plots,
                     export.formats.obj=export.formats.obj
		     )

    export.sda.plots(data,group.labels,export.dir=dir.figures,
                     max.p=max.p,
                     lda=lda,
                     qda=qda,
                     file.prefix=file.prefix,
                     data.title=data.title,
                     princomp=princomp,
                     export.formats.plot=export.formats.plots,
                     export.formats.obj=export.formats.obj
		     )

  }

  ## perform all analyses with randomized data
  rand <- randomize(data)
  for (princomp in c(F,T)) {
    export.mva.plots(rand,group.labels,export.dir=dir.figures,
                     file.prefix=paste(file.prefix,"rand",sep="_"),
                     data.title=paste("randomized", data.title),
                     princomp=princomp,
                     export.formats.plot=export.formats.plots,
                     export.formats.obj=export.formats.obj
		     )

      export.sda.plots(rand,group.labels,export.dir=dir.figures,
                     max.p=max.p,
                     lda=lda,
                     qda=qda,
                     file.prefix=paste(file.prefix,"rand",sep="_"),
                     data.title=paste("randomized", data.title),
                     princomp=princomp,
                     export.formats.plot=export.formats.plots,
                     export.formats.obj=export.formats.obj
		     )
  }

  file.prefix <- file.prefix.ori
  data.title <- data.title.ori

  ## comparison of all the stepwise discrimnant results
#  pda.eval <- compare.stepwise.pda.methods(data,group.labels,max.p=max.p,data.title=data.title)
#  file.name <- paste(dir.figures, paste(file.prefix, "stepwise_PDA_error_profiles", sep="_"),sep="/")
#  export.plot(file.prefix=file.name, export.formats=export.formats.plots)
#  export.object(pda.eval, file=file.name,export.formats=export.formats.obj)

}

export.mva.plots <- function(data,
                     	     group.labels,
                             princomp=F, # perform a principal component transformation before analysis
			     export.dir=".",
			     file.prefix,
			     data.title,
			     export.formats.plots = c("postscript","jpg"),
			     export.formats.obj = c("dput")
			     ) {

  wd <- getwd()
  setwd(export.dir)


  ## check the data
  data <- check.data(data)

  if (princomp) {
    ## principal components with prcomp
    pc <- prcomp.with.group.labels(data,group.labels,data.title=data.title, display.plot=T)
    file.name <- paste(file.prefix,'PC_prcomp',sep='_')
    export.plot(file=file.name,export.formats=export.formats.plots)
    export.object(pc, file=file.name,export.formats=export.formats.obj)

    ## principal components with princomp
    pc2 <- princomp.with.group.labels(data,group.labels,data.title=data.title,display.plot=T)
    file.name <- paste(file.prefix,'PC_princomp',sep='_')
    export.plot(file=file.name,export.formats=export.formats.plots)
    export.object(pc2, file=file.name,export.formats=export.formats.obj)

    data <- pc2$scores
    file.prefix <- paste(file.prefix,"PC",sep="_")
    data.title <- paste(data.title," (PC)")
  }


  ## variable distributions
  plot.frequency.distrib (data,main=paste('Value distributions -', data.title=data.title))
  file.name <- paste(file.prefix,'variable_distributions',sep='_')
  export.plot(file=file.name, export.formats=export.formats.plots)

  #### var pair
  i <-  1
  j <-  2
  plot.var.pair(data,group.labels,data.title=data.title,x=i,y=j)
  file.name <- paste(file.prefix,'vars',i,j,sep='_')
  export.plot(file=file.name,export.formats=export.formats.plots)

  #### group profiles
  plot.profiles.by.class(data,group.labels,data.title=data.title)
  file.name <- paste(file.prefix,'group_profiles',sep='_')
  export.plot(file=file.name,export.formats=export.formats.plots,width=8,height=10,horizontal=F)

  setwd(wd)
}





################################################################
#
# Checks the columns of a data frame, to filter out collinear
# variables.
# Example:
#   x <- matrix(rnorm(100), ncol=5)
#   x <- x[,c(1:5,2,3,1)]
#   y <- check.data(x)
#
################################################################
check.data <- function (data,
		        cor.threshold=0.9999 ### threshold on correlation for removing collinear variables
			) {
#  data <- as.data.frame(data)
  p <- dim(data)[2]
  vars.to.delete <- vector()
  c <- cor(data,use="complete.obs")
  message.verbose(c,4)
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      message.verbose(paste("checking vars", i, j, c[i,j], sep="	"),4)
      if (c[i,j] >= cor.threshold) {
         message.verbose(paste("collinear variables", i,  j,  "deleting", j,sep="	"),1)
	 vars.to.delete <- union(vars.to.delete,j)
      }
    }
  }
  message.verbose(vars.to.delete,3)
  data <- data[,setdiff(1:p,vars.to.delete)]
  return(data)
}

## ##############################################################
## Discard varibles which are constant within groups
discard.group.constants <- function (data,
                                     group.labels,
                                     ...
                                     ) {
  groups <- unique(na.omit(group.labels))

  ## Discard variables which are constant constant within some group
  vars.to.delete <- vector()
  for (g in groups) {
    colvar <- apply(data[!is.na(group.labels) & group.labels==g,], 2, var)
    vars.to.delete <- union(vars.to.delete,names(data)[colvar == 0])
  }
  data <- data[,setdiff(names(data),vars.to.delete)]

  return (data)
}

################################################################
#
# Assign coordinates to objects according
# either to linear discriminants, or to PCA.
#
calc.coordinates <- function(data,group.labels) {
  ## calculate point coordinates
  groups <- levels(as.factor(group.labels))
  if (length(groups) > 2) {
    ## use LDA to assign coordinates to the points (the two first LDs)
    tlda <- lda(data,group.labels,na.action=na.omit)
    coords <- as.data.frame(predict(tlda,data)$x)
  } else {
    ## use PCA on the labelled objects only
    pc <- prcomp.with.group.labels(data,group.labels,display.plot=F)
    coords <- as.data.frame(pc$x[,1:2])
  }
  return (coords)
}



#### obsolete
#colors.training <- c(
#                     rep(palette.families["PHO"],dim(pho.fam.profiles)[1]),
#                     rep(palette.families["MET"],dim(met.fam.profiles)[1]),
#                     rep(palette.families["CTL"],dim(ctl.fam.profiles)[1])
#                     )
#names(colors.training) <- row.names(training)
#
#colors.all <- rep(palette["grayB"],length.out=dim(data)[1])
#names(colors.all) <- row.names(data)
#colors.all[names(colors.training)] <- colors.training

################################################################
# Randomize a data frame.
# The program returns a data frame of random numbers,
# with the same characteristics as the input data frame :
# - dimensions )n,p)
# - covariance matrix
# - row and column names
################################################################
randomize <- function (data, # a data frame
                       theor.distrib="mvrnorm" ## currently, the only supported;
                       #### I should add "normal" and "binomial"
		       ) {
  #### calculate parameters from input data
  n<- dim(data)[1]
  p<- dim(data)[2]
#  means <- vector(length=p)
#  for (i in 1:p) {
#    means[i] <- mean(data[,i])
#  }
  means <- apply(data,2,mean)
  Sigma <- cov(data)

  #### generate the multivariate random numbers
  rdata <- data.frame(mvrnorm(n,means,Sigma))
  names(rdata) <- names(data)
  row.names(rdata) <- row.names(data)

  #### return the result
  return (rdata)
}


## ##############################################################
## Perform one separate randomization per group
## For each group, random numbers are generated using the mvrnorm()
## function with the covariance matrix of that group.
randomize.by.groups <- function (data, # a data frame
		                 group.labels,
				 theor.distrib="mvrnorm" ## currently, the only supported;
				 #### I should add "normal" and "binomial"
				 ) {
  groups <- levels(as.factor(group.labels))

  ## randomize labeled objects
  for (g in groups) {
    data.group <- data[group.labels==g,]
    n <- dim(data.group)[1]
    means <- apply(data.group,2,mean)
    Sigma <- cov(data.group)
    data[group.labels==g,] <- data.frame(mvrnorm(n,means,Sigma))
  }

  ## randomize non-labeled objects
  data.group <- data[is.na(group.labels),]
  n <- dim(data.group)[1]
  means <- apply(data.group,2,mean)
  Sigma <- cov(data.group)
  data[is.na(group.labels),] <- data.frame(mvrnorm(n,means,Sigma))


  return (data)
}

################################################################
# Plot a pair of variables
################################################################
plot.var.pair <- function (	data, # a data frame
				group.labels, # group.labels associated to objects
				data.title="",
				x=1, # first variable
				y=2 # second variable
			) {
  groups <- levels(as.factor(group.labels))
  group.palette <- get.palette(data,group.labels)
#  print (group.palette)
  plot(data[,x],data[,y],
       col=palette["grayB"],
       pch=20,
       main=paste(data.title,'- selected variables', x, 'and' , y),
       xlab=names(data)[x],
       ylab=names(data)[y]
#       xlab=paste("variable", x],
#       ylab=paste("variable", y)
       )

  for (g in groups) {
    points(data[group.labels==g,x],data[group.labels==g,y],col=group.palette[g],pch=19)
  }
  legend(min(data[,x],na.rm=T),max(data[,y],na.rm=T),
         names(c(palette["grayB"],group.palette)),
         c(palette["grayB"],group.palette),
         lwd=3,
         cex=0.8)
}

##  ################################################################
##  # Draw profiles
##  ################################################################
##
##  plot.group.profiles <- function (data,group.labels) {
##  	par(mfrow=(n2mfrow(length(groups))))
##          group.palette <- get.palette(data,group.labels)
##  	for (g in groups) {
##  		plot.profiles(data[group.labels==g,],main=paste(data.title, "-", g),colpal=group.palette[g])
##  	}
##  	par(mfrow=c(1,1))
##  }

################################################################
# Plot profiles of objects as a function of their predefined class
#
plot.profiles.by.class <- function (data,
                                    group.labels,
				    data.title = "",
				    plot.legend = T, # plot the legend
                                    ylim=NA
                                    ) {
  groups <- levels(as.factor(group.labels))
  group.palette <- get.palette(data,group.labels)
  par(mfrow=n2mfrow(length(groups)))
  if (is.na(ylim)) {
    ylim <- c(min(data[!is.na(group.labels),],na.rm=T),max(data[!is.na(group.labels),],na.rm=T))
  }
  for (g in groups) {
    plot.profiles(data[group.labels==g,],
                  main=paste(data.title, "- group", g),
                  ylim=ylim,
		  plot.legend=plot.legend,
                  colpal=group.palette[g])
  }
  par(mfrow=c(1,1))
}


################################################################
#### test some basic multivariate functions
################################################################

# #### calculate the variance, covariance and correlation matrices between chips
# var(data,use="complete.obs")
# cov(data,use="complete.obs")
# cor(data,use="complete.obs")

# #### eigen values
# eigen(var(data,use="complete.obs"))

################################################################
#### principal component analysis with prcomp
################################################################
prcomp.with.group.labels <- function (data,
                                      group.labels,
				      leg.pos=NA,
				      data.title="",
                                      whole.data=F,
                                      display.plot=T) {
  if (whole.data) {
    print ("Calculating principal components on the whole data set")
    pc <- prcomp(data)
    plot.main <- paste("PCA with whole data set - ",data.title, sep="")
  } else {
    ## calculate the PC on the labeled objects only
    verbose ("Calculating principal components on the labeled objects only",2)
    pc <- prcomp(data[!is.na(group.labels),])
    plot.main <- paste("PCA with labelled objects - ",data.title, sep="")
    ## rotate the whole data set
    verbose ("Rotating the whole data set",2)
    data.scaled <- scale(data,center=T,scale=T)
    pc$x <- as.matrix(data.scaled) %*% pc$rotation
  }


  if (display.plot) {
    pc1 <- pc$x[,"PC1"]
    pc2 <- pc$x[,"PC2"]
    plot(pc1,pc2,
         main        = paste("PCA - ", data.title,sep=""),
         col         = palette["grayB"],
         xlab        = "PC1",
         ylab        = "PC2",
         panel.first = c(grid(col=1),
           abline(h=0,col=1),
           abline(v=0,col=1)),
         pch=20)
    group.palette <- get.palette(data,group.labels)
    groups <-  levels(as.factor(group.labels))
    for (g in groups) {
      points(pc1[group.labels==g],pc2[group.labels==g],
             col         = group.palette[g],
             pch=19)
    }
    if (is.na(leg.pos)) {
      leg.pos <- c(min(pc1),max(pc2))
    }
    legend(leg.pos[1],leg.pos[2],
           c("all",names(group.palette)),
           c(palette["grayB"],group.palette),
           cex=0.8)
  }
  return(pc)
}


################################################################
#### principal component analysis with princomp (eigen values)
################################################################
princomp.with.group.labels <- function (data,
                                        group.labels,
					data.title="",
                                        whole.data=F,
                                        display.plot=T) {

  ## calculate the principal components on the basis of the labeled
  ## objects only, and rotates the whole data set on the basis of these
  ## loadings
  ## princomp.with.group.labels(data)
  ## princomp.with.group.labels(data,group.labels)

  if (whole.data) {
    print ("Calculating principal components on the whole data set, using correlation matrix")
    pc <- princomp(data,cor=TRUE)
    plot.main <- paste("PCA with whole data set - ",data.title, sep="")
  } else {
    ## calculate the PC on the labeled objects only
    print ("Calculating principal components on the labeled objects only, using correlation matrix")
    pc <- princomp(data[!is.na(group.labels),],cor=TRUE)
    plot.main <- paste("PCA with labelled objects - ",data.title, sep="")

    print ("Rotating and scaling the whole data set")
    data.scaled <- scale(data,center=T,scale=T)
    pc$scores <- as.matrix(data.scaled) %*% pc$loadings
  }

  if (display.plot) {
    par(mfrow=c(1,2))
    par(cex=0.8)

    ## plot component variances
    screeplot(pc, npcs=length(pc$sdev),main=paste("PCA component variances - ",data.title,sep=""))

    ## plot objects along the two principal components
    biplot(pc,
           main=plot.main,
           col=c("#0000bb","#ff0000"))
    par(cex=1)
    par(mfrow=c(1,1))
  }
  return (pc)
}

################################################################
#### hierarchical clustering
################################################################

plot.my.hclust <- function (data,
			    group.labels,
			    data.title="") {
  par(mfrow=c(2,2))
## cluster the columns (variables) on the whole data set
  hc.var <- hclust(dist(t(data)), "ave")
  plot(hc.var,
       main=paste('Hierarchical clustering of variables - ',
         data.title),
       hang=-1)

## cluster the columns (variables) on the training set
  training.hc.var <- hclust(dist(t(data[!is.na(group.labels),])), "ave")
  plot(training.hc.var,
       main=paste('Hierarchical clustering of variables - ',
         data.title),
       hang=-1)


## cluster the rows (objects)
  hc.obj <- hclust(dist(data[!is.na(group.labels),]), "ave")
  par(cex=0.5)
  par(cex.main=2)
  plot(hc.obj,
       main=paste('Hierarchical clustering of objects - ',
         data.title),
       hang=-1,
       col="#008800")
  par(cex=1)
  par(mfrow=c(1,1))
}

