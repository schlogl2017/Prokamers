################################################################
## R library for exploring RNA-seq data
##
## Author : Jacques van Helden
## First implementation: 2013-06-28
##
## Call:
##   source(file.path(dir.util, "rna-seq_lib.R"))


## Some custom utilities
source(file.path(dir.util, "util.R"))
source(file.path(dir.util, "util_plots.R"))
source(file.path(dir.util, "multitesting_corrections.R"))



################################################################
## Define colors for plots
plot.colors <- c(gene='#BBBBBB',
                 grid='#DDDDDD',
                 test='#0000BB',
                 ctrl='#BB0044',
                 up='red',
                 down='darkgreen',
                 absent='#BB00BB',
                 none='#BBBBBB'
                 )

################################################################
## Use the hypergeometric test (equivalent to Fisher' exact test) to
## detect differentially expressed genes (DEG) between two RNA-seq
## samples. This method applies only for experiments without
## replicate.
deg.hypergeo  <- function (test.counts,
                           ctrl.counts,
                           test.name="test",
                           ctrl.name="ctrl",
                           no0.value = 0.1, ## Replace zero values by a small non-zero values for some graphical representations. 
                           draw.plots=TRUE, ## Draw plots if TRUE
                           threshold.eval=0.01
                           ) {
  
  ## Prepare a data frame to store the result
  test.vs.ctrl <- data.frame(test=test.counts,
                             ctrl=ctrl.counts)
  test.vs.ctrl$sum.per.gene <- test.vs.ctrl$test + test.vs.ctrl$ctrl

  ##  names(test.vs.ctrl) <- c(test.name, ctrl.name)

  counts.per.sample <- apply(test.vs.ctrl, 2, sum)

  ## Replace NA values by zeros
  test.vs.ctrl[is.na(test.vs.ctrl)] <- 0

  ## Create separate columns with zero values 
  test.vs.ctrl$test.no0 <- test.vs.ctrl$test
  test.vs.ctrl$test.no0[test.vs.ctrl$test.no0 == 0] <- no0.value
  test.vs.ctrl$ctrl.no0 <- test.vs.ctrl$ctrl
  test.vs.ctrl$ctrl.no0[test.vs.ctrl$ctrl.no0 == 0] <- no0.value
  test.vs.ctrl$test.log2 <- log(test.vs.ctrl$test.no0,base=2)
  test.vs.ctrl$ctrl.log2 <- log(test.vs.ctrl$ctrl.no0,base=2)


  ## Compute log ratios
  absent.gene <- (test.vs.ctrl$test==0) & (test.vs.ctrl$ctrl ==0)
  test.vs.ctrl$log.ratio <- log(test.vs.ctrl$test / test.vs.ctrl$ctrl, base=2)
  test.vs.ctrl[absent.gene, "log.ratio"] = 1
  test.vs.ctrl$log.ratio.no0 <- log(test.vs.ctrl$test.no0 / test.vs.ctrl$ctrl.no0, base=2)

  ## Compute M and A
  test.vs.ctrl$M=test.vs.ctrl$test.log2 - test.vs.ctrl$ctrl.log2;
  test.vs.ctrl$A=(test.vs.ctrl$test.log2 + test.vs.ctrl$ctrl.log2)/2


  ## ##############################################################
  ## Select differentially expressed gene with a very straightforward
  ## approach: a Fisher test (showe p-value is computed with the
  ## hypergeometric distribution).
  ##
  ## Beware: this method is usually not applied for detecting RNA-seq,
  ## but it applies here because we have no replicates. In this
  ## particular case, the Fisher test and the hypergeometric
  ## distribution thus seems the most appropriate.
  ##
  ## Model:
  ##
  ##  m = number of "marked" elements, i.e. all the read counts in the
  ##      test sample
  ##
  ##  n = number of "non-marked" elements, i.e. all the read counts in
  ##      the control sample
  ##
  ##  k = "selected" elements, i.e. all the counts related to the gene
  ##       of interest (test + control conditions)
  ##
  ##  q = number of marked elements in the selection, i.e. the counts of
  ##      the gene of interest in the test sample.
  i <- 2010
  n.test <- test.vs.ctrl[i,"test"]
  n.ctrl <- test.vs.ctrl[i,"ctrl"]
  total.test <- sum(test.vs.ctrl[,"test"])
  
  ##hyper.pval.up <-phyper(q=, m=test.vs.ctrl$sum, n=
  hyper.pval.up <- phyper(q=test.vs.ctrl$test-1, m=test.vs.ctrl$test, n=test.vs.ctrl$ctrl, k=test.vs.ctrl$test + test.vs.ctrl$ctrl, lower=FALSE)
  hyper.pval.down <- phyper(q=test.vs.ctrl$ctrl-1, m=test.vs.ctrl$ctrl, n=test.vs.ctrl$test, k=test.vs.ctrl$test + test.vs.ctrl$ctrl, lower=FALSE)
  test.vs.ctrl["hyper.pval"] <- hyper.pval.up
  test.vs.ctrl[test.vs.ctrl$test < test.vs.ctrl$ctrl, "hyper.pval"] <- hyper.pval.down[test.vs.ctrl$test < test.vs.ctrl$ctrl]
  
  test.vs.ctrl$hyper.eval <- test.vs.ctrl$hyper.pval * nrow(test.vs.ctrl)
  test.vs.ctrl$hyper.sig <- -log(test.vs.ctrl$hyper.eval, base=10)

  
  ## ## Apply multiple testing corrections
  ## mtc <- multitest.corrections(pval=test.vs.ctrl$hyper.pval, 
  ##                              names=rownames(test.vs.ctrl),
  ##                              plots=TRUE, ## If true, generate illutrsative plots
  ##                              lambda=0.5, ## lambda parameter for estimating pi0 = m0/m1
  ##                              threshold.eval=0.01, ## Threshold.Eval of significance for the plots
  ##                              main=paste(test.name, "vs", ctrl.name), ## Prefix for the main title of the plots
  ##                              file.prefix=file.path(dir.figures, paste(sep="", test.name, "_vs_", ctrl.name, "hypergeo_test")))
  ## test.vs.ctrl <- cbind(test.vs.ctrl, mtc$multitest.table[c("rank","eval", "fwer","qval.Storey", "qval.0")])
  test.vs.ctrl$regul <- "none"
  test.vs.ctrl[absent.gene,"regul"] <- "absent"
  test.vs.ctrl[(test.vs.ctrl$hyper.eval <= threshold.eval) & (test.vs.ctrl$test > test.vs.ctrl$ctrl),"regul"] <- "up"
  test.vs.ctrl[(test.vs.ctrl$hyper.eval <= threshold.eval) & (test.vs.ctrl$test < test.vs.ctrl$ctrl),"regul"] <- "down"
  table(test.vs.ctrl$regul)
  test.vs.ctrl$gene.color <- plot.colors[test.vs.ctrl$regul]

  ## ##############################################################
  ## Plot histogram of log counts (zero values are ignored)
  if (draw.plots) {
    x11(width=7,height=8)
    par(mfrow=c(2,1))
    hist(test.vs.ctrl$test.log2, breaks=100, main=paste(sep="","log2(",test.name,")"), col=plot.colors["test"],
         xlab="log(counts) [0 changed to 0.1]", ylab="Number of genes")
    hist(test.vs.ctrl$ctrl.log2, breaks=100, main=paste(sep="","log2(",ctrl.name,")"), col=plot.colors["ctrl"],
         xlab="log(counts) [0 changed to 0.1]", ylab="Number of genes")
    par(mfrow=c(1,1))
  setwd(dir.figures); export.plot('log2_count_distrib_histo', width=7, height=8, export.formats=c('pdf'))
    
    ## Draw an XY-plot to compare the two conditions.
    X11(width=7,heigh=7)
    plot(test.vs.ctrl$ctrl.log2,
         test.vs.ctrl$test.log2,
         col=test.vs.ctrl$gene.color,
         main=paste(test.name, "vs", ctrl.name),
         panel.first=c(
           grid(),
       abline(h=0, col="black"),
           abline(v=0, col="black")),
         xlab=paste(sep="", "log2(",test.name,")"),
         ylab=paste(sep="", "log2(",ctrl.name,")")
         )
    abline(a=0,b=1, col='black')
    setwd(dir.figures); export.plot('test_vs_ctrl_log2', width=7, height=7, export.formats=c('pdf'))
    
    ## Draw an MA-plot to compare the two conditions. Weakness: we loose
    ## all the zero values, which probably include significant DEG.
    X11(width=7,heigh=7)
    plot(test.vs.ctrl$A,
         test.vs.ctrl$M,
     col=test.vs.ctrl$gene.color,
     main=paste(test.name, "vs", ctrl.name),
     panel.first=c(
       grid(),
       abline(h=0, col="black"),
       abline(v=0, col="black")),
     xlab=paste(sep="", "A=(log2(",test.name,") + log2(",ctrl.name,"))/2"),
     ylab=paste(sep="", "M=log2(",test.name,"/",ctrl.name,")"),
     )
    abline(h=0, col='black')
    setwd(dir.figures); export.plot('test_vs_ctrl_MA', width=7, height=7, export.formats=c('pdf'))
  }
  return(test.vs.ctrl)
}

