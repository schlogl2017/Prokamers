################################################################
##
## Starting from a contingency table, calculate
## various contingency statistics
##
## See Brohee & van Helden (2006) for an example of utilization:
## comparison between annotated protein complexes and clusters
## obtained with various graph clustering methods.
##
## The input matrix is not mandatorily square.
##
## source(file.path(dir.util, "contingency_stats.R"))

contingency.stats <- function (x, ## A matrix or data frame containing a contingency table (counts)
                               row.label = "row",
                               col.label = "col",
                               row.margins = apply(x,1,sum),
                               col.margins = apply(x,2,sum)
                               ) {

  ## ##############################################################
  ## Total counts per row and per column
#  row.margins <- apply(x,1,sum)
#  col.margins <- apply(x,2,sum)

  row.proportions <- x/row.margins
  row.proportions[is.na(row.proportions)] <- 0
  col.proportions <- t(t(x)/col.margins)
  col.proportions[is.na(col.proportions)] <- 0

  ## Positive predictive value per column
  col.PPV <- apply(x, 2, max)/col.margins
  col.PPV[is.na(col.PPV)] <- 0
  col.PPV.avg <- mean(col.PPV)
  col.PPV.avg.weighted <- sum(col.PPV*col.margins)/sum(col.margins)

  ## Sensitivity per row
  row.Sn <- apply(x, 1, max)/row.margins
  row.Sn[is.na(row.Sn)] <- 0
  row.Sn.avg <- mean(row.Sn)
  row.Sn.avg.weighted <- sum(row.Sn*row.margins)/sum(row.margins)

  ## Accuracy (arithmetic mean of PPV and Sn), weighted or not
  acc.arith <- (col.PPV.avg+row.Sn.avg)/2
  acc.arith.weighted <- (col.PPV.avg.weighted+row.Sn.avg.weighted)/2

  ## Geometric accuracy (geometric mean of PPV and Sn), weighted or not
  acc.geom <- sqrt(col.PPV.avg*row.Sn.avg)
  acc.geom.weighted <- sqrt(col.PPV.avg.weighted*row.Sn.avg.weighted)
  
  ## ##############################################################
  ## Cell accuracy (=sqrt(Compactness))
  cell.acc.geom <- sqrt(row.proportions * col.proportions)
  cell.acc.geom.rowmax <- apply (cell.acc.geom, 1, max)
  cell.acc.geom.colmax <- apply (cell.acc.geom, 2, max)
  cell.acc.geom.rowmax.colavg <- mean(cell.acc.geom.rowmax)
  cell.acc.geom.colmax.rowavg <- mean(cell.acc.geom.colmax)
  cell.acc.geom.geom.rowcol <- sqrt(cell.acc.geom.rowmax.colavg*cell.acc.geom.colmax.rowavg)
    
#   row.cell.acc.geom <- apply(cell.acc.geom, 1, sum)
#   row.cell.acc.geom.avg <- mean(row.cell.acc.geom)
#   row.cell.acc.geom.avg.weighted <- sum(row.cell.acc.geom*row.margins)/sum(row.margins)
  
#   col.cell.acc.geom <- apply(cell.acc.geom, 2, sum)
#   col.cell.acc.geom.avg <- mean(col.cell.acc.geom)
#   col.cell.acc.geom.avg.weighted <- sum(col.cell.acc.geom*col.margins)/sum(col.margins)

  ## ##############################################################
  ## Compactness
  compactness <- row.proportions * col.proportions
  
  row.compactness <- apply(compactness, 1, sum)
  row.compactness.avg <- mean(row.compactness)
  row.compactness.avg.weighted <- sum(row.compactness*row.margins)/sum(row.margins)
  
  col.compactness <- apply(compactness, 2, sum)
  col.compactness.avg <- mean(col.compactness)
  col.compactness.avg.weighted <- sum(col.compactness*col.margins)/sum(col.margins)

  geom.compactness <- sqrt(col.compactness.avg*row.compactness.avg)
  geom.compactness.weighted <- sqrt(col.compactness.avg.weighted*row.compactness.avg.weighted)

  ## ##############################################################
  ## Fill the result object
  result <- list(
                 x = x,
                 row.proportions=row.proportions,
                 col.proportions=col.proportions,
                 compactness=compactness,
                 cell.acc.geom=cell.acc.geom,

                 row.stats=data.frame(
                   row.margins,
                   row.Sn,
                   row.compactness,
                   cell.acc.geom.rowmax
                   ),
                 col.stats=data.frame(
                   col.margins,
                   col.PPV,
                   col.compactness,
                   cell.acc.geom.colmax
                   ),

                 global.stats=data.frame(
                   c(
                     row.Sn.avg=row.Sn.avg,
                     col.PPV.avg=col.PPV.avg,
                     acc.arith=acc.arith,
                     acc.geom=acc.geom,

                     row.Sn.avg.weighted=row.Sn.avg.weighted,
                     col.PPV.avg.weighted=col.PPV.avg.weighted,
                     acc.arith.weighted=acc.arith.weighted,
                     acc.geom.weighted=acc.geom.weighted,

                     row.compactness.avg=row.compactness.avg,
                     col.compactness.avg=col.compactness.avg,
                     geom.compactness=geom.compactness,

                     row.compactness.avg.weighted=row.compactness.avg.weighted,
                     col.compactness.avg.weighted=col.compactness.avg.weighted,
                     geom.compactness.weigthed=geom.compactness.weighted,

                     cell.acc.geom.rowmax.colavg=cell.acc.geom.rowmax.colavg,
                     cell.acc.geom.colmax.rowavg=cell.acc.geom.colmax.rowavg,
                     cell.acc.geom.geom.rowcol=cell.acc.geom.geom.rowcol
                     )
                   )
                 )
  names(result$global.stats) <- "value"

  return(result)
}



contingency.stats.example1 <- function() {

  ## Generate a small contingency table
  x <- t(matrix(c(200,0,0,0,0,
                50,50,0,0,0,
                0,0,20,0,0,
                0,0,0,10,30),ncol=4))

  ## Calculate and display contingency statistics
  res <- contingency.stats(x)
  return (res)
}

contingency.stats.example2 <- function(ncells=20,
                                       ncol=5,
                                       maxval=200
                                       ) {

  ## generate a table of random numbers
  x <- matrix(round(runif(n=ncells,max=maxval)),ncol=ncol)

  ## Calculate and display contingency statistics
  res <- contingency.stats(x)
  return (res)
}
