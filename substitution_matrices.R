################################################################
## Draw a hinton diagrm with the PAM250 substitution matrix
##
## Author: Jacques van Helden
##
## Running this script requires to first run the script config.R
##
## Single-command execution:
## source(file.path(dir.R.files, 'substitution_matrices.R'))

source(file.path(dir.util, 'util_plots.R'))

in.colors <- F
if (in.colors) {
  plot.colors <- c(positive='#000088', negative='#FFFF00')
} else {
  plot.colors <- c(positive='black', negative='#BBBBBB')
}

dir.data.subst <- file.path(dir.data, "sequence_alignment", "substitution_matrices")
dir.res.subst <-  file.path(dir.results,  "substitution_matrices")

X11(width=8,height=8)
pam250 <- read.delim(file.path(dir.data.subst, 'PAM250_log-odds_matrix.tab'), row.names=1, header=1)
plot.hinton.diagram(pam250[1:20,1:20], line.color=c(plot.colors['positive'], "gray"), fill.color=c(plot.colors['positive'],plot.colors['negative']), main="PAM250 log-odds",shape="rectangle")
setwd(dir.figures); export.plot(file.prefix="pam250_hinton_diagram", export.formats=export.formats.plots, width=8,height=8)

X11(width=8,height=8)
blosum62 <- read.delim(file.path(dir.data.subst, 'BLOSUM62_log-odds_matrix.tab'), row.names=1, header=1)
plot.hinton.diagram(blosum62[1:20,1:20], line.color=c(plot.colors['positive'], "gray"), fill.color=c(plot.colors['positive'],plot.colors['negative']), main="BLOSUM62 log-odds",shape="rectangle")
setwd(dir.figures); export.plot(file.prefix="blosum62_hinton_diagram", export.formats=export.formats.plots, width=8,height=8)

 
## Plot a Hinton diagram with black background for the book cover
X11(width=6,height=6)
par(bg='black')
par(col.axis='white')
par(col.lab='white')
par(col='white')
par(col.main='white')

plot.hinton.diagram(pam250[1:20,1:20], line.color=c("#FFFF88", "#FF4422"), fill.color=c("#FFFF88","#FF4422"), main="PAM 250 log-odds",shape="rectangle")
setwd(dir.figures); export.plot(file.prefix="pam250_hinton_diagram_cover", export.formats=export.formats.plots, width=6,height=6)

plot.hinton.diagram(blosum62[1:20,1:20], line.color=c("#FFFF88", "#FF4422"), fill.color=c("#FFFF88","#FF4422"), main="BLOSUM62 log-odds",shape="rectangle")
setwd(dir.figures); export.plot(file.prefix="blosum62_hinton_diagram_cover", export.formats=export.formats.plots, width=6,height=6)

par(bg='white')
par(col.axis='black')
par(col.lab='black')
par(col='black')
par(col.main='black')

################################################################
## PAM elementary probabilities

## Load residue frequencies
aa.normalized.freq <- as.matrix(read.delim(file.path(dir.data.subst, 'aa_normalized_freq_dayhoff_1978.tab'), row.names=1, header=1))

## Load residue mutability
aa.mutability <- as.matrix(read.delim(file.path(dir.data.subst, 'relative_mutabilities_dayhoff_1978.tab'), row.names=1, header=1))



## Load PAM001 substitution rates from Dayhoff articles. Beware:
## diagonal elements are missing but can be computed (the sum of each
## column must be 10,000).
pam001.freq <- read.delim(file.path(dir.data.subst, 'PAM001_freq_matrix.tab'), row.names=1, header=1)
pam001.freq <- as.matrix(pam001.freq[,2:21]) ## Truncate second column containing the residue single-letter code
subst.per.aa <- apply(pam001.freq, 2, sum, na.rm=T)

lambda <- 0.00013310 ## I don't remember where I got this value from !!!
#lambda <- 0.0001 ## I don't remember where I got this value from !!!

pam001.proba <- pam001.freq
for (i in 1:20) {
  ## I don't understand why I have to work on rows (rather than
  ## columns) to obtain the same result as in the Excel sheet
  ## "PAM_matrices.xls"
  pam001.proba[i,] <- 10000 * lambda * pam001.freq[i,] * as.vector(aa.mutability) / as.vector(subst.per.aa)
  pam001.proba[i,i] <- 10000*(1 - lambda*aa.mutability[i,1])
}
pam001.proba <- pam001.proba/10000
apply(pam001.proba, 2, sum)
range(pam001.proba)

################################################################
## Load the PAM001 probability matrix (has to be rescaled by a factor of 1/10000)
pam001.proba.read <- read.delim(file.path(dir.data.subst, 'PAM001_proba_matrix.tab'), row.names=1, header=1)
pam001.proba.read <- as.matrix(pam001.proba.read[2:21]) / 10000

## Compare the matrix read from the article and the one computed above
max(pam001.proba.read - pam001.proba)
plot(pam001.proba.read, pam001.proba, log='xy')
plot(pam001.proba.read, pam001.proba/pam001.proba.read, log='xy')

## plot a Hinton diagram to show the differences between the computed and published matrices
pam001.diff <- (pam001.proba-pam001.proba.read)
plot.hinton.diagram(pam001.diff, line.color=c(plot.colors['positive'], "gray"), fill.color=c(plot.colors['positive'],plot.colors['negative']), main="PAM250 log-odds",shape="rectangle")




dim(pam001.proba)
sum(pam001.proba)
apply(pam001.proba, 2, sum)

diag(pam001.proba)

## Compute percent of identity wth a single matrix product
percent.identity <- as.vector(as.matrix(aa.normalized.freq)) %*% diag(pam001.proba)

## Expose the detail of the computation in a data frame
identity.detail <- data.frame('freq' = as.vector(as.matrix(aa.normalized.freq)),
                              'cons.percent' = diag(pam001.proba))
identity.detail$product <- (identity.detail$freq * identity.detail$cons.percent)
percent.identity.check <- sum(identity.detail$product)

################################################################
## Extrapolate PAM matrix series from PAM001

pam002.proba <- pam001.proba %*% pam001.proba
apply(pam001.proba, 2, sum)


## TO BE DONE. I HAD THE RESULT IN THE BOOK, AND I CAN'T FIND THE R CODE ANYMORE

