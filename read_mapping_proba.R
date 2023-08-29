################################################################
##
## Illustrations of the binomial distribution, by an example inspired
## from the read mapping procedure commonly used in "Next generation
## sequencing".
##
## Author: Jacques van Helden
##
## Running this script requires to first run the script config.R
##    source('http://pedagogix-tagc.univ-mrs.fr/courses/statistics_bioinformatics/R-files/config.R')
##
## Single-command execution:
##   source(file.path(dir.R.files, 'ngs_read_mapping', 'read_mapping_proba.R'))



################################################################
## A file containing 10,000,000 raw reads 36p each is mapped against
## the Human genome (23 chromosomes, totaling ~3Gb). Before chosing the number of accepted
## mismatches for the mapping, we would like to estimate the number of
## matches expected at random as a function of the accepted number of
## mismatches. To simplify the calculation, we assume that genomic
## sequences follow a Bernoulli schema withe quiprobable nucleotides.

## Parameters of the problem
k <- 36 ## Read length
s <- 0:36 ## Number of accepted mismatches
p <- 0.25 ## Matching probability for one nucleotide, assuming equiprobability and independence
q <- 1 - p ## Mismatch probability for one nucleotide
x <- k -s ## Number of matches
N <- 10000000 ## Number of reads in the file set
chrom <- 23 ## number of chromosomes
Gb <- 3 ## Genome size in Gb
G <- Gb * 1e+9 ## Genome size in nucleotides

## Count the total number of positions against which each read will be
## matched. Discard the k-1 positions at the end of each chromosome
T <- 2 * (G - chrom *(k-1))

## Probability of perfect match
perfect.match.proba <- p^k

## Probability of at most one substitution
p.one.subst <- k*p^(k-1)*(1-p)
p.at.most.one <- p^k + p.one.subst

## Probability of exactly s substitutions
proba.exact.s.subst <- dbinom(s, prob=q, size=k)

## Matching probability for at most s substitutions. We use the increasing CDF = P (S <= s).
pval.per.read <- pbinom(s, size=k, prob=q, lower=TRUE)

## Compute the e-value for a genome-scale search of a single read
eval.per.read <- pval.per.read * T

## E-value for genome-scale search of N reads
eval.total <- eval.per.read * N

## FWER for a single read versus the genome.  FWER indicates the
## probability to find at least one match per chance in the genome.
fwer.per.read <- pbinom(q=0, prob=pval.per.read, size=T, lower.tail=FALSE)

## FWER for the whole analysis: probability to obtain by chance at
##  least one match in the whole genome for at least one read
fwer.total <- pbinom(q=0, prob=pval.per.read, size=T*N, lower.tail=FALSE)

## Export result table
results <- data.frame(s=s,x=x,
                      pval=pval.per.read,
                      eval.per.read=eval.per.read,
                      fwer.per.read=fwer.per.read,
                      eval.total = eval.total,
                      fwer.total = fwer.total
                      )
setwd(dir.results); write.table(results, file='ngs_read_mapping_proba.tab', row.names=FALSE, quote=FALSE)
setwd(dir.results); write.table(signif(results, digits=3), file='ngs_read_mapping_proba_rounded.tab', row.names=FALSE, quote=FALSE, , sep= "\t", eol="\n")
setwd(dir.results); write.table(signif(results, digits=3), file='ngs_read_mapping_proba_rounded.tex', row.names=FALSE, quote=FALSE, , sep= " & ", eol="\\\\\n")

## Draw a plot with matching probabilities for a single read
x11(width=7, height=5)
plot(s, pval.per.read, type="l", panel.first=grid(), ylim=c(0,1),
     main=paste(sep='', 'Matching probability for a ', k,'bp read'),
     xlab='Accepted substitutions (s)',
     ylab='Probabilities', lwd=2
     )
lines(s, proba.exact.s.subst, lty='solid', lwd=4, type='h', col='#888888')
legend("topleft", legend=c("exactly s substitutions", "at most s substitutions"), lwd=c(4,2),
       lty=c("solid", "solid"),
       col=c("#888888", "black"),
       bg="white", bty="o")
setwd(dir.figures); export.plot(file.prefix='read_mapping_proba', export.formats=export.formats.plots, width=7,height=5)


## Draw a plot with matching probabilities for a single read, using log scale for probabilities
x11(width=7, height=5)
plot(s, pval.per.read, type="l", panel.first=grid(), ylim=c(min(proba.exact.s.subst),1),
     main=paste(sep='', 'Matching probability for a ', k,'bp read'),
     xlab='Accepted substitutions (s)',
     ylab='Probabilities (log scale)', lwd=2, log="y"
     )
lines(s, proba.exact.s.subst, lty='dashed', lwd=2, type='l', col='#888888')
legend("topleft", legend=c("exactly s substitutions", "at most s substitutions"), lwd=c(2,2),
       lty=c("dashed", "solid"),
       col=c("#888888", "black"),
       bg="white", bty="y")
setwd(dir.figures); export.plot(file.prefix='read_mapping_proba_logy', export.formats=export.formats.plots, width=7,height=5)


## ## Draw a plot with e-values
## x11(width=7, height=5)
## plot(s, eval.total, type="l",
##      panel.first=grid(), 
##      main=paste(sep='', 'Read mapping - total e-value (N=', N, ', k=', k,', G=', Gb, 'Gb)'),
##      xlab='Accepted substitutions (s)',
##      ylab='Number of matches expected at random (log scale)', lwd=2, log="y"
##      )
## legend("topleft", legend=c("e-value"), lwd=c(2,2),
##        lty=c("solid"),
##        col=c("black"),
##        bg="white", bty="y")
## setwd(dir.figures); export.plot(file.prefix='read_mapping_eval_logy', export.formats=export.formats.plots, width=7,height=5)

## Draw a plot with p-value + read-wise evalue + total e-value
ymax <- 1e+17
ymin <- 1e-23
scalebars <- 17+23+1

x11(width=7, height=7)
par(mar=c(5,6,4,2))
plot(s, pval.per.read, type="p", pch=1, panel.first=grid(lty='solid', col='#BBBBBB'),
     ylim=c(ymin, ymax),
     main=paste(sep='', 'Random expectation for read mapping (N=', N, ', k=', k,', G=', Gb, 'Gb)'),
     xlab='Accepted substitutions (s)',
     ylab='', log="y", las=1)
title(ylab='Random expectation (log scale)', mgp=c(4,1,0))
abline(h=1, col='black', lwd=2)
axis(2, at = 10^(-23:20), labels=FALSE)
lines(s, eval.per.read, pch='+', type='p')
lines(s, eval.total, pch=2, type='p')
legend("bottomright", legend=c("p-value", "read-wise e-value", "total e-value"), 
       pch=c(1,3,2),
       bg="white", bty="y")
setwd(dir.figures); export.plot(file.prefix='read_mapping_pval_eval', export.formats=export.formats.plots, width=7,height=7)

