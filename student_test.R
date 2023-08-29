## ##############################################################
## Student test for comparing sample means
##
## Author: Jacques van Helden
##
## Running this script requires to first run the script config.R
## Running this script requires to first run the script config.R
##    source('http://pedagogix-tagc.univ-mrs.fr/courses/statistics_bioinformatics/R-files/config.R')
##
##
## Single-command execution:
## source(file.path(dir.R.files, 'student_test.R'))

source(file.path(dir.util, 'util_chip_analysis.R'))
source(file.path(dir.util, 'microarray_util.R'))
source(file.path(dir.util, 'util_student_test_multi.R'))
source(file.path(dir.util, "multitesting_corrections.R"))
## Check the requirement for some packages
pkg <- "multtest"
if (!suppressPackageStartupMessages(require(pkg, quietly=TRUE, character.only = TRUE))) {
  source("http://bioconductor.org/biocLite.R")
  biocLite();
  biocLite(pkg)
}
library(multtest)

## dir.golub <- file.path(dir.data,'gene_expression','golub_1999')
## dir.golub <- file.path(dir.data,'Golub_1999')

dir.results.student <- file.path(dir.results, 'student_multitest')
dir.create(dir.results.student, recursive=TRUE,showWarnings=FALSE)

## Storey's library for the q-value (downloaed from
## http://genomics.princeton.edu/storeylab/qvalue/), which includes an
## estimation of the proportion of true and alternative hypotheses.
library(qvalue)

## color specification (actually I don't use colors for the book
in.colors <- T
if (in.colors) {
  plot.col <- c('lines'='blue',
                'grid'='#CCCCCC',
                'points'='#AAAAAA',
                'positive'='brown'
                )
#  heatmap.palette <- blue.to.yellow()
#  heatmap.palette <- blue.to.red()
#  heatmap.palette <- green.to.red()
  heatmap.palette <- blue.white.red()

} else {
  plot.col <- c('lines'='black',
                'grid'='#CCCCCC',
                'points'='#AAAAAA',
                'positive'='black'
                )
  heatmap.palette <- gray((0:255)/255)
}


## ##############################################################
## Golub dataset. Classify cell types according to gene expression
## profiles
## load the data
data(golub)
dim(golub)
## xpr <- data.frame(cbind (golub.gnames, golub))
## names(xpr) <- c("ID", "description", "name", 1:ncol(golub))

## tissue types
## 27 patients with acute lymphoblastic leukemia (ALL)
## and 11 patients with acute myeloid leukemia (AML)
print(golub.cl)
golub.cl[golub.cl==0] <- "ALL"
golub.cl[golub.cl==1] <- "AML"

## set the names to genes in the matrix
row.names(golub) <- golub.gnames[,3]

## Perform a chip-wise standardization
golub.standard <- standardize.chips(golub, centr.est='median', disp.est='iqr')
names(golub.standard)

## ##############################################################
## t.test with a single gene
g <- 347
x <- golub.standard$z[g,]

## separate data in two vectors
x.ALL <- x[golub.cl=="ALL"]
x.AML <- x[golub.cl=="AML"]
t.test(x.ALL,x.AML)

## The same test, but use a formula instead of the two separate samples
t <- t.test(x ~ golub.cl, na.action=na.omit)
print(t) ## print short report

names(t)
attributes(t)

## ##############################################################
## Perform Weoch t-test on each row in parallel

## apply the t.test to each gene of the data set
t.table <- t.test.multi(golub.standard$z,golub.cl, volcano.plot=T, robust.est=F,main='Welch test - Golub (1999)',plot.col=plot.col)
dim(t.table)
names(t.table)
setwd(dir.results.student); export.plot(file.prefix='golub_volcano_plot', export.formats=export.formats.plots, width=12,height=12)

## ## Draw standard error as a function of Welch significance
## x11(width=6,height=6)
## plot(t.table[,c('st.err.diff','sig')],
##      xlab='Standard error of mean difference',
##      ylab='Welch significance (sig=-log10(E-value)',
##      col=plot.col['points'],
##      panel.first=grid(lty='solid',col=plot.col['grid']))
## abline(h=0,col=plot.col['lines'],lwd=2)
## points(t.table[t.table$sig > 0,c('st.err.diff','sig')],col=plot.col['positive'])
## setwd(dir.results.student); export.plot(file.prefix='golub_sd_sig', export.formats=export.formats.plots, width=6,height=6)

## ## Plot standard error versus mean differences, with circles
## ## proportional to significance
## x11(width=6,height=6)
## plot(t.table[,c('means.diff', 'st.err.diff')],
##      cex=t.table$sig-min(t.table$sig),
##      ylab='Standard error on the difference',
##      xlab='Difference between the means',
##      col=plot.col['points'],
##      panel.first=grid(lty='solid',col=plot.col['grid']))
## abline(v=0)
## points(t.table[t.table$sig > 0,c('means.diff', 'st.err.diff')],
##        cex=t.table[t.table$sig > 0,'sig']-min(t.table$sig),
##        col=plot.col['positive'])
## setwd(dir.results.student); export.plot(file.prefix='golub_diff_sd_mean_sig', export.formats=export.formats.plots, width=6,height=6)

## ################################################################
## ## Apply various multi-testing corrections to the T-test table.
golub.multitest <- multitest.corrections(pval=t.table$P.value, names=rownames(t.table), file.prefix=file.path(dir.results.student, 'golub_multitest_corrections'), lambda=0.5, main="Golub")
t.table <- cbind(t.table, golub.multitest$multitest.table)
print(as.data.frame(golub.multitest$nb.signif))

################################################################
## Sort genes according to the P-value of the t-test
sorted.genes <- row.names(t.table[order(t.table$P.value),])
golub.sorted <- golub[sorted.genes,]
t.table.sorted <- t.table[sorted.genes,]
t.table.sorted$rank <- 1:nrow(t.table.sorted)
golub.standard.sorted <- standardize.chips(golub.sorted,centr.est='median',disp.est='iqr')

## select expression profiles for the genes with differential expression
selected <- row.names(t.table.sorted[t.table.sorted$E.value < 1,]) # names of the regulated genes
x.sorted <- golub.standard.sorted$z[selected,]
dim(x.sorted)


################################################################
## Export the results

## export the filtered genes
write.table(x.sorted,file.path(dir.results.student,'golub_standardized_filtered_student_E1.tab'),sep='\t',quote=F)
write.table(golub.standard.sorted$z,file.path(dir.results.student, 'golub_standardized_3051.tab'),sep='\t',quote=F)
write.table(golub,file.path(dir.results.student, 'golub_3051.tab'),sep='\t',quote=F)

## Export tables with t-test statistics + multi-testing corrections
write.table(t.table,file.path(dir.results.student, 'golub_standardized_student_3051.tab'),sep='\t',quote=F)
write.table(t.table.sorted,file.path(dir.results.student, 'golub_standardized_student_3051_sorted.tab'),sep='\t',quote=F)

## Export all genes sorted by t.obs value. Beware, this will put on
## the top the genes over-expressed in group A versus B, and on the
## bottom the genes overexpressed in group B versus A.
names.sorted.by.tobs <- row.names(t.table[order(t.table$t.obs),])
write.table(golub.standard$z[names.sorted.by.tobs,],file.path(dir.results.student, 'golub_zscores_3051_sorted_by_tobs.tab'),sep='\t',quote=F)

################################################################
## Standardize by gene using robust estimators (median, IQR)
##
## Note: I generally do NOT recommend gene-wise standardization,
## because it looses the information about gene-specific variance. In
## this particular case, I do it for the purpose of visual
## illustration, to generate a heatmap where genes will be sorted by
## increasing value of the t.obs, in order to progressively go from
## ALL-expressed to AML-expressed genes.

m.est <- apply(golub.standard$z,1,median,na.rm=T)
iqr <- apply (golub.standard$z, 1, IQR, na.rm=T)
s.est <- iqr/(qnorm(0.75) - qnorm(0.25))
golub.z.genewise.stand <- (golub.standard$z - m.est ) / s.est

## Check that the standardized genes same median (should be 0) and IQR (should be  1.34898)
range(apply(golub.z.genewise.stand,1,median))
range(apply(golub.z.genewise.stand,1,IQR))
write.table(golub.z.genewise.stand[names.sorted.by.tobs,],file.path(dir.results.student, 'golub_zscores_3051_sorted_by_tobs_genewise_stand.tab'),sep='\t',quote=F)

## Draw a heatmap of profiles sorted by T-test statistics (t.obs)
x11(width=8,height=8)
heatmap(golub.z.genewise.stand[names.sorted.by.tobs,],
        main='Golub, genes sorted by t.obs',
        Rowv=NA, Colv=NA, zlim=c(-2.5,2.5), col=heatmap.palette,
        labRow=NA,labCol=golub.cl
        )
setwd(dir.results.student); export.plot(file.prefix='golub_profiles_sorted_by_tobs', export.formats=export.formats.plots, width=6,height=6)


################################################################
## Negative control= run same analysis with permuted values

## Permute Golub's data and apply Welch test
perm.data <- data.frame(matrix(sample(golub.standard$z, replace=F), nrow=nrow(golub.standard$z), ncol=ncol(golub.standard$z)))
x11(width=6,height=6)
t.table.perm <- t.test.multi(perm.data,golub.cl, volcano.plot=T, robust.est=F,main='Welch test - randomized data',plot.col=plot.col)
setwd(dir.results.student); export.plot(file.prefix='perm_data_volcano_plot', export.formats=export.formats.plots, width=6,height=6)


## Perform multitesting corrections on Welch result of permuted data
perm.multitest <- multitest.corrections(pval=t.table.perm$P.value, file.prefix=file.path(dir.results.student, 'permuted_data'), main='Permuted data')
t.table.perm <- cbind(t.table.perm, perm.multitest$multitest.table)
print(as.data.frame(perm.multitest$nb.signif))

## Gene-wise standardization of permuted profiles
perm.m.est <- apply(perm.data,1,median,na.rm=T)
perm.iqr <- apply (perm.data, 1, IQR, na.rm=T)
perm.s.est <- perm.iqr/(qnorm(0.75) - qnorm(0.25))
perm.genewise.stand <- (perm.data - perm.m.est ) / perm.s.est

## Check that the standardized genes same median (should be 0) and PERM.IQR (should be  1.34898)
range(apply(perm.genewise.stand,1,median))
range(apply(perm.genewise.stand,1,IQR))

## Export the table sorted by t.obs
write.table(perm.genewise.stand[order(t.table.perm$t.obs),],file.path(dir.results.student, 'perm_data_3051_sorted_by_tobs_genewise_stand.tab'),sep='\t',quote=F)

## Draw a heatmap of profiles sorted by T-test statistics (t.obs)
x11(width=8,height=8)
heatmap(as.matrix(perm.genewise.stand[order(t.table.perm$t.obs),]),
        main='Permuted data, genes sorted by t.obs',
        Rowv=NA, Colv=NA, zlim=c(-2.5,2.5), col=heatmap.palette,
        labRow=NA,labCol=golub.cl
        )
setwd(dir.results.student); export.plot(file.prefix='perm_profiles_sorted_by_tobs', export.formats=export.formats.plots, width=6,height=6)

