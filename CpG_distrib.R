################################################################
## Compute the distribution of CpG dinucleotides in all human
## promoters, to highlight the bimodal distribution due to general CpG
## avoidance except in CpG islands.
##
## Author: Jacques van Helden
##
## Running this script requires to first run the script config.R
##    source('http://pedagogix-tagc.univ-mrs.fr/courses/statistics_bioinformatics/R-files/config.R')
##
## Single-command execution:
##   source(file.path(dir.R.files, 'CpG_distrib.R'))

source(file.path(dir.R.files, "util", "util_central_tendency.R"))

## The dinucleotide occurrence file was computed with the RSAT programs retrieve-seq and
## oligo-analysis. It contains 1 row per promoter, and 1 column per
## dinucleotide.
oligo.file <- file.path(dir.data, 'oligo_frequencies','human_CpG','2nt-1str-noov_Homo_sapiens_EnsEMBL_allup500_mRNA.tab')
oligos <- read.delim(file=oligo.file,sep='\t',comment.char=';', header=1,row.names=NULL)

oligo <- "aa"

X11(width=7,height=5)

for (oligo in names(oligos)[2:ncol(oligos)]) {
  occ <- oligos[,oligo]

  hist(occ,breaks=50,col='#BBBBBB',
       main=paste('Occurrences of',toupper(oligo),'dinucleotides in human promoters (500bp)'),
       xlab='Occurrences per promoter',
       ylab='Number of promoters',
       color='#BBBBBB',border='#888888'
       )

  ## Export the plot
  setwd(dir.figures); export.plot(file.prefix=paste("histo_human_promoters_",oligo,"_occ",sep=""), export.formats=export.formats.plots, width=7,height=5)
  
  ##############################################################
  ## Draw histogram with central tendency estimators
  central.tendency.histo(x=occ,
                         class.interval=2,
                         main=paste('Occurrences of',toupper(oligo),'dinucleotides in human promoters (500bp)'),
                         xlab='Occurrences per promoter',
                         ylab='Number of promoters',
                         col='#DDDDDD', border='#888888',
                         stat.colors=c(mean="black",median='black',mode='black'),
                         stat.lty=c(mean="solid",median='dashed',mode='dotted'),
                         legend.box=TRUE
                         )

  ## Export the plot
  setwd(dir.figures); export.plot(file.prefix=paste("central_tendency_human_promoters_",oligo,"_occ",sep=""), export.formats=export.formats.plots, width=7,height=5)

}
