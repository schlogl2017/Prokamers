################################################################
## Illustration of the Markov chains
source(file.path(dir.util, 'util_plots.R'))

################################################################
## Define a function to compute one transition matrix and plot the
## corresponding Hinton diagram

one.plot <- function(organism,
                     k = 4, ## k-mer length
                     type = 'genomic', ## Sequence type (genomic | upstream)
                     from = -500, ## Distal limit for the sequences
                     to = -1, ## Proximal limit for upstream sequences
                     shape='bars', ## Shape fr the Hinton diagram
                     ... ## All other parameters are passed to plot.hinton.doagram()
                     ) {

  m <- k -1 ## Markov order

  ## Directories and files
  dir.markov <- file.path(dir.data, 'genome_analysis', 'oligo_frequencies')
  if (type == 'genomic') {
    file.markov <- file.path(dir.markov, paste(k, 'nt_genome_',organism,'-ovlp-1str.freq', sep=''))
  } else {
    file.markov <- file.path(dir.markov, paste(k, 'nt_upstream-noorf_',organism,'_up',from,'to',to,'-ovlp-1str.freq', sep=''))
  }

  ## Load the k-mer frequencies
  x <- read.delim(file.markov, comment.char=';', header=1,row.names=NULL,as.is=T)
  names(x) <- c('pattern', 'id', 'freq', 'occ')

  ## Compute prefixes and suffixes
  prefixes <- sort(unique(substring(x$pattern,1,last=k-1)))
  suffixes <- sort(unique(substring(x$pattern,k)))

  ## Fill a table with pattern frequencies
  frequencies <- data.frame(matrix(nrow=length(prefixes), ncol=length(suffixes), NA))
  names(frequencies) = suffixes
  row.names(frequencies) = prefixes
  for (i in 1:nrow(x)) {
    pattern <- x[i, "pattern"]
    prefix <- substring(pattern, 1, last=k-1)
    suffix <- substring(pattern, k)
    frequencies[prefix,suffix] <- x[i, "freq"]
  }

  ## Compute marginal frequencies
  prefix.proba <- apply(frequencies, 1, sum)
  suffix.proba <- apply(frequencies, 2, sum)
  transitions <- frequencies / prefix.proba

                                        # transitions$prefix.proba <- prefix.proba
#  transitions["suffix",] <- suffix.proba

  ## Draw a Hinton diagram with the transition frequencies
  plot.height=min(12,4^(k-1)+2)
  X11(width=4,height=plot.height)
  plot.hinton.diagram(transitions,main=paste(organism, ";", type),xlab="suffix",ylab="prefix",shape=shape,fill.colors=c("grey","white"), ...)

  ## Export the plot
  file.name <- paste(organism, type, sep='')
  if (type == "upstream") {
    file.name <- paste(file.name, "_from", from, "_to", to, sep='')
  }
  file.name <- paste(file.name, "_markov", m, sep='')
  setwd(dir.figures); export.plot(file.prefix=file.name, export.formats=c("png","eps"), width=4,height=plot.height)

  return(transitions)
}


## Parameters
organisms <- c('Saccharomyces_cerevisiae','Escherichia_coli_K12')
for (org in organisms) {
  for (k in 1:4) {
    one.plot(organism = org, k=k, type='genomic')
    dev.off()
    one.plot(organism = 'Saccharomyces_cerevisiae', k=3, type='upstream', from=-500, to=-1)
    dev.off()
  }
}



