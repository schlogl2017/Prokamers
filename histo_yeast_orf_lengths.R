################################################################
#
# Figure: histogram of yeast orf lengths
#

## Load the data
length.file <- file.path(dir.data, "orf_lengths", "yeast_orf_length_enum.txt")
lengths	<- scan(length.file)

## Compute class number
class.interval <- 300
class.number <- as.integer(max(lengths)/class.interval) + 1

## Draw the histogram
X11(width=6,height=4)
yeast.hist <- hist(lengths,
                   breaks = class.interval*(0:class.number),
                   col='#DDDDDD',
                   main='Length distribution of yeast protein-coding genes',
                   xlab=paste("Gene length (class interval = ", class.interval, ")"),
                   ylab="Number of genes"
                   )

## Export the plot
setwd(dir.figures); export.plot(file.prefix="histo_yeast_coding_gene_lengths", export.formats=export.formats.plots, width=6,height=4)
