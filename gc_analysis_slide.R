library(Biostrings) # read and manipulate DNA data
library(ggplot2) #data vis

#download.file('ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/reference/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz',
#              'genome.fna.gz')
filename <- "Genomes_new/Archaea/Asgard_group/C_Heimdallarchaeota/GCA_020348965.1/CHR/CP070760.1.fna"
genome = readDNAStringSet(filename)

# searching DNAStringSet for matches
matchPattern("AAAAAAAAA", genome[[1]], max.mismatch = 2)

# sliding window search
slidingGC = letterFrequencyInSlidingView(genome[[1]], 
                                         view.width=100000,
                                         letters=c('G', 'C'),
                                         as.prob = T)

sampleGC = slidingGC[seq(1, nrow(slidingGC), by=10000),]
sampleGC = as.data.frame(sampleGC)
sampleGC$position = seq(1, nrow(slidingGC), by=10000)
sampleGC$diff = sampleGC$G - sampleGC$C

# plot
ggplot(sampleGC) +
  geom_path(aes(x = position, y = diff)) +
  theme_minimal() + 
  xlab('Genome Position') +
  ylab('P(G) - P(C)')
