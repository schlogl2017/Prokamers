
################################################################
##
## The table below provides the frequencies of start and stop codons
## in genomic, coding and intergenic sequences respectively.
## Calculate the probability to observe, in each sequence type, an
## open reading frame of
## a) at least 30 bp
## b at least 300 bp
## c) at least 999 bp
##
## Seq	Codon	Genomic	Coding	Intergenic
## ATG	Start	0.01825	0.01868	0.01706
## TAA	Stop	0.02238	0.01991	0.02900
## TAG	Stop	0.01289	0.01246	0.01408
## TGA	Stop	0.02012	0.02118	0.01738
##
## If we observe an open reading frame of at least 300bp, what is the
## probability to be in a coding region, knowing that 72% of the
## genome is coding ?

## Specify the frequencies of start and stop codons in different sequence types
p.start <- c(genomic=0.01825,coding=0.01868,intergenic=0.01706)
p.TAA <- c(genomic=0.02238,coding=0.01991,intergenic=0.02900)
p.TAG <- c(genomic=0.01289,coding=0.01246,intergenic=0.01408)
p.TGA <- c(genomic=0.02012,coding=0.02118,intergenic=0.01738)
p.stop <- p.TAA + p.TGA + p.TAG

## Calcualte the probability to find an ORF of at least L bp at a given position.
L <- c(30,300,999) ## ORF length in nucleotides
n <- L/3 ## Number of codons (including the start and stop)

result <- cbind(L=L,
                n=n,
                "n-2"=n-2,
                t(data.frame(p.ORF.30 = p.start*(1-p.stop)^(L[1]/3-2),
                             p.ORF.300 = p.start*(1-p.stop)^(L[2]/3-2),
                             p.ORF.999 = p.start*(1-p.stop)^(L[3]/3-2)
                             ))
                )
result <- data.frame(result)
row.names(result)=L
names(result) <- c("L", "n", "n-2", "P(ORF|genome)", "P(ORF|coding)", "P(ORF|inter)")

setwd(dir.results); export.object(result, "ORF_len_proba", export.format='table')

################################################################
## Knowing that coding sequences represent 72% of the genome, estimate
## the probability to be in a coding region, if we observe an ORF
## covering at least 300bp?

p.coding <- c(prior=0.72)
p.intergenic = 1 - p.coding

p.ORF.30["prior"] <- p.ORF.30["coding"]*p.coding["prior"] + p.ORF.30["intergenic"]*p.intergenic["prior"]
p.ORF.300["prior"] <- p.ORF.300["coding"]*p.coding["prior"] + p.ORF.300["intergenic"]*p.intergenic["prior"]
p.ORF.999["prior"] <- p.ORF.999["coding"]*p.coding["prior"] + p.ORF.999["intergenic"]*p.intergenic["prior"]

p.coding["ORF30"] <- p.ORF.30["coding"]*p.coding["prior"]/p.ORF.30["prior"]
p.coding["ORF300"] <- p.ORF.300["coding"]*p.coding["prior"]/p.ORF.300["prior"]
p.coding["ORF999"] <- p.ORF.999["coding"]*p.coding["prior"]/p.ORF.999["prior"]

print(p.coding)
