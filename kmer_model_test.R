source("kmer_model_functions.R")

filename = "Jelly/Bacteria/Pseudomonadota/Acidithiobacillia/CHR/GCA_900174455.1_k1_10.csv"
name <- strsplit(filename, split="/")[[1]][6]
name <- substr(name, start = 1, stop = 15)

#len_seqs <- read.csv("Results/Lengths/Archaea/Asgard_chr_lengths.tsv", 
                     #header=FALSE, sep = "\t", col.names = c("ID", 'length'))

#len_seq <- len_seqs[len_seqs$ID == name, 'length']

kmer_data <- read.csv(filename, 
                      header=TRUE, 
                      sep = ",", 
                      )
head(kmer_data)
tail(kmer_data)

# bases_freq_from_counts <-function(kmer_data, 
#                                   k=1, 
#                                   cols=c("kmer", 
#                                          "observed_freq")){
#   nucs <- c(subset(kmer_data, 
#                    nchar(kmer) == 1, 
#                    select=cols[1]))[[1]]
#   freq <- c(subset(kmer_data, 
#                    nchar(kmer) == 1, 
#                    select=cols[2]))[[1]]
#   names(freq) <- nucs
#   return(freq)
# }

b <- bases_freq_from_counts(kmer_data)


# kmer_priori <- function(base_freqs, kmer){
#   k <- nchar(kmer)
#   ## calculate word probability
#   pr <- vector()
#   pr[kmer] <- 1 # initialize
#   for (i in 1:k) {
#     nuc <- substr(kmer,i,i) ## Select the letter at position i of word W
#     pr[kmer] <- pr[kmer] * base_freqs[nuc] ## update word probability
#   }
#   return(pr[kmer])
# }

aatg <- kmer_priori(b, "AATG")

#which(kmer_data$kmer == "AATG")

df <- calculate_kmer_probs(kmer_data, 
                           columns=c("kmer", 
                                     "observed",
                                     "observed_freq"),
                           k=4)

sum(df$priori_prob)
sum(df$observed_freq)
