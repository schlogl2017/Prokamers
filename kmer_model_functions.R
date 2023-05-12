# Functions to use in kmer analysis
# Author: Paulo Sérgio Schlögl
# last modified May, 2023
# first written Apr, 2022

# Reading the data and subseting the data
select_data <- function(filename, k){
  # select the kmer data of length k from a big csv
  # with two columns: kmers, counts and frequencies
  data <- read.csv(filename, 
                   header=TRUE, 
                   stringsAsFactors=FALSE)
  return(subset(data, nchar(kmer) == k))
}
  
select_data_counts <- function(counts, columns, k){
  df_obs <- subset(counts, 
                   nchar(kmer) == k, 
                   select=columns)
  return(df_obs)
}
  
# Calculate the Model depending on kmer length (k)
# kmer lemgth <= 2
expected_zero_order <- function(kmer, kmer_counts){
  # A simple Bernoulli model
  n <- nchar(kmer)
  # subsetting the dataframe where A,C,G,T ave counts and sum the values
  # thats thats is the sequence length
  seq_len <- sum(subset(kmer_counts$observed, nchar(kmer_counts$kmer) == 1))
  N <- seq_len - n + 1
  k_prob <- kmer_prob(kmer_counts, kmer)
  kmer_exp <- N * k_prob
  kmer_var <- kmer_exp * (1 - kmer_exp / N)
  return(list(kmer = kmer, expected = kmer_exp, variance = kmer_var))
}

get_kmer_expected_var_ZOM <- function(kmer_counts, columns, k){
  # Calculates the expected kmer values based in a HMO of k - 2
  df <- select_data_counts(kmer_counts, columns, k)
  n <- nrow(df)
  expected <- data.frame(kmer = character(), expected = numeric(0), variance = numeric(0))
  for(i in 1:n){
    expected[i, ] <- expected_zero_order(df$kmer[i], kmer_counts)
  }
  merged_df <- merge(df, expected, by = "kmer")
  return(merged_df)
}
    
# Higher Markove Orders
expected_high_markov_var <- function(kmer, kmer_counts){
  # Function to calculate the expected values of counted kmers
  # the model is a higher order markov k - 2
  n <- nchar(kmer)
  seq_len <- sum(subset(kmer_counts$observed, nchar(kmer_counts$kmer) == 1))
  N <- seq_len - n + 1
  pref <- substr(kmer, 1, n - 1)
  p <- as.double(kmer_counts[kmer_counts$kmer == pref, 'observed'])
  suf <- substr(kmer, 2, n)
  s <- as.double(kmer_counts[kmer_counts$kmer == suf, 'observed'])
  mid <- substr(kmer, 2, n - 1)
  m <- as.double(kmer_counts[kmer_counts$kmer == mid, 'observed'])
  kmer_exp <- (p * s) / m
  kmer_var <- kmer_exp * (1 - kmer_exp / N)
  return(list(kmer = kmer, expected = kmer_exp, variance = kmer_var))
}
    
get_kmer_expected_var_HMO <- function(kmer_counts, columns, k){
  # Calculates the expected kmer values based in a HMO of k - 2
  df <- select_data_counts(kmer_counts, columns, k)
  n <- nrow(df)
  expected <- data.frame(kmer = character(), 
                         expected = numeric(0), 
                         variance = numeric(0))
  for(i in 1:n){
    expected[i, ] <- expected_high_markov_var(df$kmer[i], kmer_counts)
  }
  merged_df <- merge(df, expected, by = "kmer")
  return(merged_df)
}

# Utilities for Statistics
kmer_frequency <- function(df, k, col){
  # Calculates the kmer frquency of counted and expected kmers values
  # create the column names to the new column
  new_name <- paste(col, "_", "freq", sep = "", collapse = "")
  num_rows <- nrow(df)
  total <- sum(df[, col])
  df[, "frequency"] <- df[, col] / total
  names(df)[names(df) == "frequency"] <- new_name
  return(df)
}

bases_freq_from_counts <-function(kmer_data, 
                                  k=1, 
                                  cols=c("kmer", 
                                         "observed_freq")){
  nucs <- c(subset(kmer_data, 
                   nchar(kmer) == 1, 
                   select=cols[1]))[[1]]
  freq <- c(subset(kmer_data, 
                   nchar(kmer) == 1, 
                   select=cols[2]))[[1]]
  names(freq) <- nucs
  return(freq)
}
   

kmer_priori <- function(base_freqs, kmer){
  k <- nchar(kmer)
  ## calculate word probability
  pr <- vector()
  pr[kmer] <- 1 # initialize
  for (i in 1:k) {
    nuc <- substr(kmer,i,i) ## Select the letter at position i of word W
    pr[kmer] <- pr[kmer] * base_freqs[nuc] ## update word probability
  }
  return(pr[kmer])
}

calculate_kmer_probs <- function(kmer_data,
                                 columns, k){
  bases <- bases_freq_from_counts(kmer_data)
  data <- select_data_counts(kmer_data, 
                             columns = columns, 
                             k)
  n <- nrow(data)
  
  priori <- data.frame(kmer = character(),
                       prob = numeric(0))
  for(i in 1:n){
    kmer <- data$kmer[i]
    pri <- kmer_priori(bases, kmer)
    #cat(kmer, " ", pri,"\n")
    priori[i, ] <- c("kmer" = kmer, prob = pri)
  }
  colnames(priori) <- c("kmer", "priori_prob")
  df <- merge(data, priori, by = "kmer")
  df$priori_prob <- as.numeric(df$priori_prob) 
  return(df)
}

logs_ratio <- function(kmer_data){
  obs <- kmer_data[, "observed"]
  exp <- kmer_data[, "expected"]
  kmer_data[, "log_ratio"] <- log(obs / exp, base=10)
  return(kmer_data)
}
    
log_likelihood_score <- function(kmer_data){
  exp <- kmer_data[, "expected"]
  lgr <- kmer_data[, "log_ratio"]
  kmer_data[, "log_likeood"] <- exp * log(lgr)
  return(kmer_data)
}
    
poisson_pval <- function(kmer_counts, kmer_data, k){
  N <- sum(subset(kmer_counts$observed, nchar(kmer_counts$kmer) == 1))
  k_prob <- kmer_prob(kmer_counts, kme)
  obs <- kmer_data[, "observed"]
  kmer_data[, "pval_pois"] <- ppois(obs -1, 
                                    lambda=kmer_pr * N,
                                    lower.tail=F)
  return(kmer_data)
}
    
get_kmer_sd <- function(kmer_data){
  # calculate the data standard deviation
  kmer_data$sd <- sqrt(kmer_data$variance)
  return(kmer_data)
}
    
    
get_kmer_Zscores <- function(kmer_data){
  # Calculates the z scores from the kmer data to select 
  # over/under kmers
  obs <- kmer_data[, "observed"]
  exp <- kmer_data[, "expected"]
  sd <- kmer_data[, "sd"]
  kmer_data[, "Zscore"] <- (obs - exp) / sd
  return(kmer_data)
}
    
pval <- function(kmer_data){
  z <- kmer_data[, "Zscore"]
  kmer_data[, "Pval"] <- pnorm(-abs(z)) * 2
  return(kmer_data)
}
    
get_e_vals <- function(kmer_data, k){
  p <- kmer_data[, "Pval"]
  kmer_data[, "Eval"] <- (4^k) * p
  return(kmer_data)
}

# Selecting the final kmer data with statistics
exceptional_words <-function(filename, columns, eval=0.001, k){
  kmer_counts <- read.csv(filename, 
                          header=FALSE, 
                          sep = ",", 
                          col.names = c("kmer", 
                                        "observed", 
                                        "observed_freq"))
  km <- get_kmer_expected_var_HMO(kmer_counts, columns, k)
  km <- kmer_frequency(km, k, col="expected")
  km <- get_kmer_sd(km)
  km <- get_kmer_Zscores(km)
  km <- pval(km)
  km <- get_e_vals(km, k)
  words <- km[which(km$Eval < eval), ]
  words <-words[order(words$Zscore), ]
  words$rank <- rank(words$Zscore)
  return(words)
}
    