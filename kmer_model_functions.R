# Functions to use in kmer analysis
# Author: Paulo Sérgio Schlögl
# last modified set, 2023
# first written Apr, 2022

# Reading the data and subseting the data
select_data <- function(filename, k){
  # select the kmer data of length k from a big csv
  # with two columns: kmers, counts and frequencies
  data <- read.csv(filename, 
                   header = TRUE,
                   col.names = c("kmer", "observed"),
                   colClasses=c("kmer"="character", 
                                "observed"="numeric"),
                   stringsAsFactors = FALSE)
  return(subset(data, nchar(kmer) == k))
}
  
select_data_counts <- function(counts, columns, k){
  df_obs <- subset(counts, 
                   nchar(kmer) == k, 
                   select = columns)
  return(df_obs)
}

# k3 <- select_data_counts(kmer_counts, columns = c("kmer", "observed"), k = 3)
# nrow(k3)

bases_freq_from_counts <- function(counts, 
                                   k=1, 
                                   cols = c("kmer", 
                                            "observed_freq")){
  # get the bases characters as A C G T
  nucs <- subset(counts,
                 nchar(kmer) == 1, 
                 select = cols[1])[[1]]
  # get the frequencies of A C G T from the data
  freq <- as.numeric(c(subset(counts, 
                              nchar(kmer) == 1, 
                              select = cols[2]))[[1]])
  # naming the vector data
  names(freq) <- nucs
  return(freq)
}

# bases_freq_from_counts(kmer_counts)

# Calculate the Model depending on kmer length (k)
# kmer lemgth <= 2
kmer_priori <- function(counts, kmer){
  n <- nchar(kmer)
  ## get the bases frequencies
  base_freqs <- bases_freq_from_counts(counts)
  ## calculate word probability
  pr <- vector()
  pr[kmer] <- 1
  # Select the letter at position i of word W
  for (i in 1:n) {
    nuc <- substr(kmer,i,i)
    # update word probability
    pr[kmer] <- pr[kmer] * base_freqs[nuc]
  }
  return(pr[kmer])
}

# test AATG 0.002901841
#   A            A           T           G
#  0.2186290 * 0.2186290 * 0.2168582 * 0.2799512 = 0.002901842
# kmer_priori(kmer_counts, "AATG") 0.002901842

calculate_kmer_probs <- function(counts,
                                 columns = c("kmer", "observed"),
                                 k){
  data <- select_data_counts(counts, columns = columns, k = k)
  n <- nrow(data)
  priori <- data.frame(kmer = character(),
                       prob = numeric(0))
  for(i in 1:n){
    kmer <- data$kmer[i]
    pri <- kmer_priori(counts, kmer)
    priori[i, ] <- c("kmer" = kmer, prob = pri)
  }
  colnames(priori) <- c("kmer", "priori_prob")
  df <- merge(data, priori, by = "kmer")
  df$priori_prob <- as.numeric(df$priori_prob) 
  return(df)
}

# calculate_kmer_probs(kmer_counts, 
#                      columns = c("kmer", "observed"), 
#                      k = 3)

expected_zero_order <- function(kmer, counts, seq_len){
  # A simple Bernoulli model
  n <- nchar(kmer)
  # subsetting the dataframe where A,C,G,T ave counts and sum the values
  # thats thats is the sequence length
  # nucs <- as.numeric(subset(counts$observed, nchar(counts$kmer) == 1))
  N <- seq_len - n + 1
  k_prob <- kmer_priori(counts, kmer)
  kmer_exp <- N * k_prob
  kmer_var <- kmer_exp * (1 - kmer_exp / N)
  return(list(kmer = kmer, expected = kmer_exp, variance = kmer_var))
}

#expected_zero_order("AATG", kmer_counts)

get_kmer_expected_var_ZOM <- function(counts, columns, seq_len, k){
  # Calculates the expected kmer values based in a HMO of k - 2
  df <- select_data_counts(counts, columns, k)
  n <- nrow(df)
  expected <- data.frame(kmer = character(), 
                         expected = numeric(0), 
                         variance = numeric(0))
  for(i in 1:n){
    expected[i, ] <- expected_zero_order(df$kmer[i], 
                                         counts, seq_len)
  }
  merged_df <- merge(df, expected, by = "kmer")
  return(merged_df)
}

# get_kmer_expected_var_ZOM(kmer_counts, columns = c("kmer", "observed"), k = 3)


myvar <- function(kmer, counts, data){
  # it results in a higher variantion than
  # kmer_exp * (m - p) * (m - s) / (m^2)
  N <- sum(as.numeric(subset(counts$observed, 
                             nchar(counts$kmer) == 1)))
  ke <- data[which(data$kmer == kmer), "expected"]
  v <- ke * (1 - ke / N)
  return(v)
 }


# Higher Markove Orders
get_kmer_expected_var_HMO_new <- function(kmer_data, kmer_counts){
  # Calculates the expected kmer values based in a HMO of k - 2
  n <- nrow(kmer_data)
  kmers <- kmer_data$kmer
  expected <- data.frame(kmer = character(), variance = numeric(0))
  for(i in 1:n){
    kmer <- kmers[i]
    v <- expected_kmer_HMO_var(kmer, kmer_data, kmer_counts)
    expected[i, ] <- c("kmer" = kmer, variance = v)
  }
  merged_df <- merge(kmer_data, expected, by = "kmer")
  return(merged_df)
}

# expected_high_markov_var("TTG", kmer_counts)

    
# get_kmer_expected_var_HMO <- function(kmer_data, kmer_counts, columns, k){
#   # Calculates the expected kmer values based in a HMO of k - 2
#   df <- select_data_counts(kmer_counts, columns, k)
#   n <- nrow(df)
#   expected <- data.frame(kmer = character(), 
#                          expected = numeric(0), 
#                          variance = numeric(0))
#   for(i in 1:n){
#     expected[i, ] <- expected_kmer_HMO_var(df$kmer[i], kmer_data, kmer_counts)
#   }
#   merged_df <- merge(df, expected, by = "kmer")
#   return(merged_df)
# }

# get two data frames and kmer to calculate the variantio
variance_high_markov <- function(kmer, kmer_data, kmer_counts){
  kmer_exp <- numeric(0)
  n <- nchar(kmer)
  pref <- substr(kmer, 1, n - 1)
  p <- as.double(kmer_counts[which(kmer_counts$kmer == pref),
                             'observed'])
  suf <- substr(kmer, 2, n)
  s <- as.double(kmer_counts[which(kmer_counts$kmer == suf),
                             'observed'])
  mid <- substr(kmer, 2, n - 1)
  m <- as.double(kmer_counts[which(kmer_counts$kmer == mid),
                             'observed'])
  exp <- kmer_data[which(kmer_data$kmer == kmer), "expected"]
  kmer_exp <- exp
  kmer_var <- exp * (m - p) * (m - s) / (m^2)
  return(list(kmer = kmer, expected = kmer_exp, variance = kmer_var))
}


get_kmer_expected_var <- function(kmer_data, kmer_count){
  kmer_list <- kmer_data$kmer
  n <- length(kmer_list)
  expected <- data.frame(kmer = character(), expected = numeric(0), variance = numeric(0))
  for(i in 1:n){
    kmer <- kmer_list[i]
    exp_var <- variance_high_markov(kmer, kmer_data, kmer_counts)
    expected[i, ] <- exp_var
  }
  merged_df <- merge(kmer_data, expected, by = "kmer")
  return(merged_df)
}


get_kmer_expected_var_HMO_new <- function(kmer_data, kmer_counts){
  # Calculates the expected kmer values based in a HMO of k - 2
  n <- nrow(kmer_data)
  kmers <- kmer_data$kmer
  expected <- data.frame(kmer = character(), variance = numeric(0))
  for(i in 1:n){
    kmer <- kmers[i]
    v <- expected_kmer_HMO_var(kmer, kmer_data, kmer_counts)
    expected[i, ] <- c("kmer" = kmer, variance = v)
  }
  merged_df <- merge(kmer_data, expected, by = "kmer")
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

logs_ratio <- function(kmer_data){
  obs <- kmer_data[, "observed"]
  exp <- kmer_data[, "expected"]
  kmer_data[, "log_ratio"] <- log(obs / exp, base=10)
  return(kmer_data)
}

# data is already with kmer expected calculated
# data <- logs_ratio(kmer_data)
    
# logs_ratio <- function(kmer_data){
#   obs <- kmer_data[, "observed"]
#   exp <- kmer_data[, "expected"]
#   kmer_data[, "log_ratio"] <- log(obs / exp, base=10)
#   return(kmer_data)
# }

# data <-logs_ratio(data)

log_likelihood_score <- function(kmer_data){
  exp <- kmer_data[, "expected"]
  lgr <- kmer_data[, "log_ratio"]
  kmer_data[, "log_likelihood"] <- exp * log(abs(lgr))
  return(kmer_data)
}


poisson_pval <- function(counts, kmer_data, columns, k){
  N <- sum(subset(counts$observed, 
                  nchar(counts$kmer) == 1))
  k_probs <- kmer_data$priori_prob
  obs <- kmer_data$observed
  kmer_data[, "pval_pois"] <- ppois(obs - 1, 
                                    lambda = k_probs * N,
                                    lower.tail = FALSE)
  return(kmer_data)
}


variance <- function(data, N){
  # N = length sequence
  n <- nrow(data)
  for (i in 1:n){
    e <- data[which(data$kmer == data$kmer[i]), "expected"]
    var <- e * (1 - e/N)
    data[i, "variance"] <- e * (1 - e/N)
  }
  return(data)
}

myvar2 <- function(pref, suf, mid, expec){
  v <- expec * (mid - pref) * (mid - suf) / mid^2
  return(v)
}

get_kmer_sd <- function(kmer_data){
  # calculate the data standard deviation
  kmer_data$sd <- sqrt(as.numeric(kmer_data$variance))
  return(kmer_data)
}
    
    
get_kmer_Zscores <- function(kmer_data){
  # Calculates the z scores from the kmer data to select 
  # over/under kmers
  obs <- as.numeric(kmer_data[, "observed"])
  exp <- as.numeric(kmer_data[, "expected"])
  sd <- as.numeric(kmer_data[, "sd"])
  kmer_data[, "Zscore"] <- as.double((obs - exp) / sd)
  return(kmer_data)
}
    
pval <- function(kmer_data){
  z <- kmer_data[, "Zscore"]
  kmer_data[, "Pval"] <- as.double(pnorm(-abs(z)) * 2)
  return(kmer_data)
}
    
get_e_vals <- function(kmer_data, k){
  p <- kmer_data[, "Pval"]
  kmer_data[, "Eval"] <- (4^k) * p
  return(kmer_data)
}

# Selecting the final kmer data with statistics
exceptional_words_ <- function(filename1,
                              filename2,
                              eval = 0.001, 
                              k){
  kmer_counts <- read.csv(filename, 
                   header = FALSE, 
                   sep = ",", 
                   stringsAsFactors = FALSE)
  data <- read.csv(filename, 
                   header = FALSE, 
                   sep = ",", 
                   stringsAsFactors = FALSE)
  
  data <- get_kmer_expected_var_HMO_new(data, kmer_counts)
  data <- get_kmer_sd(data)
  data <- get_kmer_Zscores(data)
  data <- pval(data)
  data <- get_e_vals(data, k)
  words <- data[which(data$Eval < eval), ]
  words <- words[order(words$Zscore), ]
  words$rank <- rank(words$Zscore)
  return(words)
}

exceptional_words_in_genomes <-function(filename, 
                                        columns = c("kmer",
                                                    "observed",
                                                    "expected"), 
                                        eval = 0.001, 
                                        k){
  data <- read.csv(filename, 
                   header = TRUE, 
                   stringsAsFactors = FALSE)
  folder_spl <- strsplit(filename, split = "/")[[1]]
  name_spl <- strsplit(filename, split = "/")[[1]][7]
  name <- substr(name, start = 1, stop = 15)
  data <-variance(data, N)
  data <- get_kmer_sd(data)
  data <- get_kmer_Zscores(data)
  data <- pval(data)
  data <- get_e_vals(data, k)
  words <- data[which(data$Eval < eval), ]
  words <- words[order(words$Zscore), ]
  words$rank <- rank(words$Zscore)
  outfolder <- paste0("Results/Final_Data/",
                      folder_spl[3], 
                      "/",
                      folder_spl[4],
                      "/",
                      folder_spl[5] ,
                      "/", 
                      folder_spl[6], 
                      "/")
  path <- file.path(outfolder)
  if(!dir.exists(path)){
    dir.create(path, recursive = TRUE)
  }
  filename <- paste0(name, "_k", k, ".tsv")
  outname <- file.path(path, filename)
  print(outname)
  write.table(words,
  file = outname, sep = "\t",
            row.names = FALSE,
            quote = FALSE)
}    

dinuc_ab <- function(counts){
  di_ab <- vector()
  # a vector with bases frequencies
  bases <- bases_freq_from_counts(counts, 
                                     k=1, 
                                     cols = c("kmer", 
                                              "observed_freq"))
  # a dataframe with dinucleotides frequencies
  dinucs <- select_data_counts(counts, 
                               columns = c("kmer",
                                           "observed_freq"), 
                               k=2)
  # get the names to the vector with dinucleotides abundance
  dinames <- names(dinucs)
  n <- nrow(dinucs)
  for (i in 1:n) {
    # split the dinucleotides in the bases
    base_spl <- strsplit(dinucs$kmer[i], "")[[1]]
    kmer <- dinucs[dinucs$kmer[i], "observed_freq"]
    ab <- kmer / (bases[base_spl[1]] * bases[base_spl[2]])
    di_ab <- c(di_ab[i], ab)
  }
  names(di_ab) <- dinames
  return(di_ab)
}


exceptional_words_new <- function(filename1,
                                  filename2,
                                  eval = 0.001, 
                                  k){
  name <- tail(strsplit(filename1, split="/")[[1]], n=1)
  name <- substr(name, start = 1, stop = 15)
  cat("Load data in...")
  kmer_counts <- read.csv(filename1, 
                          header = TRUE, 
                          sep = ",", 
                          stringsAsFactors = FALSE)
  kdata <- read.csv(filename2, 
                   header = TRUE, 
                   sep = ",", 
                   stringsAsFactors = FALSE)
  cat("Data loaded!!! Done", "\n")
  cat("\n", "Calculating the exceptional words from Genome ID: ", name, "\n")
  kmer_data <- get_kmer_expected_var_HMO_new(kdata, kmer_counts)
  data <- get_kmer_sd(kmer_data)
  data <- get_kmer_Zscores(data)
  data <- pval(data)
  data <- get_e_vals(data, k)
  words <- data[which(data$Eval < eval), ]
  words <- words[order(words$Zscore), ]
  words$rank <- rank(words$Zscore)
  return(words)
}
