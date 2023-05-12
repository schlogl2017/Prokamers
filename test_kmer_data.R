select_data <- function(filename, k){
  # select the kmer data of length k from a big csv
  # with two columns: kmers, counts and frequencies
  data <- read.csv(filename, 
                   header=TRUE, 
                   stringsAsFactors=FALSE)
  return(subset(data, nchar(kmer) == k))
}

select_data_counts <- function(counts, k, columns){
  df_obs <- subset(counts, nchar(kmer) == k, select=columns)
  return(df_obs)
}

expected_zero_order <- function(kmer, kmer_counts){
  # A simple Bernulli model
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


kmer_prob <- function(kmer_counts, kmer, columns=c("kmer", "observed_freq")){
  splt_kmer <- unlist(strsplit(kmer, ""))
  bases <- select_data_counts(kmer_counts, 1, columns)
  kmer_prob <- 1.0
  for(base in splt_kmer){
    val <- subset(bases$observed_freq, bases$kmer == base)
    kmer_prob = kmer_prob * val
  }
  return(kmer_prob)
}

get_kmer_expected_var_ZOM <- function(kmer_counts, columns, k){
  # Calculates the expected kmer values based in a HMO of k - 2
  df <- df <- select_data_counts(kmer_counts, k, columns)
  n <- nrow(df)
  expected <- data.frame(kmer = character(), expected = numeric(0), variance = numeric(0))
  for(i in 1:n){
    expected[i, ] <- expected_zero_order(df$kmer[i], kmer_counts)
  }
  merged_df <- merge(df, expected, by = "kmer")
  return(merged_df)
}

expected_high_markov_var <- function(kmer, kmer_counts){
  # Function to calculate the expected values of counted kmers
  # the model is a higher order markov k - 2
  n <- nchar(kmer)
  seq_len <- sum(subset(kmer_counts$observed, nchar(kmer_counts$kmer) == 1))
  N <- seq_len - n + 1
  pref <- substr(kmer, 1, n - 1)
  p <- as.double(kmer_counts[kmer_counts$kmer == pref,
                             'observed'])
  suf <- substr(kmer, 2, n)
  s <- as.double(kmer_counts[kmer_counts$kmer == suf,
                             'observed'])
  mid <- substr(kmer, 2, n - 1)
  m <- as.double(kmer_counts[kmer_counts$kmer == mid,
                             'observed'])
  kmer_exp <- (p * s) / m
  kmer_var <- kmer_exp * (1 - kmer_exp / N)
  return(list(kmer = kmer, expected = kmer_exp, variance = kmer_var))
}

get_kmer_expected_var_HMO <- function(kmer_counts, columns, k){
  # Calculates the expected kmer values based in a HMO of k - 2
  df <- select_data_counts(kmer_counts, k, columns)
  n <- nrow(df)
  expected <- data.frame(kmer = character(), expected = numeric(0), variance = numeric(0))
  for(i in 1:n){
    expected[i, ] <- expected_high_markov_var(df$kmer[i], kmer_counts)
  }
  merged_df <- merge(df, expected, by = "kmer")
  return(merged_df)
}



filename <- "Jelly/Bacteria/Pseudomonadota/Acidithiobacillia/CHR/GCA_000559045.1_k1_10.csv"
k <- 6

# Data 
kmer_counts <- read.csv(filename, 
                        header=TRUE, 
                        stringsAsFactors=FALSE)
data <- select_data(filename, k)

# Analysis
kmer <- "TACCGG"

# Bernulli/ZOM
obs <- data[data$kmer == kmer, "observed"]
exp <- expected_zero_order(kmer, kmer_counts)[[2]]

# higher order
exp_h <- expected_high_markov_var(kmer, kmer_counts)[[2]]

obs/exp
obs/exp_h

hmo <- get_kmer_expected_var_HMO(kmer_counts,
                          k, 
                          columns=c("kmer", "observed"))

zom <- get_kmer_expected_var_ZOM(kmer_counts, 
                                 columns=c("kmer", "observed", "observed_freq"), 
                                 k)

splt_kmer <- unlist(strsplit(kmer, ""))
bases <- select_data_counts(kmer_counts, 1, c("kmer", "observed_freq"))
bases
kmer_prob <- 1.0
for(base in splt_kmer){
  val <- subset(bases$observed_freq, bases$kmer == base)
  kmer_prob = kmer_prob * val
}
print(kmer_prob)
0.2355125*0.2335181*0.2619436*0.2619436*0.2690258*0.2690258

