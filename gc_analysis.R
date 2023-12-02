library("ggplot2")
library("tidyverse")

filename <- "./Results/GC_sldw/Viruses/Unclassified_viruses/GCA_000838185.1.fna_nuc_100.txt"
col_names <- c("ID", 
               "start", 
               "end", 
               "AT", 
               "CG", 
               "A", 
               "C", 
               "G", 
               "T", 
               "N", 
               "other",
               "len_seq")
data = read.table(filename,
                 header=FALSE, 
                 sep = "\t", 
                 col.names = col_names)
data <- data[which(data$len_seq == 100), ]
data["diff"] <- data[["CG"]] - cg_mean

cg_total <- sum(data[["CG"]])
cg_mean <- mean(data[["CG"]]) # mean(data$CG)
cg_skew <- (data$C - data$G) / (data$C + data$G)
min(cg_skew)
gc_var <- log(sum(abs(data$diff)) / nrow(data))

ggplot(data, aes(start, diff, colour = "red")) + 
  geom_line()

ggplot(data, aes(start, cg_skew, colour = "red")) + 
  geom_line()

ggplot(data, aes(start, AT)) + 
  geom_line()

plot(data$start, data$CG, 
     type = "l",
     ylim = c(0, 1.0),
     col = "red") 
lines(data$start, 
      data$AT, 
      type = "l",
      col = "blue")


plot(data$start, data$A, type = "l", col = 1, ylim = c(0, 100))  # Plot with Base R
lines(data$start, data$C, type = "l", col = 2)
lines(data$start, data$G, type = "l", col = 3)
lines(data$start, data$T, type = "l", col = 4)
lines(data$start, data$N, type = "l", col = 5)
lines(data$start, data$other, type = "l", col = 6)

