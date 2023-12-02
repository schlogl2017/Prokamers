indir <- "Oligo_plots/Data"
pat <- paste0("_", 2, "k_", "chr", ".data.csv")
filename <- list.files(indir, pattern = pat, full.names = TRUE)
# read the data
oligo_data <-  read.csv(filename)
# separate the data
arc <- oligo_data[oligo_data$Domain == "Archaea", ]
bac <- oligo_data[oligo_data$Domain == "Bacteria", ]


plot(x = 1,
     type = "n",
     xlim = c(0, 16), 
     ylim = c(0, 1),
     pch = 16,
     xlab = "Height", 
     ylab = "Weight",
     main = "Adding points to a plot with points()")
# Add gridlines to a plot with grid()
plot(arc$kmer,
     as.numeric(arc$frequency),
     pch = 16,
     col = gray(.1, .2))
plot(bac$kmer,
     as.numeric(bac$frequency),
     pch = 16,
     col = blue(.1, .2))
# Add gridlines
grid()


