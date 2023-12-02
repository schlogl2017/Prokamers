library(dplyr)
library(knitr)
library(tidyverse)
library("PerformanceAnalytics")
library(reshape2)
library(ggplot2)
library(corrplot)
library(Hmisc)

## bom pro script
filename <- "Results/oligo_usage/oligo_usage_all_chr_kmer4.tb"
t<- read.table(filename, 
               sep = '\t', 
               header = T, 
               as.is = TRUE)

x <- t %>% select_if(is.numeric)
# data frame
corr_data <- cor(x)

# can get the correlation and its pvalues
corr_data1 <- rcorr(as.matrix(x))

# Use heatmap()
heatmap(x = corr_data, symm = TRUE)

melted_cormat <- melt(corr_data)
head(melted_cormat)
# Var1         Var2     value
# 1           Acidiphilium Acidiphilium 1.0000000
# 2 Acidipropionibacterium Acidiphilium 0.7335886
# 3      Acidithiobacillus Acidiphilium 0.4892160
# 4         Acidobacterium Acidiphilium 0.6162587
# 5           Acidothermus Acidiphilium 0.7865565
# 6             Acidovorax Acidiphilium 0.7705472

ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

## testes
# heat map 2 - same as above
heatmap(x = corr_data1$r, symm = TRUE)

## To compute the matrix of p-value, a custom R function is used :
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

p.mat <- cor.mtest(x)
head(p.mat[, 1:5])


melted_cormat <- melt(corr_data)
head(melted_cormat)
# Var1         Var2     value
# 1           Acidiphilium Acidiphilium 1.0000000
# 2 Acidipropionibacterium Acidiphilium 0.7335886
# 3      Acidithiobacillus Acidiphilium 0.4892160
# 4         Acidobacterium Acidiphilium 0.6162587
# 5           Acidothermus Acidiphilium 0.7865565
# 6             Acidovorax Acidiphilium 0.7705472

ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

## a function that formats a data frame correctly for 
# Hmisc::rcorr() and then turns each of the three 
# elements of the list ( r, n and P)
cors <- function(df) { 
  # turn all three matrices (r, n, and P into a data frame)
  M <- Hmisc::rcorr(as.matrix(df))
  # return the three data frames in a list return(Mdf)
  Mdf <- map(M, ~data.frame(.x))

cors(x) %>% first() %>% head() %>% kable()


library(heatmaply)
### Let's Plot
heatmaply_cor(x = cor(x),
              xlab = "Features",
              ylab = "Features",
              k_col = 2,
              k_row = 2)
              
              
              
              
              
              
              
df.scaled <- scale(x)              
dist.eucl <- dist(df.scaled, method = "euclidean")              
              
dist.cor <- get_dist(df.scaled, method = "pearson")
# Red: high similarity (ie: low dissimilarity) | Blue: low similarity
fviz_dist(dist.eucl)             
              
# Compute the dissimilarity matrix
# df = the standardized data
res.dist <- dist(x, method = "euclidean")              
res.hc <- hclust(d = res.dist, method = "complete")              
fviz_dend(res.hc, cex = 0.5)              
              
# Compute cophentic distance
res.coph <- cophenetic(res.hc)
# Correlation between cophenetic distance and
# the original distance
cor(res.dist, res.coph)              
              
res.hc1 <- hclust(d = res.dist, method = "average")              
fviz_dend(res.hc1, cex = 0.5)              

# Compute cophentic distance
res.coph1 <- cophenetic(res.hc1)
# Correlation between cophenetic distance and
# the original distance
cor(res.dist, res.coph1)              

# The correlation coefficient shows that using a different linkage method creates a tree
# that represents the original distances slightly better.

# Cut tree into 4 groups
grp <- cutree(res.hc, k = 8)
head(grp, n = 8)

# Number of members in each cluster
table(grp)

# Using the function fviz_cluster() [in factoextra], we can also visualize 
# the result in a scatter plot. Observations are represented by points in 
# the plot, using principal components. A frame is drawn around each cluster.
# Cut in 4 groups and color by groups
fviz_dend(res.hc, k = 8, # Cut in four groups
          cex = 0.5, # label size
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
)
options(ggrepel.max.overlaps = Inf)
fviz_cluster(list(data = x, cluster = grp),
             ellipse.type = "convex", # Concentration ellipse
             repel = TRUE, # Avoid label overplotting (slow)
             show.clust.cent = FALSE, ggtheme = theme_minimal())


library("cluster")

# Agglomerative Nesting (Hierarchical Clustering)
res.agnes <- agnes(x = x, # data matrix
                   stand = TRUE, # Standardize the data
                   metric = "euclidean", # metric for distance matrix
                   method = "ward" # Linkage method
)
# DIvisive ANAlysis Clustering
res.diana <- diana(x = x, # data matrix
                   stand = TRUE, # standardize the data
                   metric = "euclidean" # metric for distance matrix
)

fviz_dend(res.agnes, cex = 0.6, k=4)

library(dendextend)
# Compute distance matrix
res.dist <- dist(x, method = "euclidean")
# Compute 2 hierarchical clusterings
hc1 <- hclust(res.dist, method = "average")
hc2 <- hclust(res.dist, method = "ward.D2")
# Create two dendrograms
dend1 <- as.dendrogram (hc1)
dend2 <- as.dendrogram (hc2)
# Create a list to hold dendrograms
dend_list <- dendlist(dend1, dend2)

tanglegram(dend1, dend2)

tanglegram(dend1, dend2,
           highlight_distinct_edges = FALSE, # Turn-off dashed lines
           common_subtrees_color_lines = FALSE, # Turn-off line colors
           common_subtrees_color_branches = TRUE, # Color common branches
           main = paste("entanglement =", round(entanglement(dend_list), 2))
)
# Note that "unique" nodes, with a combination of labels/items not present in 
# the other tree, are highlighted with dashed lines.

# Cophenetic correlation matrix
cor.dendlist(dend_list, method = "cophenetic")

# Baker correlation matrix
cor.dendlist(dend_list, method = "baker")

# Create multiple dendrograms by chaining
dend1 <- x %>% dist %>% hclust("complete") %>% as.dendrogram
dend2 <- x %>% dist %>% hclust("single") %>% as.dendrogram
dend3 <- x %>% dist %>% hclust("average") %>% as.dendrogram
dend4 <- x %>% dist %>% hclust("centroid") %>% as.dendrogram
# Compute correlation matrix
dend_list <- dendlist("Complete" = dend1, "Single" = dend2,
                      "Average" = dend3, "Centroid" = dend4)
cors <- cor.dendlist(dend_list)
# Print correlation matrix
round(cors, 2)

# Visualize the correlation matrix using corrplot package
library(corrplot)
corrplot(cors, "pie", "lower")


# Compute distances and hierarchical clustering
dd <- dist(scale(x), method = "euclidean")
hc <- hclust(dd, method = "complete")

library(factoextra)
fviz_dend(hc, cex = 0.5)

fviz_dend(hc, cex = 0.5,
          main = "Dendrogram - complete",
          xlab = "Objects", ylab = "Distance", sub = "")

# To draw a horizontal dendrogram, type this:
fviz_dend(hc, cex = 0.5, horiz = TRUE)

fviz_dend(hc, k = 4, # Cut in four groups
          cex = 0.5, # label size
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE, # Add rectangle around groups
          rect_border = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          rect_fill = TRUE)

# To change the plot theme, use the argument ggtheme
fviz_dend(hc, k = 4,
          # Cut in four groups
          cex = 0.5,
          # label size
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, # color labels by groups
          ggtheme = theme_gray()
          # Change theme
)

# Allowed values for k_color include brewer palettes from RColorBrewer Package
fviz_dend(hc, cex = 0.5, k = 4, # Cut in four groups
          k_colors = "jco")

#If you want to draw a horizontal dendrogram with rectangle around clusters, use this:
fviz_dend(hc, k = 4, cex = 0.4, horiz = TRUE, k_colors = "jco",
            rect = TRUE, rect_border = "jco", rect_fill = TRUE)

# Additionally, you can plot a circular dendrogram using the option type = “circular”.
fviz_dend(hc, cex = 0.5, k = 4,
          k_colors = "jco", type = "circular")

# To plot a phylogenic-like tree, use type = “phylogenic” and 
# repel = TRUE (to avoid labels overplotting). This functionality 
# requires the R package igraph. Make sure that
# it’s installed before typing the following R code
require("igraph")
fviz_dend(hc, k = 4, k_colors = "jco",
          type = "phylogenic", repel = TRUE)

fviz_dend(hc, k = 4, # Cut in four groups
          k_colors = "jco",
          type = "phylogenic", repel = TRUE,
          phylo_layout = "layout.gem")

# for large data
# Zooming in the dendrogram
fviz_dend(hc, xlim = c(1, 10), ylim = c(1, 8))


# Create a plot of the whole dendrogram,
# and extract the dendrogram data
dend_plot <- fviz_dend(hc, k = 4, # Cut in four groups
                       cex = 0.5, # label size
                       k_colors = "jco"
)
dend_data <- attr(dend_plot, "dendrogram") # Extract dendrogram data
# Cut the dendrogram at height h = 10
dend_cuts <- cut(dend_data, h = 10)
# Visualize the truncated version containing
# two branches
fviz_dend(dend_cuts$upper)

# Plot the whole dendrogram
print(dend_plot)

# not working
# Plot subtree 1
fviz_dend(dend_cuts$lower[[1]], main = "Subtree 1")
# Plot subtree 2
fviz_dend(dend_cuts$lower[[2]], main = "Subtree 2")
fviz_dend(dend_cuts$lower[[2]], type = "circular")


# install.packages("ape")
library("ape")
# Default plot
hc2 <- hclust(dist(x), method = "complete")
plot(as.phylo(hc2), cex = 0.6, label.offset = 0.5)

library("ggdendro")
# Visualization using the default theme named theme_dendro()
ggdendrogram(hc)
# Rotate the plot and remove default theme
ggdendrogram(hc, rotate = TRUE, theme_dendro = FALSE)

## Extract dendrogram plot data
# Build dendrogram object from hclust results
dend <- as.dendrogram(hc)
# Extract the data (for rectangular lines)
# Type can be "rectangle" or "triangle"
dend_data <- dendro_data(dend, type = "rectangle")
# What contains dend_data
names(dend_data)
# Extract data for line segments
head(dend_data$segments)
# Extract data for labels
head(dend_data$labels)
# Plot line segments and add labels
p <- ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label),
            hjust = 1, angle = 90, size = 3)+
  ylim(-3, 15)
print(p)

dend2 <- as.dendrogram(hc)
plot(dend2)

# Create a dendrogram and plot it
dendx <- x[1:5,] %>%  scale %>% 
  dist %>% hclust %>% as.dendrogram
dendx %>% plot

# Get the labels of the tree
labels(dendx)
# Change the labels, and then plot:
dendx %>% set("labels", c("a", "b", "c", "d", "e")) %>% plot

# pvclust and dendextend
# The package dendextend can be used to enhance many packages including pvclust. 
# Recall that, pvclust is for calculating p-values for hierarchical clustering.
# pvclust can be used as follow:
library(pvclust)
#data(lung) # 916 genes for 73 subjects
#set.seed(1234)
result <- pvclust(x[1:100, 1:10], method.dist="cor", 
                  method.hclust="average", nboot=10)

# Default plot of the result
plot(result)
pvrect(result)

resultall <- pvclust(x, method.dist="cor", 
                  method.hclust="average", nboot=10)

# Default plot of the result
plot(resultall)
pvrect(resultall)
