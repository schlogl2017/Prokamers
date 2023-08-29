## ##############################################################
## Compare various clustering methods on the data set from Golub, 1999

library(e1071) ## For matchCLasses

cluster.analysis <- function (x,data.name,dir.results, dir.figures) {
  
  dir.create(dir.results, recursive=T, showWarnings=F)
  dir.create(dir.figures, recursive=T, showWarnings=F)
  ## ##############################################################
  ## Calculate dissimilarity matrix
  
  ## Dimensions per row
  dim.per.row <- apply(!is.na(x),1,sum) ## Dimensions per row
  dim.per.pair <- as.matrix(dim.per.row) %*% t(as.matrix(dim.per.row))
  
  ## Euclidian distance
  gene.dist.eu <- dist(x, method="euclidian")
  gene.dist.epd <- gene.dist.eu/as.dist(dim.per.pair) ## Euclidian distance per dimension
  sample.dist.eu <- dist(t(x), method="euclidian")
  sample.dist.epd <- sample.dist.eu/as.dist(ncol(x)) ## Euclidian distance per dimension
  
  ## dissimilarity based on Pearson's correlation coefficient
  gene.dist.cor <- as.dist(1-cor(t(x),use="pairwise.complete.obs"))
  sample.dist.cor <- as.dist(1-cor(x,use="pairwise.complete.obs"))
  
  ## average dot product
  ## Warning: NA values are not treated properly !!!!!!!!!!!!
  gene.dist.dp <- as.matrix(x) %*% t(as.matrix(x))
  gene.dist.dp <- max(gene.dist.dp,na.rm=T) -gene.dist.dp ## Transform the similarity into a dissimilarity
  gene.dist.dp <- as.dist(gene.dist.dp)
  sample.dist.dp <- t(as.matrix(x)) %*% as.matrix(x)
  sample.dist.dp <- max(sample.dist.dp,na.rm=T) - sample.dist.dp ## Transform the similarity into a dissimilarity
  sample.dist.dp <- as.dist(sample.dist.dp)
  
  ## ##############################################################
  ## Compare the metrics
  
  ## Euclidian versus correlation
  x11(width=8,height=8)
  plot(gene.dist.epd,gene.dist.cor,
       main='metrics comparisons',
       xlab='Euclidian distance per dimension (d/p)',
       ylab='1 -correlation',
       panel.first=grid(col="black")
       )
  setwd(dir.figures); export.plot(file.prefix=paste(data.name, "eucl_vs_corr",sep="_"), export.formats=export.formats.plots, width=8,height=8)
  
  ## Euclidian versus dot product
  x11(width=8,height=8)
  plot(gene.dist.eu,
       gene.dist.dp,
       main='metrics comparisons',
       xlab='Euclidian distance',
       ylab='max(dot product) -(dot product)',
       panel.first=grid(col="black")
       )
  setwd(dir.figures); export.plot(file.prefix=paste(data.name, "eucl_vs_dotp", sep="_"), export.formats=export.formats.plots, width=8,height=8)
  
  ## correlation versus dot product
x11(width=8,height=8)
  plot(gene.dist.cor,gene.dist.dp,
       main='metrics comparisons',
       xlab='1 - (correlation)',
       ylab='max(dot product) -(dot product)',
       panel.first=grid(col="black")
       )
  setwd(dir.figures); export.plot(file.prefix=paste(data.name, "cor_vs_dotp", sep="_"), export.formats=export.formats.plots, width=8,height=8)
  
  ## ##############################################################
  ## Calculate a tree on the basis of the matrix
  
  ## plot parameters
  par(cex.main=2)
  par(cex=0.65)

################################################################
  ## effect of the linkage method

  ## With the Euclidian distance
  x11(width=16,height=12)
  par(mfrow=c(4,1))
  for (m in c("single","average","complete","ward")) {
    tree <- hclust(gene.dist.eu,method=m)
    plot(tree,main=paste(data.name, ";", ,m,"linkage", "; Euclidian distance"), xlab="")
  }
  par(mfrow=c(1,1))
  setwd(dir.figures); export.plot(file.prefix=paste(data.name, "hierarchical_eucl_linkage_method", sep="_"), export.formats=export.formats.plots, width=16,height=12)

  ## With the dot product
  x11(width=16,height=12)
  par(mfrow=c(4,1))
  for (m in c("single","average","complete","ward")) {
    tree <- hclust(gene.dist.dp,method=m)
    plot(tree,main=paste(data.name, ";",m,"linkage", "; Dot product"), xlab="")
  }
  par(mfrow=c(1,1))
  setwd(dir.figures); export.plot(file.prefix=paste(data.name, "hierarchical_dotp_linkage_method",sep="_"), export.formats=export.formats.plots, width=16,height=12)

################################################################
  ## effect of the dissimilarity metric
  m <- "complete"
  for (m in c('single','average','complete', 'ward')) {
    x11(width=16,height=12)
    par(mfrow=c(3,1))

    ## Euclidian distance
    tree.eu <- hclust(gene.dist.eu,method=m)
    plot(tree.eu,main=paste(data.name, "; ",m,"linkage", "; Euclidian distance"), xlab="")
    
    ## Dot product
    tree.dp <- hclust(gene.dist.dp,method=m)
    plot(tree.dp,main=paste(data.name, ";",m,"linkage", "; Dot product"), xlab="")

    ## Coefficient of correlation
    tree.cor <- hclust(gene.dist.cor,method=m)
    plot(tree.cor,main=paste(data.name, ";",m,"linkage", "; Correlation"), xlab="")
    
    setwd(dir.figures); export.plot(file.prefix=paste(data.name, m ,
                                      "hierarchical_distance_metrics",
                                      sep='_'),
                                    export.formats=export.formats.plots, width=16,height=12)
    par(mfrow=c(1,1)) 
  }


  max.level <- max(abs(range(x)))
  breaks <- 2*max.level*(0:255)/255-max.level


################################################################
  ## Generate heat maps with the different metrics and agglomeration
  ## rules
  for (m in c('single','average','complete', 'ward')) {
    for (dist in c("eu", "cor", "dp")) {
      gene.dist.matrix <- get (paste("gene.dist.",dist,sep=""))
      gene.tree <- hclust(gene.dist.matrix,method=m)
      sample.dist.matrix <- get (paste("sample.dist.",dist,sep=""))
      sample.tree <- hclust(sample.dist.matrix,method=m)

      ## Plot heat map with Euclidian distance
      x11(width=8,height=9)
      heatmap(as.matrix(x), scale="none", col=green.to.red(),breaks=breaks,
              Rowv=as.dendrogram(gene.tree),
              Colv=as.dendrogram(sample.tree),
              main=paste(data.name, ";",dist,"distance;", m, "linkage", sep=" "))
      setwd(dir.figures);export.plot(file=paste(data.name, dist,m, "heatmap", sep='_'),export.format="pdf",width=8,height=9)
    }
  }


  ## ##############################################################
  ## Cut and prune the tree

  x11(width=16,height=8)
  par(cex.main=2)
  par(cex=0.65)
  k <- 8

  par(mfrow=c(2,1))
  gene.tree.eu.ward <- hclust(gene.dist.eu,method="ward")
  plot(gene.tree.eu.ward, main=paste(data.name, "; Euclidian distance; Ward linkage"))
  gene.clusters.hier.eu.ward <- cutree(gene.tree.eu.ward,k=k)
  table(gene.clusters.hier.eu.ward)

  par(cex=1.3)
  par(cex=1)
  library(maptree)
  pr <- prune.clust(gene.tree.eu.ward,k=k)
  plot(pr, main=paste("pruned tree, k=", k),xlab="")
  par(mfrow=c(1,1))

  setwd(dir.figures); export.plot(file.prefix=paste(data.name, "clustering_hierarchical_pruning",sep="_"), export.formats=export.formats.plots, width=16,height=8)

  ## restore plot parameters
  par(cex=1)
  par(cex.main=1)

  ## ##############################################################
  ## K-means clustering
  x11(width=16,height=8)
  clusters.k <- kmean.profiles(x, k)
  setwd(dir.figures); export.plot(file.prefix=paste(data.name, "_clustering_kmean_profiles_k", k, sep=""), export.formats=export.formats.plots, width=16,height=8)

  ## ##############################################################
  ## Compare results between different clustering results
  ## hierarchical vs k-means (both with Euclidian distance

  ## Contingency table
  tab <- table(gene.clusters.hier.eu.ward,clusters.k$cluster)
  tab <- as.data.frame(tab[,])
  names(tab) <- paste('k', names(tab),sep='')
  row.names(tab) <- paste('h', row.names(tab),sep='')
  print(tab)
  setwd(dir.results); write.table(tab, file.path(dir.results,paste(data.name, 'hierarchical_vs_kmeans.tab', sep="_")), sep='\t',quote=F)

  ## Find optimal match between the two classifications
  class.match <- matchClasses(as.matrix(tab),method="exact")
  tab.matched <- as.data.frame(tab[,class.match])
  print(tab.matched)
  setwd(dir.results); write.table(tab.matched, file.path(dir.results,paste(data.name, 'hierarchical_vs_kmeans_matched.tab', sep="_")), sep='\t',quote=F)

  ## Calculate statistics on the matched classes
  ## This function permutes the columns and rows of the contingency
  ## table in order to maximize the diagonal
  match.comp <- compareMatchedClasses(gene.clusters.hier.eu.ward,clusters.k$cluster, method="exact")
  hit.rate <- match.comp$diag
  kappa <- match.comp$kappa
  setwd(dir.results); sink(paste(data.name, 'hierarchical_vs_kmeans.txt', sep="_"))
  print("confusion table")
  print(tab)
  print(class.match)
  print("matched classes")
  print(tab.matched)
  print(c("hit.rate", hit.rate))
  print(c("kappa", kappa))
  sink()

  ## ##############################################################
  ## Plot K-means clustering result with expression profiles

  x11(width=12, height=8)
  k <- 11
  km <- kmean.profiles(x,k=k)
  setwd(dir.figures); export.plot(file.prefix=paste(data.name, "kmean_profiles", sep="_"), export.formats=export.formats.plots, width=12,height=8)


################################################################
                                        # SVD visualisation with coloring from kmeans clusters
  palette.rep	<- rep(palette,length.out=k)

  x11(width=10,height=10)
  plot.svd(x,col=palette.rep[km$cluster],main=paste(data.name, " - SVD plot"),cex.main=2)
  setwd(dir.figures); export.plot(file.prefix=paste(data.name,"SVD",sep="_"), export.formats=export.formats.plots, width=10,height=10)

#### Transformation in principal component space 
  x11(width=10,height=10)
  pcomp	 <- prcomp(x)
  plot(pcomp$x[,"PC1"],pcomp$x[,"PC2"],
       col = palette.rep[km$cluster],
       xlab = "PC1", 
       ylab = "PC2",
       pch=19,
       panel.first = c(grid(col=1),
         abline(h=0,col=1),
         abline(v=0,col=1)),
       main=paste( data.name, "; PCA plot"),
       cex.main=2
       )
  setwd(dir.figures); export.plot(file.prefix=paste(data.name,"PCA",sep="_"), export.formats=export.formats.plots, width=10,height=10)


################################################################
  ## Compare the results of two runs of the K-means clustering

  k <- 20
  ## First run of K-means clustering
  clusters.k.1 <- kmeans (x, k)
  table(clusters.k.1$cluster)

  ## Second run of K-means clustering
  clusters.k.2 <- kmeans (x, k)
  table(clusters.k.2$cluster)

  table(clusters.k.1$cluster,clusters.k.2$cluster)

  ## Compare hierarchical and k-means clustersing
  (contingency <- table(clusters.k.1$cluster, clusters.k.2$cluster))
  (class.match <- matchClasses(as.matrix(contingency),method="exact"))
  (confusion.table <- as.data.frame(contingency[,class.match]))
  (match.comp <- compareMatchedClasses(clusters.k.1$cluster, clusters.k.2$cluster, method="exact"))

################################################################
  ## Compare k-means clustering results with selected genes or all
  ## genes, respectively

  k.all <- 21
  clusters.k.all <- kmeans(na.omit(x),k.all,iter=100)

  k.select <- 20
  clusters.k.select <- kmeans(x,k.select)

  selected <- rep(F,nrow(x))
  selected[match(row.names(x),row.names(x))] <- T
  not.selected <- !selected

  clusters.selected <- rep("not.sel", length(clusters.k.all$cluster))
  names(clusters.selected) <- names(clusters.k.all$cluster)
  length(clusters.selected)
  clusters.selected[names(clusters.k.select$cluster)] <- clusters.k.select$cluster
  table(clusters.k.all$cluster,clusters.selected)

  k.all.2 <- kmean.profiles(na.omit(x), k.all)

  ## Hierarchical clustering
  tree.eucl.complete <- hclust(e, method='complete')
  plot(tree.eucl.complete)

  ## Prune the tree and diplsay the pruned tree
  pr <- prune(tree.eucl.complete,k=k)
  pr <- prune.clust(tree.eucl.complete,k=k)
  plot(pr, main=paste("pruned tree, k=", k),xlab="")

  ## Cut the tree
  gene.clusters.hier.eucl.complete <- cutree(tree.eucl.complete, k=k)

  contingency <- table(gene.clusters.hier.eucl.complete, clusters.k.select$cluster)
  class.match <- matchClasses(as.matrix(contingency),method="exact")
  confusion.table <- as.data.frame(contingency[,class.match])

  ## Compare hierarchical and k-means clustersing
  match.comp <- compareMatchedClasses(gene.clusters.hier.eucl.complete, clusters.k.select$cluster, method="exact")
  (match.comp$diag)

}
