################################################################
#
# Testing multivariate analysis with gene expression data
#
################################################################

#print.files <- F # save graphics in postscript files rather than displaying them on the screen

#### load the multivariate analysis libraries
#library(mva)
library(mda)
library(MASS)

################################################################
#### Read and prepare data for analysis
################################################################
source(file.path(dir.util, 'util_descr_stats.R'))
source(file.path(dir.util, 'util_chip_analysis.R'))

data.type <- "expression";
#data.type <- "regulatory_elements";
#data.type <- "random";
source(file.path(dir.R.files, "ogawa1998", "load_phosphate.R"))
source(file.path(dir.util, 'util_sda.R'))

################################################################
# Plot distributions of the variables
################################################################
x11(width=7,height=5)
par(mai=c(1,1,1,0.2))
plot.frequency.distrib (data,main=paste('Value distributions -', data.title))
setwd(dir.figures); export.plot(file.prefix=paste(sep='', dir.figures,data.type,'_variable_distributions'), export.formats=export.formats.plots, width=7,height=5)


################################################################
# Plot the first two columns 
################################################################
x11(width=6,height=6)
plot(data[,xcol],data[,ycol],col=palette.families["ALL"],pch=20,main=paste(data.title,' - two selected variables',sep=""),xlab=xcol, ylab=ycol, panel.first=grid())
points(ctl.fam.profiles[,xcol],ctl.fam.profiles[,ycol],col=palette.families["CTL"],pch=20)
#points(met.fam.profiles[,xcol],met.fam.profiles[,ycol],col=palette.families["MET"],pch=19)
points(pho.fam.profiles[,xcol],pho.fam.profiles[,ycol],col=palette.families["PHO"],pch=12)
legend("bottomright", ,legend=c("All genes", "Control", "Pho4p regulon"),
       col=palette.families[c("ALL", "CTL", "PHO")],
       pch=c(19,19,12),
       lwd=0,bg='white',bty='o')
##points(training[,1],training[,2],col=colors.training,pch=22)
setwd(dir.figures); export.plot(file.prefix=paste(data.type,'_var1_var2.ps',sep=''), export.formats=export.formats.plots, width=6,height=6)


################################################################
# Plot all pairs of coordinates
################################################################
#	if (print.files) {postscript(paste(dir.figures,data.type,'_variable_pairs.ps',sep=''))}
#	pairs(data)
#	cor(data,use="complete.obs")
#	## if (print.files) {dev.off()}

################################################################
# Draw profiles
################################################################
x11(width=10, height=9)
par(mfrow=c(2,2))
plot.profiles(data,row.names(pho.fam.profiles),main=paste(data.title, "- PHO genes"),col.profile=palette.families["PHO"])
plot.profiles(data,row.names(met.fam.profiles),main=paste(data.title, "- MET genes"),col.profiles=palette.families["MET"])
plot.profiles(data,row.names(ctl.fam.profiles),main=paste(data.title, "- CTL genes"),col.profiles=palette.families["CTL"])
plot.profiles(data,main=paste(data.title, "- ALL genes"),col.profiles=palette.families["ALL"])
par(mfrow=c(1,1))
setwd(dir.figures); export.plot(file.prefix=paste(data.type,'_profiles_per_group.ps',sep=''), export.formats=export.formats.plots, width=10,height=9)

################################################################
#### test some basic multivariate functions
################################################################

#### calculate the variance, covariance and correlation matrices between chips
var(data,use="complete.obs")
cov(data,use="complete.obs")
cor(data,use="complete.obs")

#### eigen values
eigen(var(data,use="complete.obs"))

################################################################
#### principal component analysis with prcomp
################################################################
##if (print.files) {postscript(paste(dir.figures,data.type,'_pca_two_component_plot.ps',sep=''))}
x11(width=6, height=6)
pc <- prcomp(data)
pc1 <- pc$x[,"PC1"]
pc2 <- pc$x[,"PC2"]
plot(pc1,pc2,
     main        = paste("PCA - ", data.title,sep=""),
     col         = palette.families["ALL"],
     xlab        = "PC1", 
     ylab        = "PC2",
     panel.first = c(grid(col=1),
       abline(h=0,col=1),
       abline(v=0,col=1)),
     pch=20)
points(pc$x[row.names(ctl.fam.profiles),"PC1"],pc$x[row.names(ctl.fam.profiles),"PC2"],
       col         = palette.families["CTL"],
       pch=19)
points(pc$x[row.names(met.fam.profiles),"PC1"],pc$x[row.names(met.fam.profiles),"PC2"],
       col         = palette.families["MET"],
       pch=19)
points(pc$x[row.names(pho.fam.profiles),"PC1"],pc$x[row.names(pho.fam.profiles),"PC2"],
       col         = palette.families["PHO"],
       pch=19)
legend("topright",c("PHO","MET","CTL"), col=palette.families,lwd=3,cex=0.8,bg='white',bty='o')
setwd(dir.figures); export.plot(file.prefix=paste("phosphate_",data.type,'_pca_two_component_plot',sep=''), export.formats=export.formats.plots, width=6,height=6,bg='white',bty='o')

################################################################
#### principal component analysis with princomp
################################################################

pc2 <- princomp(data,cor=TRUE)
##summary.princomp(pc2,loadings = TRUE, cutoff = 0.2, digits = 2)
summary(pc2,loadings = TRUE, cutoff = 0.2, digits = 2)
tpc2 <- princomp(training[,-length(training)],cor=TRUE)

#### plot component variances 
x11(width=7,height=5)
plot(pc2, main=paste("PCA component variances - ",data.title,sep=""))  
setwd(dir.figures); export.plot(file.prefix=paste("phosphate_",data.type,'_var_per_comp',sep=''), export.formats=export.formats.plots, width=7,height=5)

#### princomp with training data only
if (print.files) {postscript(paste(dir.figures,data.type,'_pca_biplot.ps',sep=""))}
par(mfrow=c(1,2))
par(cex=0.8)
biplot(tpc2,
       main=paste("PCA plot - training set - ",data.title, sep=""),
       col=c("#0000bb","#ff0000"))
par(cex=1)

#### plot objects along the two principal components
par(cex=0.8)
biplot(pc2,
       main=paste("PCA plot - whole data set - ",data.title, sep=""),
       col=c("#0000bb","#ff0000"))
par(cex=1)
par(mfrow=c(1,1))
## if (print.files) {dev.off()}

################################################################
#### hierarchical clustering
################################################################


#### cluster the columns (variables) on the whole data set
hc.var <- hclust(dist(t(data)), "ave")
plot(hc.var,
     main=paste('Hierarchical clustering of variables - ',
       data.title),
     hang=-1)

#### cluster the columns (variables) on the training set
training.hc.var <- hclust(dist(t(training[,-length(training)])), "ave")
plot(training.hc.var,
     main=paste('Hierarchical clustering of variables - ',
       data.title),
     hang=-1)

#### cluster the rows (objects)
if (print.files) {postscript(paste(dir.figures,data.type,'_hierarchical_clustering.ps',sep=""))}
hc.obj <- hclust(dist(training[,-length(training)]), "ave")
par(cex=0.5)
par(cex.main=2)
plot(hc.obj,
     main=paste('Hierarchical clustering of objects - ',
       data.title),
     hang=-1,
     col="#008800")
par(cex=1)

## if (print.files) {dev.off()}


################################################################
#### Flexible Discriminant Analysis
################################################################

#### apparently, fda requires at least 3 families
tfda <- fda(family ~ .,data=training)
plot(tfda)
confusion(tfda)
fda.pred <- predict(tfda,training)
tfda.validation <- data.frame(
                              family=training$family,
                              predict=fda.pred,
                              row.names=row.names(training)
                              )
confusion(fda.pred,training$family)

################################################################
#### Linear discriminant analysis
################################################################
tlda <- lda(family ~ .,training,na.action=na.omit)
tlda.loo <- lda(family ~ .,training,na.action=na.omit,CV=TRUE)

lda.compa <- data.frame(
                        row.names=row.names(training),
                        family=training$family,
                        loo=tlda.loo$class,
                        predict(tlda,training)
                        )
lda.loo.pred.pho <- row.names(lda.compa[lda.compa$loo == "PHO",])
lda.loo.pred.met <- row.names(lda.compa[lda.compa$loo == "MET",])
lda.loo.pred.ctl <- row.names(lda.compa[lda.compa$loo == "CTL",])

#### confusion matrices
confusion(lda.compa$class,lda.compa$family)
confusion(lda.compa$loo,lda.compa$family)

loo.failures <- lda.compa[(lda.compa$family != lda.compa$loo),]
loo.failures
failures <- lda.compa[(lda.compa$family != lda.compa$class),]
failures
training[row.names(failures),]


################################################################
# Plot profiles of misclassified units
if (print.files) {postscript(paste(dir.figures,data.type,'_LDA_loo_profiles_misclassified.ps',sep=""))}
par(mfrow=c(2,2))
failures.pho <- row.names(failures[failures$class=="PHO",])
if (length(failures.pho > 0)) {plot.profiles(data,failures.pho,main="Misclassifications - classified as PHO",colors.training[failures.pho])}
failures.met <- row.names(failures[failures$class=="MET",])
if (length(failures.met > 0)) {plot.profiles(data,failures.met,main="Misclassifications - classified as MET",colors.training[failures.met])}
failures.ctl <- row.names(failures[failures$class=="CTL",])
if (length(failures.ctl > 0)) {plot.profiles(data,failures.ctl,main="Misclassifications - classified as CTL",colors.training[failures.ctl])}
par(mfrow=c(1,1))
## if (print.files) {dev.off()}

################################################################
# Plot profiles of genes as a function of their predicted class

if (print.files) {postscript(paste(dir.figures,data.type,'_LDA_loo_profiles.ps',sep=""))}
par(mai=c(0.7,0.7,0.7,0.1))
par(mfrow=c(2,2))
if (length(lda.loo.pred.pho > 0)) {
  plot.profiles(data,lda.loo.pred.pho,
                main="LDA - Leave-one-out - Predicted as PHO",
                colors.training[lda.loo.pred.pho])
}
if (length(lda.loo.pred.met > 0)) {
  plot.profiles(data,lda.loo.pred.met,
                main="LDA - Leave-one-out - Predicted as MET",
                colors.training[lda.loo.pred.met])
}
if (length(lda.loo.pred.ctl > 0)) {
  plot.profiles(data,lda.loo.pred.ctl,
                main="LDA - Predicted as CTL",
                colors.training[lda.loo.pred.ctl])
}
par(mfrow=c(1,1))
## if (print.files) {dev.off()}


#### WORKS UNDER WINDOWS BUT NOT UNDER LINUX
#biplot(tlda, main=paste("linear discriminant analysis - ", data.title))

#### plot the calibration result
if (print.files) {postscript(paste(dir.figures,data.type,'_LDA_internal.ps',sep=""))}
par(mai=c(0.7,0.7,0.7,0.1))
par(font=2)
par(cex=0.7)
par(cex.main=2)
plot(lda.compa$x.LD1,lda.compa$x.LD2,
     col=colors.training,
     main=paste("LDA Calibration -", data.title),
     pch=substr(as.vector(lda.compa$class),1,1),
     xlab="Linear discriminant function 1",ylab="Linear discriminant function 2"
     )
## if (print.files) {dev.off()}

if (print.files) {postscript(paste(dir.figures,data.type,'_LDA_leave_one_out.ps',sep=""))}
plot(lda.compa$x.LD1,lda.compa$x.LD2,
     col=colors.training,
     main=paste("LDA calibration - Leave-one-out - ", data.title),
     pch=substr(as.vector(lda.compa$loo),1,1),
     xlab="Linear discriminant function 1",ylab="Linear discriminant function 2"
     )
par(cex=1)
par(cex.main=1)
par(font=1)
par(mai=c(1,1,1,1))
## if (print.files) {dev.off()}

#### define reasonable prior for whole-genome predictions
my.prior <- c(0.01,0.01,0.98)
names(my.prior) <- c("PHO","MET","CTL")
tlda$prior <- my.prior
lda.pred<- data.frame(predict(tlda,data))
lda.pred.pho <- row.names(lda.pred[lda.pred$class=="PHO",])
lda.pred.met <- row.names(lda.pred[lda.pred$class=="MET",])

#### predict in the whole gene set
if (print.files) {postscript(paste(dir.figures,data.type,'_LDA_prediction.ps',sep=""))}
plot(lda.pred$x.LD1,lda.pred$x.LD2,
     main=paste("Whole data set predictions - ", data.title),
     col="#999999",
     xlab="Linear discriminant function 1",ylab="Linear discriminant function 2"
     )
if (length(lda.pred.pho) > 0) {
  points(lda.pred[lda.pred.pho,"x.LD1"],lda.pred[lda.pred.pho,"x.LD2"],col=palette.families["PHO"],pch=19)
  text(lda.pred[lda.pred.pho,"x.LD1"],lda.pred[lda.pred.pho,"x.LD2"],labels=lda.pred.pho,pos=1,font=2,cex=0.8,col=palette.families["PHO"])
}
if (length(lda.pred.met) > 0) {
  points(lda.pred[lda.pred.met,"x.LD1"],lda.pred[lda.pred.met,"x.LD2"],col=palette.families["MET"],pch=19)
  text(lda.pred[lda.pred.met,"x.LD1"],lda.pred[lda.pred.met,"x.LD2"],labels=lda.pred.met,pos=1,font=2,cex=0.8,col=palette.families["MET"])
}

## if (print.files) {dev.off()}

#### count predictions by class
dim(lda.pred[lda.pred$class=="PHO",])
dim(lda.pred[lda.pred$class=="MET",])
dim(lda.pred[lda.pred$class=="CTL",])

################################################################
#### Quadratic discriminant analysis
################################################################
tqda <- qda(family ~ .,data=training,na.action=na.omit)
tqda.loo <- qda(family ~ .,training,na.action=na.omit,CV=TRUE)

qda.compa <- data.frame(
                        family=training$family,
                        loo=tqda.loo$class,
                        class=predict(tqda,training)$class,
                        predict(tqda,training)$posterior,
                        row.names=row.names(training)
                        )

#### confusion matrices
confusion(predict(tqda)$class,training$family)
confusion(tqda.loo$class,training$family)

loo.failures <- qda.compa[(qda.compa$family != qda.compa$loo),]
loo.failures
failures <- qda.compa[(qda.compa$family != qda.compa$class),]
failures
training[row.names(failures),]
#biplot(tqda$PHO,tqda$MET, main=paste("linear discriminant analysis - ", data.title))

#### plot the calibration result
if (print.files) {postscript(paste(dir.figures,data.type,'_QDA_internal_analysis_plot.ps',sep=""))}
pca.training <- prcomp(training[,-length(training)])
par(mai=c(0.7,0.7,0.7,0.1))
par(font=2)
par(cex=0.7)
par(cex.main=2)
par(cex.axis=1)
plot(pca.training$x[,1:2],
     col=colors.training,
     main=paste("QDA - internal analysis - ", data.title),
     pch=substr(as.vector(qda.compa$class),1,1),
     xlab="PC1",ylab="PC2"
     )
## if (print.files) {dev.off()}

if (print.files) {postscript(paste(dir.figures,data.type,'_QDA_leave_one_out_plot.ps',sep=""))}
plot(pca.training$x[,1:2],,
     col=colors.training,
     main=paste("QDA - Leave-one-out - ", data.title),
     pch=substr(as.vector(qda.compa$loo),1,1),
     xlab="PC1",ylab="PC2"
     )
par(cex=1)
par(cex.main=1)
par(cex.axis=1)
par(font=1)
## if (print.files) {dev.off()}

#### define reasonable prior for whole-genome predictions
my.prior <- c(0.01,0.01,0.98)
names(my.prior) <- c("PHO","MET","CTL")
tqda$prior <- my.prior

#### predict in the whole gene set
qda.pred <- data.frame(predict(tqda,data))
qda.pred.pho <- row.names(qda.pred[qda.pred$class=="PHO",])
qda.pred.met <- row.names(qda.pred[qda.pred$class=="MET",])
qda.pred.ctl <- row.names(qda.pred[qda.pred$class=="CTL",])

length(qda.pred.pho)
length(qda.pred.met)

#### plot the predictions on the PCA plane
if (print.files) {postscript(paste(dir.figures,data.type,'_QDA_prediction.ps',sep=""))}
pca.data <- prcomp(data)
par(mai=c(1,1,0.7,0.1))
plot(pca.data$x[,1:2],
     col=palette.families["ALL"],
     main=paste("QDA - Whole data set predictions - ", data.title),
     xlab="PC1",ylab="PC2"
     )

if (length(qda.pred.pho) > 0) {
  points(pca.data$x[qda.pred.pho,1:2],col=palette.families["PHO"],pch=19)
  text(pca.data$x[qda.pred.pho,1:2],labels=qda.pred.pho,pos=1,font=2,cex=0.8,col=palette.families["PHO"])
}
if (length(qda.pred.met) > 0) {
  points(pca.data$x[qda.pred.met,1:2],col=palette.families["MET"],pch=19)
  text(pca.data$x[qda.pred.met,1:2],labels=qda.pred.met,pos=1,font=2,cex=0.8,col=palette.families["MET"])
}
## if (print.files) {dev.off()}

################################################################
# Plot profiles of predicted PHO and MET genes
if (print.files) {postscript(paste(dir.figures,data.type,'_QDA_pred_profiles.ps',sep=""))}
par(mai=c(0.7,0.7,0.7,0.1))
par(mfrow=c(2,2))
if (length(qda.pred.pho > 0)) {
  plot.profiles(data,qda.pred.pho,
                main="QDA - Predicted as PHO",
                col=colors.all[qda.pred.pho])
}
if (length(qda.pred.met > 0)) {
  plot.profiles(data,qda.pred.met,
                main="QDA - Predicted as MET",
                col=colors.all[qda.pred.met])
}
if (length(qda.pred.ctl > 0)) {
  plot.profiles(data,qda.pred.ctl,
                main="QDA - Predicted as CTL",
                col=colors.all[qda.pred.ctl])
#	colors.training[qda.pred.ctl])
}
par(mfrow=c(1,1))
## if (print.files) {dev.off()}

#### count predictions by class
length(qda.pred.pho)
length(qda.pred.met)
length(qda.pred.ctl)

################################################################
#
# Variable ordering
#
source('sda.R') # load the library
x<- training[,-length(training)]
y <- training$family
variable.ordering(x,y,method="qda")
