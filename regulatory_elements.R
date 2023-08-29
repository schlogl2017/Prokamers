################################################################
#
# Analysis of regulatory elements
#

fig.dir <- '../figures/'
print.files <- T

#### load the data
data.type <- "regulatory_elements";
source("load_phosphate.R")



################################################################
#
# Plot some XY graphs
#
plot2var <- function (xcol=NA,ycol=NA,print.all=F) {
  if (is.na(xcol)) {xcol <- 1}
  if (is.na(ycol)) {ycol <- 1}
  par(font.lab=2)
  par(cex.lab=1.5)
  par(cex.main=1.5)
  if (print.all) {
    plot(data[,xcol],data[,ycol],col=palette.families["ALL"],type="p",pch=20,main=paste(data.title),xlab=xcol, ylab=ycol)
    points(ctl.fam.profiles[,xcol],ctl.fam.profiles[,ycol],col=palette.families["CTL"],pch=19)
  } else {
    plot(ctl.fam.profiles[,xcol],ctl.fam.profiles[,ycol],col=palette.families["CTL"],type="p",pch=20,main=paste(data.title),xlab=xcol, ylab=ycol)
  }
  points(met.fam.profiles[,xcol],met.fam.profiles[,ycol],col=palette.families["MET"],pch=19)
  points(pho.fam.profiles[,xcol],pho.fam.profiles[,ycol],col=palette.families["PHO"],pch=19)
  par(cex.main=1)
  par(cex.lab=1)
  par(font.lab=1)
}

if (print.files) {postscript(paste(fig.dir,data.type,'_selected_variable_pairs.ps',sep=''))}
par(mai=c(0.6,0.6,0.4,0.1))
par(mfrow=c(3,3))
plot2var("Pho4p.cacgtk.1","Met4p.1")
plot2var("Pho4p.cacgtk.1","Met31p.1")
plot2var("Met4p.1","Met31p.1")
plot2var("Pho4p.cacgtk.1","Pho4p.cacgtt.1")
plot2var("Pho4p.cacgtg.1","Pho4p.cacgtt.1")
plot2var("Pho4p.cacgtk.1","Pho4p.cacgtk.2")
plot2var("Pho4p.cacgtk.1","Pho4p.cacgtk.3")
plot2var("Met4p.1","Met4p.2")
plot2var("Met31p.1","Met31p.2")
par(mfrow=c(1,1))
par(mai=c(1,1,1,0.3))
if (print.files) {dev.off()}
