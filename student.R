## ##############################################################
##
## Draw a series of Student distributions with various degrees of
## freedom
##
## Author: Jacques van Helden
##
## Running this script requires to first run the script config.R
##
## Single-command execution:
## source(file.path(dir.R.files, 'student.R'))
##
## NOTE: THIS SCRIPT IS PROBABLY PARTLY REDUNDANT WITH
## golub/t-test_golub.R, BUT THERE ARE INTERESTING PLOTS HERE. TO
## CHECK.

x <- seq(from = -5, to=5, by=0.01)

in.colors <- TRUE ## Activate this to generate the color drawing for the book cover
freedom <- c(30,10,5,2,1)

if (in.colors) {
  palette <-   rainbow(10)
  bg <- '#000000'
  fg <- '#FFFFFF'
  file.prefix <- "student_distrib_cover"
} else {
  palette <- gray((0:16)/16)
  bg <- '#FFFFFF'
  fg <- '#000000'
  file.prefix <- "student_distrib"
}
lty <- c(1,6,5,4,2,3,1,6,5,4,2,3)

X11(width=7,height=5)
par(bg=bg)
par(fg=fg)

i <- 1
plot(x,dnorm(x),
     type="l",
     col=palette[i],
     lwd=2,
     lty=1,
     main="Student distributions",
     panel.first=grid(),
     xlab="x",
     ylab="density",
     col.axis=fg,
     col.lab=fg
     )
for (n in freedom) {
  i <- i+1
  lines(x,dt(x,n),
        type="l",
        col=palette[i],
        lwd=2,
        lty=lty[i]
        )
}

legend("topright",
       legend=c("normal",paste("student, df=",freedom)),
       col=palette[1:(length(freedom)+1)],
       lty=lty,
       lwd=2,
       bty="o",
       bg=bg)

setwd(dir.figures); export.plot(file.prefix=file.prefix, export.formats=export.formats.plots, width=7,height=5)
