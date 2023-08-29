#' @title Generate a sigmoid curve with a Hill function
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Generate a sigmoid curve with a Hill function
#'
#' @param x.min=0    Min value for x
#' @param x.max=1    Max value for x
#' @param s=0.4      Threshold value
#' @param lambda=10  Exponent of the hill function
#' @param n=101 ## Number of points 
#' 
#' @details
#' First version: ?
#' Last modification: 2015-10. 
#'
#' @examples
#' ha <- hill.distrib(x.min=0, x.max=1, n=101, s=0.48, lambda=10)
#' hb <- hill.distrib(x.min=0, x.max=1, n=101, s=0.52, lambda=10)
#' x <- ha$x
#' h1 <- 1-pmax(ha$Hill.incr, 0.6)
#' h2 <- 1-pmax(hb$Hill.decr, 0.6)
#' h3 <- h1+h2-0.4
#' range(h3)
#' par(mfrow=c(3,1))
#' plot(x, h1, type="l", col="blue", panel.first=grid())
#' plot(x, h2,type="l", col="red", panel.first=grid())
#' plot(x, h3, type="l", col="red", panel.first=grid())
#' 
#' hc <- hill.distrib(n=101, s=0.3)
#' plot.hill(x, h1, n=101, s=0.3)
#' par(mfrow=c(1,1))
#' @export
hill.distrib <- function (
                          x.min=0, ## Min value for x
                          x.max=1, ## Max value for x
                          s=0.4, ## Threshold value
                          lambda=10, ## Exponent of the hill function
                          n=101 ## Number of points 
                          ) {

  x <- seq(from=x.min, to=x.max, length.out=n)
  
  ## Compute the increasing Hill function
  return(data.frame(x=x,
                    Hill.incr=x^lambda / (s^lambda + x^lambda),
                    Hill.decr=s^lambda / (s^lambda + x^lambda)
                    )
         )
}

plot.hill <- function (x,y,
                       s, ## Threshold value
                       n, ## Number of points 
                       xlab='concentration of A',
                       ylab='synthesis rate of B',
                       main = 'Hill function', ...) {
  plot(x,y, ylim=c(0,1), type='l', lwd=2, 
       xlab=xlab, 
       ylab=ylab, 
       main=main, panel.first=grid(), ...)
  arrows(s, 0, s, 0.5, angle = 30, length=0.05, code = 1, lwd=1)
  arrows(s, 0.5, 0, 0.5, angle = 30, length=0.05, code = 1, lwd=1)
  text(x=s,y=0.05, label="s", pos=4)
  text(x=0,y=0.55, label="0.5", pos=4)
  text(x=s*1.1, y=0.4, label=paste("f(x) = x^",n,"/(",s,"^",n,"+x^",n,")",sep=""), pos=4)
}



