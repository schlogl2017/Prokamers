# ============================================================
# Tutorial on ggplot2 basics
# by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)
# =============================================================

# library(ggplot2)
# head(diamonds)
# A dataset containing the prices and other attributes of almost 54,000 diamonds. The variables are as follows:
# price   = price in US dollars ($326–$18,823)
# carat   = weight of the diamond (0.2–5.01)
# cut     = quality of the cut (Fair, Good, Very Good, Premium, Ideal)
# colour  = diamond colour, from J (worst) to D (best)
# clarity = a measurement of how clear the diamond is (I1 (worst), SI1, SI2, VS1, VS2, VVS1, VVS2, IF (best))
# x       = length in mm (0–10.74)
# y       = width in mm (0–58.9)
# z       = depth in mm (0–31.8)
# depth   = total depth percentage = z / mean(x, y) = 2 * z / (x + y) (43–79)
# table   = width of top of diamond relative to widest point (43–95)


# head(mtcars)
# A dataset containing the prices and other attributes of almost 54,000 diamonds. The variables are as follows:
# mpg   =	Miles/(US) gallon
# cyl	  = Number of cylinders
# disp	= Displacement (cu.in.)
# hp	  = Gross horsepower
# drat	= Rear axle ratio
# wt	  = Weight (lb/1000)
# qsec	= 1/4 mile time
# vs	  = V/S
# am	  = Transmission (0 = automatic, 1 = manual)
# gear	= Number of forward gears
# carb	= Number of carburetors

# qplot(carat, price, data = diamonds)
# qplot(log(carat), log(price), data = diamonds)
# qplot(carat, x * y * z, data = diamonds)
# qplot(carat, price, data = diamonds[1:50,], colour = color)
# qplot(carat, price, data = diamonds[1:50,], shape = cut)
# qplot(carat, price, data = diamonds[1:50,], size = price)

# library(scales)
# qplot(carat, price, data = diamonds, colour = I(alpha("black", 1/200)))

# qplot(carat, price, data = diamonds, geom = c("point", "smooth"))
# qplot(carat, price, data = diamonds[1:100,], geom = c("point", "smooth"),span=0.2,se=TRUE)

# library(mgcv)
# qplot(carat, price, data = diamonds[1:100,], geom = c("point", "smooth"),method="gam", formula= y ~ s(x))
# qplot(carat, price, data = diamonds, geom = c("point", "smooth"),method="gam", formula= y ~ s(x,bs="cs"))

# library(splines)
# qplot(carat, price, data = diamonds[1:100,], geom=c("point", "smooth"),method = "lm")
# qplot(carat, price, data = diamonds[1:100,], geom=c("point", "smooth"),method = "lm", formula=y ~ poly(x,2))
# qplot(carat, price, data = diamonds[1:100,], geom=c("point", "smooth"),method = "lm", formula=y ~ ns(x,3))
# library(MASS)
# qplot(carat, price, data = diamonds[1:100,], geom=c("point", "smooth"),method = "rlm")

# qplot(color, price / carat, data = diamonds, geom = "jitter", colour = I(alpha("black", 1 / 10)))
# qplot(color, price / carat, data = diamonds, geom = "boxplot") 

# qplot(carat, data = diamonds, geom = "histogram")
# qplot(carat, data = diamonds, geom = "density")

# qplot(carat, data = diamonds, geom = "histogram", fill = color)
# qplot(carat, data = diamonds, geom = "density", colour = color)

# qplot(carat, data = diamonds, facets = . ~ color, geom = "histogram", binwidth = 0.1, xlim = c(0, 3))

# p<-qplot(carat, price, data = diamonds[1:50,], colour = color)
# summary(p)

# p <- ggplot(diamonds, aes(carat, price, colour = cut))
# p <- p + layer(geom = "point")
# p


# p<-ggplot(diamonds,aes(carat,price))+geom_point(colour="darkblue")
# p
# p<-ggplot(diamonds,aes(carat,price))+geom_point(aes(colour="darkblue"))
# p

# head(Oxboys)
# These data are described in Goldstein (1987) as data on the height of a selection of boys from Oxford, England versus a standardized age.
# Subject = an ordered factor giving a unique identifier for each boy in the experiment
# age     = a numeric vector giving the standardized age (dimensionless)
# height  = a numeric vector giving the height of the boy (cm)
# Occasion= an ordered factor - the result of converting age from a continuous variable to a count so these slightly unbalanced data can be analyzed as balanced.

# p <- ggplot(Oxboys, aes(age, height)) + geom_line()
# p
# p <- ggplot(Oxboys, aes(age, height, group = Subject)) + geom_line()
# p

# p <- ggplot(Oxboys, aes(age, height, group = Subject)) + geom_line()
# p + geom_smooth(method="lm", size=2, se=F)
# p <- ggplot(Oxboys, aes(age, height, group = Subject)) + geom_line()
# p + geom_smooth(aes(group="dummy"),method="lm", size=2, se=F)

# d <- ggplot(diamonds, aes(carat)) + xlim(0, 3)
# d + stat_bin(aes(ymax = ..count..), binwidth = 0.1, geom = "area")
# d + stat_bin(
#   aes(size = ..density..), binwidth = 0.1,
#   geom = "point", position="identity"
# )
# d + stat_bin(
#   aes(y = 1, fill = ..count..), binwidth = 0.1,
#   geom = "tile", position="identity"
# )

# df <- data.frame(x = c(3, 1, 5), y = c(2, 4, 6), label = c("a","b","c"))
# p <- ggplot(df, aes(x, y, label = label)) + xlab(NULL) + ylab(NULL)
# p + geom_point() + ggtitle("geom_point")
# p + geom_bar(stat="identity") +ggtitle("geom_bar(stat=\"identity\")")
# p + geom_line() + ggtitle("geom_line")
# p + geom_area() + ggtitle("geom_area")
# p + geom_path() + ggtitle("geom_path")
# p + geom_text() + ggtitle("geom_text")
# p + geom_tile() + ggtitle("geom_tile")
# p + geom_polygon() + ggtitle("geom_polygon")

# qplot(carat, price, data = diamonds[1:100,], colour = color) 
# qplot(carat, price, data = diamonds[1:100,], colour = color)  + scale_color_hue("Diamond Colour")
# qplot(carat, price, data = diamonds[1:100,], colour = color)  + 
#   scale_color_hue("Diamond Colour", breaks=c("D","E","F"))
# qplot(carat, price, data = diamonds[1:100,], colour = color)  + 
#   scale_color_hue("Diamond Colour", breaks=c("D","E","F"), labels=c("D grade","E grade","F grade"))
# qplot(carat, price, data = diamonds[1:100,], colour = color)  + 
#   scale_color_hue("Diamond Colour", limits=c("D","E","F"), labels=c("D grade","E grade","F grade"))

# qplot(log10(carat), log10(price), data = diamonds)
# qplot(carat, price, data = diamonds) + scale_x_log10() + scale_y_log10()
# qplot(carat, price, data = diamonds) + scale_x_continuous(trans="log10") + scale_y_log10()

# The surface of a 2d density estimate of the faithful dataset (Azzalini and Bowman, 1990), 
# which records the waiting time between eruptions and during of each eruption for the Old Faithful geyser in Yellowstone Park.
# f2d <- with(faithful, MASS::kde2d(eruptions, waiting, h = c(1, 10), n = 50))
# df <- with(f2d, cbind(expand.grid(x, y), as.vector(z)))
# names(df) <- c("eruptions", "waiting", "density")
# erupt <- ggplot(df, aes(waiting, eruptions, fill = density)) +
#   geom_tile() +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0))
# erupt + scale_fill_gradient(limits = c(0, 0.04))
# erupt + scale_fill_gradient(limits = c(0, 0.04), low="white", high="black")
# erupt + scale_fill_gradient2(limits = c(-0.04, 0.04),
#                              midpoint = mean(df$density))


# library(colorspace)
# fill_gradn <- function(pal) {
#   scale_fill_gradientn(colours = pal(7), limits = c(0, 0.04))
# }
# erupt + fill_gradn(rainbow_hcl)
# erupt + fill_gradn(diverge_hcl)
# erupt + fill_gradn(heat_hcl)

# library(RColorBrewer)
# RColorBrewer::display.brewer.all()
# p<-qplot(price, data = diamonds, colour = color,bindwidth=10,fill=color)
# p+scale_fill_brewer(pal = "Set1")
# p+scale_fill_brewer(pal = "Set2")
# p+scale_fill_brewer(pal = "Pastel1")