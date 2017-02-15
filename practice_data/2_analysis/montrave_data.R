## 03 Feb 2017
## Exploring function for detection function using "Distance" packages and examples online
library(Distance)
library(knitr)

### 1- using point transect data of montrave songbird case study

setwd("/Users/SantiagoD/Desktop/Distance sampling in R")
birds.point <- read.csv("montrave-point.csv")

### Following PDF file, we have to "clean the data", ie filter only species of interest and include effort data
# we have to create a function first

findFirstLast <- function(myDF, findFirst=TRUE) {
  # myDF should be a data frame or matrix
  # By default, this function finds the first occurence of each unique value in a column
  # If instead we want to find last, set findFirst to FALSE.
  # This will give `maxOrMin` a value of -1 finding the min of the negative indexes
  # is the same as finding the max of the positive indexes
  maxOrMin <- ifelse(findFirst, 1, -1)
  # For each column in myDF, make a list of all unique values (`levs`) and
  # iterate over that list, finding the min (or max) of all the indexes
  # of where that given value appears within the column
  apply(myDF, 2, function(colm) {
    levs <- unique(colm)
    sapply(levs, function(lev) {
      inds <- which(colm==lev)
      ifelse(length(inds)==0, NA, maxOrMin*min(inds*maxOrMin) )
    })
  })
}

#then clean data for robins only

robins.pt <- subset(birds.point, species=="r")
pt.with.robins <- unique(robins.pt$Sample.Label) # point IDs with robins
all.points <- seq(1: length(unique(birds.point$Sample.Label))) # point IDs all transects
# empty points will occupy the last few elements of the union statement
empty.points <- union(pt.with.robins, all.points)[(length(pt.with.robins)+1):length(all.points)]
# what was effort on empty points?
all.point.effort <- findFirstLast(birds.point)$Sample.Label
effort.empty.points <- birds.point$Effort[all.point.effort[empty.points]]
# append 'blank' rows to bottom of species-specific data set
# retain effort of transects with no sightings
empty <- NULL
for (i in 1:length(empty.points)) {
  blank <- cbind(birds.point[1,1:2], empty.points[i], effort.empty.points[i])
  empty <- rbind(empty, blank)
}
empty[,c("A","B","C")] <- NA # pad out columns that are blank (distance, species, visit)
names(empty) <- names(birds.point)
robins.pt <- rbind(robins.pt, empty)
robins.pt <- robins.pt[order(robins.pt$Sample.Label), ]

## now we have a dataset that includes both detections and no detections effort, for the abundance
# estimates... in my case with PH, I might not need to do this, since JJ will send me a more complete
# file with the total number of transects, so I can include both det, and no det from the beginning

##### NOW - fitting the detection function
robin.pt.hn.herm <- ds(robins.pt, truncation=110, transect="point",
                       key="hn", adjustment="herm", convert.units=.01) # change in conversion
robin.pt.uni.cos <- ds(robins.pt, truncation=110, transect="point",
                       key="unif", adjustment="cos", convert.units=.01)
robin.pt.haz.simp <- ds(robins.pt, truncation=110, transect="point",
                        key="hr", adjustment="poly", convert.units=.01)
pt.model.results <- rbind(robin.pt.uni.cos$dht$individuals$D, robin.pt.haz.simp$dht$individuals$D,
                          robin.pt.hn.herm$dht$individuals$D)

#goodness of fitness of model selected
robin.breaks <- c(0,22.5,32.5,42.5,52.5,62.5,77.5,110) # from Section 5.2.3.3
fit.pt.haz.simp <- ddf.gof(robin.pt.haz.simp$ddf, breaks=robin.breaks,
                           main="Montrave robin point transect data, hazard rate detection function")

#another visualization

hazard <- function(y, sigma, shape) {
  key <- 1-exp(-(y/sigma)^(-shape))
  return(key)
}
haz <- function(distances, ddfobj, point=TRUE) {
  sigma <- exp(ddfobj$par[2]) # contrary to hn models
  shape <- exp(ddfobj$par[1])
  pr.detect <- hazard(distances, sigma, shape)
  if (point) pr.detect <- pr.detect * distances
  return(pr.detect)
}
pdf.point <- function(ddf.obj, mybreaks, ...) {
  # ddf.obj is produced by a call to ds()
  # result is a plot
  # sort out the pdf
  upperbnd <- ddf.obj$meta.data$int.range[2]
  distances <- seq(0,upperbnd,length.out = 75)
  if (ddf.obj$ds$aux$ddfobj$type=="hn") {
    detfn.line <- hnherm(distances, point=TRUE, ddfobj = ddf.obj)
    df.integral <- integrate(hnherm, lower=0, upper=upperbnd, point=TRUE,
                             ddfobj=ddf.obj)[1]$value
  }
  if (ddf.obj$ds$aux$ddfobj$type=="hr") {
    detfn.line <- haz(distances, point=TRUE, ddfobj=ddf.obj)
    df.integral <- integrate(haz, lower=0, upper=upperbnd, point=TRUE,
                             ddfobj=ddf.obj)[1]$value
  }
  detfn.line <- detfn.line / df.integral
  # now for the bars
  hist.dist <- hist(ddf.obj$data$distance, breaks=mybreaks, plot=FALSE)
  # picture
  plot(hist.dist, freq=FALSE, xlab="Distance", ylab="Probability density",...)
  lines(x=distances, y=detfn.line,...)
  box()
  return(hist.dist)
}
point.plot <- pdf.point(robin.pt.haz.simp$ddf, mybreaks=robin.breaks,
                        main="Montrave robin point counts, hazard pdf", lwd=2)


chirow <- c(fit.pt.haz.simp$chisquare$chi1$chisq, fit.pt.haz.simp$chisquare$chi1$p)
ksrow <- c(fit.pt.haz.simp$dsgof$ks$Dn, fit.pt.haz.simp$dsgof$ks$p)
cvmrow <- c(fit.pt.haz.simp$dsgof$CvM$W, fit.pt.haz.simp$dsgof$CvM$p)
mytable <- rbind(chirow, ksrow, cvmrow)
rownames(mytable) <- c("Chi-square test", "K-S test", "CvM test")
kable(mytable, col.names = c("Test statistic", "P-value"), digits=3,
      caption="Goodness of fit statistics, Montrave robin point transects, hazard rate model")

# Inelegant way to build first column model names
pt.model.results[,1] <- as.character(model.results[,1])
pt.model.results[1,1] <- "Unif.cosine"
pt.model.results[2,1] <- "Hazard rate"
pt.model.results[3,1] <- "Half-norm. Hermite"
kable(pt.model.results[,1:6], digits=3,
      caption="Density estimates under Uniform/cos, Hazard rate, and half-normal Hermite models (Table)")

