## new files to explore functions in Distance package for point count

## 10 Feb 2017  
## Continued on 14 Feb 2017

## Exploring function for detection function using "Distance" packages and examples online
library(Distance)
library(knitr)

### NOW using "amakihi" data and examples from webpage, this one should use also covariates to fit the detection function
### From paper miller etal 2016

data("amakihi")

### first model Hazard-rate with one covariate
amakihi_hr_obs <- ds(amakihi, truncation = 82.5, transect = "point", 
                     key = "hr", formula = ~obs)

### second model Hazard-rate with two covariates
amakihi_hr_obs_mas <- ds(amakihi, truncation = 82.5, transect = "point", 
                     key = "hr", formula = ~obs + mas)

### I want to visualize????, but they dont explain how to build graph with visualization on paper

# Summary of model
summary(amakihi_hr_obs)

# Goodness of fit
amakihi_hn <- ds(amakihi, truncation = 82.5, transect = "point", key = "hn", adjustment = NULL)
gof_ds(amakihi_hn)                 
gof_ds(amakihi_hr_obs)

# model selection
summarize_ds_models(amakihi_hn, amakihi_hr_obs, amakihi_hr_obs_mas)
# only shows a table comparing the models and the diffence in AIC

# Estimating abundance and variance
# Paper inly provided example for line transect... ???



############################################################################
##### now following the other PDf for the "same?" data
scalemin <- sd(amakihi$mas, na.rm = TRUE)

### plots to explore behavior and distribution of data and covariates
par(mfrow=c(3,2))
hist(amakihi$distance, breaks=seq(0,260), main="",xlab = "Distance") 
boxplot(amakihi$distance[!amakihi$obs == ""]~amakihi$obs[!amakihi$obs == ""],
        xlab="Observer", ylab="Distance(m)")
hist(amakihi$distance[amakihi$distance<82.5], breaks=33,
     main="",xlab = "Distance") 
plot(amakihi$mas, amakihi$distance,
     xlab="Minutes after sunrise", ylab="Distance(m)",
     pch=19, cex=0.6) 
abline(reg=lm(amakihi$distance~amakihi$MAS),lwd=2) 
hist(amakihi$distance[amakihi$distance<82.5], breaks=10,
     main="",xlab = "Distance") 
boxplot(amakihi$distance~amakihi$has,
        xlab="Hour", ylab="Distance(m)")
par(mfrow=c(1,1))

#Adjusting raw covariates
str(amakihi$obs)

amakihi$has <- as.factor(amakihi$has)
amakihi$obs <- relevel(amakihi$obs, ref="TKP")
amakihi$has <- relevel(amakihi$has, ref="5")
amakihi$mas <- amakihi$mas/sd(amakihi$mas, na.rm=TRUE) #divides measure in sd to obtain a scale comparable to other covariates

# FITTING detection function and goodness of fit and plotting
# wrapping function

fit.and.assess <- function(data=amakihi, key="hn",
                           adj="cos", order=NULL,
                           cov=c("mas"), truncation=82.5,
                           plotddf = FALSE, ...) {
  # wrapper function to combine detection function fitting as well as
  # goodness of fit testing into a single function
  #
  #Input - name of dataframe, key/adjustment combination, order of adjustment terms,
  # and a vector of covariates to include in detection function
  # specifically for this example, I limit number of covariates to 2,
  # truncation distance (default value for amakihi),
  # flag plotddf to signal whether a plot of the detfn is required
  #
  #Output - named list containing the columns of Table 5.7 (with the exception
  # of the effective detection radius). List will subsequently be
  # post-processed to duplicate the table.
  if (length(cov)>2) stop("Function can take no more than 2 covariates")
  mystring <- ifelse(length(cov)==2,
                     paste(cov[1],cov[2],sep="+"),
                     cov[1])
  covariates <- as.formula(paste("~", mystring))
  model.fit <- ds(data, key=key, adjustment=adj,
                  formula=covariates, transect="point", order=order,
                  truncation = truncation)
  if (plotddf) plot(model.fit, pl.den=0, pch=19,
                    cex=0.7, main=paste(key, adj, cov))
  aic <- model.fit$ddf$criterion
  pars <- length(model.fit$ddf$par)
  gof <- ddf.gof(model.fit$ddf)
  chisq.p <- gof$chisquare$chi1$p
  cvm.p <- gof$dsgof$CvM$p
  ks.p <- gof$dsgof$ks$p
  findings <- list(formula=mystring, key=key, aic=aic, pars=pars,
                   chisq.p=chisq.p, cvm.p=cvm.p, ks.p=ks.p)
  return(findings)
}


## Candidate models
hour.hn <- fit.and.assess(key="hn", order=2, cov=c("has"))
hour.hr <- fit.and.assess(key="hr", adj=NULL, cov=c("has"))
obs.hour.hn <- fit.and.assess(key="hn", order=0, cov=c("obs","has")) # fails if order>0
obs.hour.hr <- fit.and.assess(key="hr", adj=NULL, cov=c("obs","has"))
obs.min.hn <- fit.and.assess(key="hn",order=0, cov=c("obs","mas")) # fails if order>0
obs.min.hr <- fit.and.assess(key="hr", adj=NULL, cov=c("obs","mas"))
min.hn <- fit.and.assess(key="hn", order=0, cov=c("mas")) # fails if order>0
min.hr <- fit.and.assess(key="hr", adj=NULL, cov=c("mas"))
obs.hn <- fit.and.assess(key="hn", order=0, cov=c("obs")) # fails if order>0
obs.hr <- fit.and.assess(key="hr", adj=NULL, cov=c("obs"))
hn <- fit.and.assess(key="hn", order=(2:5), cov=c(1))
hr <- fit.and.assess(key="hr", adj=NULL, cov=c(1))
uni <- fit.and.assess(key="uni", cov=c(1))

## summarize model selection
library(data.table)
library(knitr)
mymodels <- list(hour.hn, hour.hr, obs.hour.hn, obs.hour.hr,
                 obs.min.hn, obs.min.hr, min.hn, min.hr,
                 obs.hn, obs.hr, hn, hr, uni)
results.frame<- rbindlist(mymodels)
results.frame <- results.frame[order(results.frame$aic), ] # sort by increasing AIC
smallest <- results.frame[1, aic] # note syntax for data.table
results.frame$deltaaic <- results.frame$aic - smallest
setcolorder(results.frame, c("formula","key","pars","deltaaic","chisq.p","cvm.p","ks.p","aic"))
kable(results.frame[,.(formula, key, pars, deltaaic, chisq.p, cvm.p, ks.p)], digits=c(0,0,0,1,0,3,3))


### Missing last part, but no estimation of population size was done... 
### it was only detection function with different models
