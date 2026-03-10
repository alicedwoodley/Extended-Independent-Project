### Arrowtooth flounder Ricker model

library(FSA)
library(car) # Before dplyr to reduce conflicts with MASS 
library(dplyr)
library(magrittr)
library(plotrix)
library(nlstools)
library(lsmeans)

# No extra variables

load("~/GitHub/Extended-Independent-Project/Gulf_Of_Alaska.RData")

aflounder <- GoA.data[GoA.data$species == "Arrowtooth flounder",] # extracts AF data
aflounder <- aflounder[,-2] # removes species name

start <- srStarts(recruits~ssb, data = aflounder, type = "Ricker") # calculate starting values for parameters
start

ricker <- srFuns("Ricker") # extract ricker function

afmodel <- nls(log.recruits~log(ricker(ssb,a,b)), data = aflounder, start = start) # fits model using nls

cbind(estimates = coef(afmodel), confint(afmodel))

bootR <- nlsBoot(afmodel) # calculates confidence intervals with bootstrap method
cbind(estimates = coef(afmodel), confint(bootR))

x <- seq(min(aflounder$ssb), max(aflounder$ssb), length.out = 199) # produce values of S for prediction
predR <- ricker(x, a = coef(afmodel)) # predicted mean R
LCI <- UCI <- numeric(length(x)) # create empty vectors length x
for(i in 1:length(x)) {
  tmp <- apply(bootR$coefboot, MARGIN = 1, FUN = ricker, S = x[i])
  LCI[i] <- quantile(tmp, 0.025)
  UCI[i] <- quantile(tmp, 0.975)
}
ylmts <- range(c(predR, LCI, UCI, aflounder$recruits))
xlmts <- range(c(x, aflounder$ssb))

plot(recruits~ssb, data = aflounder, xlim = xlmts, ylim = ylmts, col = "white", ylab = "Recruits", xlab = "SSB")
polygon(c(x, rev(x)), c(LCI,rev(UCI)), col = "gray80", border = NA)
points(recruits~ssb, data = aflounder, pch = 19, col = rgb(0,0,0,1/2))
lines(predR~x, lwd = 2)

# With additional variables

pcod <- GoA.data[GoA.data$species == "Pacific cod",] # extract pacific cod data
wpollock <- GoA.data[GoA.data$species == "Walleye pollock",] # extract walleye pollock data

aflounder <- cbind(aflounder, pcod$ssb, wpollock$ssb) # adds stock size of other species to aflounder
colnames(aflounder)[6] <- "pcodssb" ; colnames(aflounder)[7] <- "wpollockssb"

newricker <- function(S, X1, X2, a, b = NULL, c = NULL, d = NULL){ # S = af stock, X1 = walleye pollock, X2 = pacific cod
  if(length(a)>1){
    d <- a[4]
    c <- a[3]
    b <- a[2]
    a <- a[1]
  }
  a*S*exp(-b*S+c*X1-d*X2) 
}

# find starting values through linearising and using lm()
linear <- lm(log(recruits/ssb)~ ssb + wpollockssb + pcodssb, data = aflounder)
start2 <- coef(linear) # extract coefficients
start2 <- list(a = exp(start2[[1]]), b = start2[[2]], c = start2[[3]], d = start2[[4]]) # prep to put in model

multiafmodel <- nls(log.recruits~log(newricker(ssb, wpollockssb, pcodssb, a, b, c, d)), data = aflounder, start = start2)
bootR2 <- nlsBoot(multiafmodel) # calculates confidence intervals with bootstrap method - fewer warnings?
cbind(estimates = coef(multiafmodel), confint(bootR2)) # odd that c is negative?

x2 <- seq(min(aflounder$ssb), max(aflounder$ssb), length.out = 38) # produce values of S for prediction
predRnew <- newricker(x2, aflounder$wpollockssb, aflounder$pcodssb, a = coef(multiafmodel)[[1]], b = coef(multiafmodel)[[2]], c = coef(multiafmodel)[[3]], d = coef(multiafmodel)[[4]] ) # predicted mean R
LCI <- UCI <- numeric(length(x2)) # create empty vectors length x
for(i in 1:length(x2)) {
  tmp <- apply(bootR2$coefboot, MARGIN = 1, function(coefs) {
    newricker(S = x2[i], X1 = aflounder$wpollockssb[i], X2 = aflounder$pcodssb[i], a = coefs)
  })
  LCI[i] <- quantile(tmp, 0.025)
  UCI[i] <- quantile(tmp, 0.975)
}
ylmts <- range(c(predRnew, LCI, UCI, aflounder$recruits))
xlmts <- range(c(x2, aflounder$ssb))

plot(recruits~ssb, data = aflounder, xlim = xlmts, ylim = ylmts, col = "white", main = "Ricker model of Arrowtooth Flounder stock-recruitment including SSB of predators and prey", cex.main = 0.8, ylab = "Recruits", xlab = "SSB")
polygon(c(x2, rev(x2)), c(LCI,rev(UCI)), col = "gray80", border = NA)
points(recruits~ssb, data = aflounder, pch = 19, col = rgb(0,0,0,1/2))
lines(predRnew~x2, lwd = 2)

# Model test - multispecies one vs single species







