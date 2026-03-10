### Stock-recruitment of single species

library(FSA)
library(car) # Before dplyr to reduce conflicts with MASS 
library(dplyr)
library(magrittr)
library(plotrix)
library(nlstools)
library(lsmeans)

load("~/GitHub/Extended-Independent-Project/Gulf_Of_Alaska.RData")

wpollock <- GoA.data[GoA.data$species == "Walleye pollock",] # extracts walleye pollock data
wpollock <- wpollock[,-2] # removes species name column

start <- srStarts(recruits~ssb, data = wpollock, type = "Ricker") # calculate starting values for parameters
start # warning: b not positive, likely poor starting value

ricker <- srFuns("Ricker") # extract ricker function from package

wpmodel <- nls(log.recruits~log(ricker(ssb,a,b)), data = wpollock, start = start) # fit ricker model using log-log

cbind(estimates = coef(wpmodel), confint(wpmodel)) # shows parameter estimates and confidence intervals

bootR <- nlsBoot(wpmodel) # calculates confidence intervals with bootstrap method - lots of warnings
cbind(estimates = coef(wpmodel), confint(bootR))

ind <- srFuns("independence") # extracts function for model with no density dependence (b = 0)
indstart <- srStarts(recruits~ssb, data = wpollock, type = "independence") # extract starting values 
indmodel <- nls(log.recruits ~ log(ind(ssb,a)), data = wpollock, start = indstart) # fit model with no density dependence

# could you do weighted nls so that one annoying point doesn't have as much of an effect?

lrt(indmodel, com = wpmodel) # likelihood ratio test
# results suggest simpler model is still an adequate fit of the data - high p value

x <- seq(min(wpollock$ssb), max(wpollock$ssb), length.out = 199) # produce values of S for prediction
predR <- ricker(x, a = coef(wpmodel)) # predicted mean R
LCI <- UCI <- numeric(length(x)) # create empty vectors length x
for(i in 1:length(x)) {
  tmp <- apply(bootR$coefboot, MARGIN = 1, FUN = ricker, S = x[i])
  LCI[i] <- quantile(tmp, 0.025)
  UCI[i] <- quantile(tmp, 0.975)
}
ylmts <- range(c(predR, LCI, UCI, wpollock$recruits))
xlmts <- range(c(x, wpollock$ssb))

plot(recruits~ssb, data = wpollock, xlim = xlmts, ylim = ylmts, col = "white", ylab = "Recruits", xlab = "SSB")
polygon(c(x, rev(x)), c(LCI,rev(UCI)), col = "gray80", border = NA)
points(recruits~ssb, data = wpollock, pch = 19, col = rgb(0,0,0,1/2))
lines(predR~x, lwd = 2)

# incorporating pacific ocean perch stock as an explanatory variable

poperch <- GoA.data[GoA.data$species == "Pacific ocean perch",] # find pop data
poperch <- poperch[,-2] 

wpollock <- cbind(wpollock, poperch$ssb) # adds stock size of pop to the wpollock dataframe
colnames(wpollock)[6] <- "popssb"

newricker <- function(S, X, a, b = NULL, c = NULL){ # S = wp stock, X = pop stock
  if(length(a)>1){
    c <- a[3]
    b <- a[2]
    a <- a[1]
  }
  a*S*exp(-b*S+c*X) # creates new ricker function
}

# find starting values through linearising and using lm()
linear <- lm(log(recruits/ssb)~ ssb + popssb, data = wpollock)
start2 <- coef(linear) # extract coefficients
start2 <- list(a = exp(start2[[1]]), b = start2[[2]], c = start2[[3]]) # prep to put in model

modelwpop <- nls(log.recruits~log(newricker(ssb, popssb, a, b, c)), data = wpollock, start = start2)

bootR2 <- nlsBoot(modelwpop) # calculates confidence intervals with bootstrap method - lots of warnings
cbind(estimates = coef(modelwpop), confint(bootR2))

# prediction with existing pop ssb

x2 <- seq(min(wpollock$ssb), max(wpollock$ssb), length.out = 38) # produce values of S for prediction
predRnew <- newricker(x2, wpollock$popssb, a = coef(modelwpop)[[1]], b = coef(modelwpop)[[2]], c = coef(modelwpop)[[3]] ) # predicted mean R
LCI <- UCI <- numeric(length(x2)) # create empty vectors length x
for(i in 1:length(x2)) {
  tmp <- apply(bootR2$coefboot, MARGIN = 1, FUN = newricker, S = x2[i], a = coef(modelwpop)[[1]], b = coef(modelwpop)[[2]], c = coef(modelwpop)[[3]])
  LCI[i] <- quantile(tmp, 0.025)
  UCI[i] <- quantile(tmp, 0.975)
}
ylmts <- range(c(predRnew, LCI, UCI, wpollock$recruits))
xlmts <- range(c(x2, wpollock$ssb))

plot(recruits~ssb, data = wpollock, xlim = xlmts, ylim = ylmts, col = "white", main = "Ricker model of Walleye Pollock stock-recruitment including Pacific Ocean Perch SSB", cex.main = 0.8, ylab = "Recruits", xlab = "SSB")
polygon(c(x2, rev(x2)), c(LCI,rev(UCI)), col = "gray80", border = NA)
points(recruits~ssb, data = wpollock, pch = 19, col = rgb(0,0,0,1/2))
lines(predRnew~x2, lwd = 2)

plot(ssb~year, data = poperch)



# 