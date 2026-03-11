### Gulf of Alaska Ricker models

# Packages:
library(FSA)
library(car) # Before dplyr to reduce conflicts with MASS 
library(dplyr)
library(magrittr)
library(plotrix)
library(nlstools)
library(lsmeans)

# Load in data
load("~/GitHub/Extended-Independent-Project/Gulf_Of_Alaska.RData")

# Extract Ricker function
ricker <- srFuns("Ricker")

# Arrowtooth flounder

aflounder <- GoA.data[GoA.data$species == "Arrowtooth flounder",] # extracts arrowtooth flounder data
aflounder <- aflounder[,-2] # removes species name

AFstart <- srStarts(recruits~ssb, data = aflounder, type = "Ricker") # calculate starting values for parameters
AFstart

AFmodel <- nls(log.recruits~log(ricker(ssb,a,b)), data = aflounder, start = AFstart) # fits model using nls

cbind(estimates = coef(AFmodel), confint(AFmodel))

AFbootR <- nlsBoot(AFmodel) # calculates confidence intervals with bootstrap method
cbind(estimates = coef(AFmodel), confint(AFbootR))

AFx <- seq(min(aflounder$ssb), max(aflounder$ssb), length.out = 199) # produce values of S for prediction
AFpredR <- ricker(AFx, a = coef(AFmodel)) # predicted mean R
AFLCI <- AFUCI <- numeric(length(AFx)) # create empty vectors length x
for(i in 1:length(AFx)) {
  tmp <- apply(AFbootR$coefboot, MARGIN = 1, FUN = ricker, S = AFx[i])
  AFLCI[i] <- quantile(tmp, 0.025)
  AFUCI[i] <- quantile(tmp, 0.975)
}

AFylmts <- range(c(AFpredR, AFLCI, AFUCI, aflounder$recruits))
AFxlmts <- range(c(AFx, aflounder$ssb))

# Flathead sole

fsole <- GoA.data[GoA.data$species == "Flathead sole",] # extracts flathead sole data
fsole <- fsole[,-2] # removes species name

FSstart <- srStarts(recruits~ssb, data = fsole, type = "Ricker") # calculate starting values for parameters
FSstart

FSmodel <- nls(log.recruits~log(ricker(ssb,a,b)), data = fsole, start = FSstart) # fits model using nls

cbind(estimates = coef(FSmodel), confint(FSmodel))

FSbootR <- nlsBoot(FSmodel) # calculates confidence intervals with bootstrap method
cbind(estimates = coef(FSmodel), confint(FSbootR)) # fit did not converge twice during bootstrapping?

FSx <- seq(min(fsole$ssb), max(fsole$ssb), length.out = 199) # produce values of S for prediction
FSpredR <- ricker(FSx, a = coef(FSmodel)) # predicted mean R
FSLCI <- FSUCI <- numeric(length(FSx)) # create empty vectors length x
for(i in 1:length(FSx)) {
  tmp <- apply(FSbootR$coefboot, MARGIN = 1, FUN = ricker, S = FSx[i])
  FSLCI[i] <- quantile(tmp, 0.025)
  FSUCI[i] <- quantile(tmp, 0.975)
}

FSylmts <- range(c(FSpredR, FSLCI, FSUCI, fsole$recruits))
FSxlmts <- range(c(FSx, fsole$ssb))

# Pacific cod

pcod <- GoA.data[GoA.data$species == "Pacific cod",] # extracts pacific cod data
fsole <- fsole[,-2] # removes species name

FSstart <- srStarts(recruits~ssb, data = fsole, type = "Ricker") # calculate starting values for parameters
FSstart

FSmodel <- nls(log.recruits~log(ricker(ssb,a,b)), data = fsole, start = FSstart) # fits model using nls

cbind(estimates = coef(FSmodel), confint(FSmodel))

FSbootR <- nlsBoot(FSmodel) # calculates confidence intervals with bootstrap method
cbind(estimates = coef(FSmodel), confint(FSbootR)) # fit did not converge twice during bootstrapping?

FSx <- seq(min(fsole$ssb), max(fsole$ssb), length.out = 199) # produce values of S for prediction
FSpredR <- ricker(FSx, a = coef(FSmodel)) # predicted mean R
FSLCI <- FSUCI <- numeric(length(FSx)) # create empty vectors length x
for(i in 1:length(FSx)) {
  tmp <- apply(FSbootR$coefboot, MARGIN = 1, FUN = ricker, S = FSx[i])
  FSLCI[i] <- quantile(tmp, 0.025)
  FSUCI[i] <- quantile(tmp, 0.975)
}

FSylmts <- range(c(FSpredR, FSLCI, FSUCI, fsole$recruits))
FSxlmts <- range(c(FSx, fsole$ssb))

# Pacific ocean perch

# Rock sole

# Walleye pollock

### Combined plot:

# Set up plot so all 6 can be displayed
par(mfrow = c(2,3)) 

# Arrowtooth flounder

plot(recruits~ssb, data = aflounder, 
     xlim = AFxlmts, ylim = AFylmts, 
     col = "white", 
     ylab = "Recruits", xlab = "SSB",
     main = "Ricker stock-recruitment curve for arrowtooth flounder")
polygon(c(AFx, rev(AFx)), c(AFLCI,rev(AFUCI)), col = "gray80", border = NA)
points(recruits~ssb, data = aflounder, pch = 19, col = rgb(0,0,0,1/2))
lines(AFpredR~AFx, lwd = 2)

# Flathead sole

plot(recruits~ssb, data = fsole, 
     xlim = FSxlmts, ylim = FSylmts, 
     col = "white", 
     ylab = "Recruits", xlab = "SSB",
     main = "Ricker stock-recruitment curve for flathead sole")
polygon(c(FSx, rev(FSx)), c(FSLCI,rev(FSUCI)), col = "gray80", border = NA)
points(recruits~ssb, data = fsole, pch = 19, col = rgb(0,0,0,1/2))
lines(FSpredR~FSx, lwd = 2)

# Pacific cod

# Pacific ocean perch

# Rock sole

# Walleye pollock
