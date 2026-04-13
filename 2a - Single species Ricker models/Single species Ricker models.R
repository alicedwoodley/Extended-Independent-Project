##### Gulf of Alaska Ricker models #####

### Packages:
library(FSA)
library(car) # Before dplyr to reduce conflicts with MASS 
library(dplyr)
library(magrittr)
library(plotrix)
library(nlstools)
library(lsmeans)
library(scales)

### Load in data
load("~/GitHub/Extended-Independent-Project/Gulf_Of_Alaska.RData")

### Extract Ricker function
ricker <- srFuns("Ricker")

### Set seed for bootstrapping
set.seed(123)

### Basic plots of a Beverton-Holt and Ricker curve for report section 2

# Extract Beverton-Holt function
bh <- srFuns("BevertonHolt", msg = T)
example <- seq(from = 0, to = 1000, by = 0.5) # Generate example SSB

par(mfrow = c(1,2))
plot(bh(example, a = 0.5, b = 0.005)~example,
     type = "l",
     lwd = 2,
     main = "(a) Beverton-Holt",
     col = "red",
     xlab = "Stock size",
     ylab = "Recruitment")

plot(ricker(example, a = 0.5, b = 0.005)~example,
     type = "l",
     lwd = 2,
     main = "(b) Ricker",
     col = "blue",
     xlab = "Stock size",
     ylab = "Recruitment")

### Arrowtooth flounder

# Extract arrowtooth flounder data
aflounder <- GoA.data[GoA.data$species == "Arrowtooth flounder",]
aflounder <- aflounder[,-2] # removes species name

# Calculate starting values for parameter estimation
AFstart <- srStarts(recruits~ssb, data = aflounder, type = "Ricker")

# Try some "bad" starting values
AFstart_bad <- list(a = 1000, b = 0.00065)

# Fit the model using nls
AFmodel <- nls(log.recruits~log(ricker(ssb,a,b)), data = aflounder, start = AFstart_bad)
summary(AFmodel) # Note: nls has performed 5 iterations from the "bad" starting values

# Show how the calculated starting values produce the almost ideal parameters before nls!
cbind(AFstart, "nls estimates" = coef(AFmodel))
# Note: all other nls models will now use the calculated starting values and so will perform 0 iterations

# Calculate 95% confidence intervals for parameter estimates using the bootstrap method
AFbootR <- nlsBoot(AFmodel)
cbind(estimates = coef(AFmodel), confint(AFbootR))

# Produce values of S to predict new values of R
AFx <- seq(min(aflounder$ssb), max(aflounder$ssb), length.out = 199)
AFpredR <- ricker(AFx, a = coef(AFmodel))
AFLCI <- AFUCI <- numeric(length(AFx))
for(i in 1:length(AFx)) { # stores a 95% confidence interval for each predicted value of R
  tmp <- apply(AFbootR$coefboot, MARGIN = 1, FUN = ricker, S = AFx[i])
  AFLCI[i] <- quantile(tmp, 0.025)
  AFUCI[i] <- quantile(tmp, 0.975)
}

# Create axis limits for plot
AFylmts <- range(c(AFpredR, AFLCI, AFUCI, aflounder$recruits))
AFxlmts <- range(c(AFx, aflounder$ssb))

### Flathead sole

# Extract flathead sole data
fsole <- GoA.data[GoA.data$species == "Flathead sole",]
fsole <- fsole[,-2] # removes species name

# Calculate starting values for parameter estimation
FSstart <- srStarts(recruits~ssb, data = fsole, type = "Ricker")
FSstart

# Fit the model using nls
FSmodel <- nls(log.recruits~log(ricker(ssb,a,b)), data = fsole, start = FSstart)

# Calculate 95% confidence intervals for parameter estimates using the bootstrap method
FSbootR <- nlsBoot(FSmodel)
cbind(estimates = coef(FSmodel), confint(FSbootR))

# Produce values of S to predict new values of R
FSx <- seq(min(fsole$ssb), max(fsole$ssb), length.out = 199)
FSpredR <- ricker(FSx, a = coef(FSmodel))
FSLCI <- FSUCI <- numeric(length(FSx))
for(i in 1:length(FSx)) { # stores a 95% confidence interval for each predicted value of R
  tmp <- apply(FSbootR$coefboot, MARGIN = 1, FUN = ricker, S = FSx[i])
  FSLCI[i] <- quantile(tmp, 0.025)
  FSUCI[i] <- quantile(tmp, 0.975)
}

# Create axis limits for plot
FSylmts <- range(c(FSpredR, FSLCI, FSUCI, fsole$recruits))
FSxlmts <- range(c(FSx, fsole$ssb))

### Pacific cod

# Extract Pacific cod data
pcod <- GoA.data[GoA.data$species == "Pacific cod",]
pcod <- pcod[,-2] # removes species name

# Calculate starting values for parameters
PCstart <- srStarts(recruits~ssb, data = pcod, type = "Ricker")
PCstart

# Fit model using nls
PCmodel <- nls(log.recruits~log(ricker(ssb,a,b)), data = pcod, start = PCstart)

# Calculate 95% confidence intervals for parameter estimates using the bootstrap method
PCbootR <- nlsBoot(PCmodel)
cbind(estimates = coef(PCmodel), confint(PCbootR))

# Produce values of S to predict new values of R
PCx <- seq(min(pcod$ssb), max(pcod$ssb), length.out = 199) 
PCpredR <- ricker(PCx, a = coef(PCmodel)) 
PCLCI <- PCUCI <- numeric(length(PCx))
for(i in 1:length(PCx)) { # stores a 95% confidence interval for each predicted value of R
  tmp <- apply(PCbootR$coefboot, MARGIN = 1, FUN = ricker, S = PCx[i])
  PCLCI[i] <- quantile(tmp, 0.025)
  PCUCI[i] <- quantile(tmp, 0.975)
}

# Creates axis limits for plot
PCylmts <- range(c(PCpredR, PCLCI, PCUCI, pcod$recruits))
PCxlmts <- range(c(PCx, pcod$ssb))

### Pacific ocean perch

# Extract Pacific ocean perch data
poperch <- GoA.data[GoA.data$species == "Pacific ocean perch",]
poperch <- poperch[,-2] # removes species name

# Calculate starting values for parameter estimation
POPstart <- srStarts(recruits~ssb, data = poperch, type = "Ricker")
POPstart

# Fit model using nls
POPmodel <- nls(log.recruits~log(ricker(ssb,a,b)), data = poperch, start = POPstart) 

# Calculate 95% confidence intervals for parameter estimates using the bootstrap method
POPbootR <- nlsBoot(POPmodel)
cbind(estimates = coef(POPmodel), confint(POPbootR))

# Produce values of S to predict new values of R
POPx <- seq(min(poperch$ssb), max(poperch$ssb), length.out = 199)
POPpredR <- ricker(POPx, a = coef(POPmodel))
POPLCI <- POPUCI <- numeric(length(POPx))
for(i in 1:length(POPx)) { # stores a 95% confidence interval for each predicted value of R
  tmp <- apply(POPbootR$coefboot, MARGIN = 1, FUN = ricker, S = POPx[i])
  POPLCI[i] <- quantile(tmp, 0.025)
  POPUCI[i] <- quantile(tmp, 0.975)
}

# Create axis limits for plot
POPylmts <- range(c(POPpredR, POPLCI, POPUCI, poperch$recruits))
POPxlmts <- range(c(POPx, poperch$ssb))

### Rock sole

# Extract rock sole data
rsole <- GoA.data[GoA.data$species == "Rock sole",]
rsole <- rsole[,-2] # removes species name

# Calculate starting values for parameter estimation
RSstart <- srStarts(recruits~ssb, data = rsole, type = "Ricker")
RSstart

# Fit the model using nls
RSmodel <- nls(log.recruits~log(ricker(ssb,a,b)), data = rsole, start = RSstart)

# Calculate 95% confidence intervals for parameter estimates using the bootstrap method
RSbootR <- nlsBoot(RSmodel) # 18 warnings!
nrow(RSbootR$coefboot)
cbind(estimates = coef(RSmodel), confint(RSbootR)) 

# Produce values of S to predict new values of R
RSx <- seq(min(rsole$ssb), max(rsole$ssb), length.out = 199)
RSpredR <- ricker(RSx, a = coef(RSmodel)) 
RSLCI <- RSUCI <- numeric(length(RSx))
for(i in 1:length(RSx)) { # stores a 95% confidence interval for each predicted value of R
  tmp <- apply(RSbootR$coefboot, MARGIN = 1, FUN = ricker, S = RSx[i])
  RSLCI[i] <- quantile(tmp, 0.025)
  RSUCI[i] <- quantile(tmp, 0.975)
}

# Create axis limits for plot
RSylmts <- range(c(RSpredR, RSLCI, RSUCI, rsole$recruits))
RSxlmts <- range(c(RSx, rsole$ssb))

### Walleye pollock

# Extract walleye pollock data
wpollock <- GoA.data[GoA.data$species == "Walleye pollock",]
wpollock <- wpollock[,-2] # removes species name

# Calculate starting values for parameter estimation
WPstart <- srStarts(recruits~ssb, data = wpollock, type = "Ricker")
WPstart # "Warning - b negative likely poor starting value" - nls still does not improve on it!

# Fit the model using nls
WPmodel <- nls(log.recruits~log(ricker(ssb,a,b)), data = wpollock, start = WPstart)

# Calculate 95% confidence intervals for parameter estimates using the bootstrap method
WPbootR <- nlsBoot(WPmodel) # 50 or more warnings!
nrow(WPbootR$coefboot) # Number of successful iterations
cbind(estimates = coef(WPmodel), confint(WPbootR))

# Produce values of S to predict new values of R
WPx <- seq(min(wpollock$ssb), max(wpollock$ssb), length.out = 199)
WPpredR <- ricker(WPx, a = coef(WPmodel))
WPLCI <- WPUCI <- numeric(length(WPx))
for(i in 1:length(WPx)) { # stores a 95% confidence interval for each predicted value of R
  tmp <- apply(WPbootR$coefboot, MARGIN = 1, FUN = ricker, S = WPx[i])
  WPLCI[i] <- quantile(tmp, 0.025)
  WPUCI[i] <- quantile(tmp, 0.975)
}

# Create axis limits for plot
WPylmts <- range(c(WPpredR, WPLCI, WPUCI, wpollock$recruits))
WPxlmts <- range(c(WPx, wpollock$ssb))

##### Combined plot #####

# Divide plot into 6
par(mfrow = c(3,2), mar = c(5.1, 4.5, 4.1, 2)) 

### Arrowtooth flounder

plot(recruits~ssb, data = aflounder, 
     xlim = AFxlmts, ylim = AFylmts, 
     col = "white", 
     ylab = "Recruits (in millions)", xlab = "SSB in (thousand) tonnes", cex.lab = 1.4,
     main = "Arrowtooth flounder", cex.main = 1.75,
     yaxt = "n", xaxt = "n")

# Add axis in thousands and millions
axis(1, at = pretty(aflounder$ssb), labels = label_number(scale = 1e-3)(pretty(aflounder$ssb)))
axis(2, at = pretty(aflounder$recruits), labels = label_number(scale = 1e-6, suffix = "M")(pretty(aflounder$recruits)))

# Add 95% confidence interval for predictions onto plot
polygon(c(AFx, rev(AFx)), c(AFLCI,rev(AFUCI)), col = palette.colors(7)[2], border = NA)

# Add existing data points and Ricker curve
points(recruits~ssb, data = aflounder, pch = 19, col = rgb(0,0,0,1/2))
lines(AFpredR~AFx, lwd = 2)

### Flathead sole

plot(recruits~ssb, data = fsole, 
     xlim = FSxlmts, ylim = FSylmts, 
     col = "white", 
     ylab = "Recruits (in millions)", xlab = "SSB in (thousand) tonnes", cex.lab = 1.4,
     main = "Flathead sole", cex.main = 1.75,
     yaxt = "n", xaxt = "n")

# Add axis in thousands and millions
axis(1, at = pretty(fsole$ssb), labels = label_number(scale = 1e-3)(pretty(fsole$ssb)))
axis(2, at = pretty(fsole$recruits), labels = label_number(scale = 1e-6)(pretty(fsole$recruits)))

# Add 95% confidence interval for predictions onto plot
polygon(c(FSx, rev(FSx)), c(FSLCI,rev(FSUCI)), col = palette.colors(7)[3], border = NA)

# Add existing data points and Ricker curve
points(recruits~ssb, data = fsole, pch = 19, col = rgb(0,0,0,1/2))
lines(FSpredR~FSx, lwd = 2)

### Pacific cod

plot(recruits~ssb, data = pcod, 
     xlim = PCxlmts, ylim = PCylmts, 
     col = "white", 
     ylab = "Recruits (in millions)", xlab = "SSB in (thousand) tonnes", cex.lab = 1.4,
     main = "Pacific cod", cex.main = 1.75,
     yaxt = "n", xaxt = "n")

# Add axis in thousands and millions
axis(1, at = pretty(pcod$ssb), labels = label_number(scale = 1e-3)(pretty(pcod$ssb)))
axis(2, at = pretty(pcod$recruits), labels = label_number(scale = 1e-6)(pretty(pcod$recruits)))

# Add 95% confidence interval for predictions onto plot
polygon(c(PCx, rev(PCx)), c(PCLCI,rev(PCUCI)), col = palette.colors(7)[4], border = NA)

# Add existing data points and Ricker curve
points(recruits~ssb, data = pcod, pch = 19, col = rgb(0,0,0,1/2))
lines(PCpredR~PCx, lwd = 2)

### Pacific ocean perch

plot(recruits~ssb, data = poperch, 
     xlim = POPxlmts, ylim = POPylmts, 
     col = "white", 
     ylab = "Recruits (in millions)", xlab = "SSB in (thousand) tonnes", cex.lab = 1.4,
     main = "Pacific ocean perch", cex.main = 1.75,
     yaxt = "n", xaxt = "n")

# Add axis in thousands and millions
axis(1, at = pretty(poperch$ssb), labels = label_number(scale = 1e-3)(pretty(poperch$ssb)))
axis(2, at = pretty(poperch$recruits), labels = label_number(scale = 1e-6)(pretty(poperch$recruits)))

# Add 95% confidence intervals for predictions onto plot
polygon(c(POPx, rev(POPx)), c(POPLCI,rev(POPUCI)), col = palette.colors(7)[5], border = NA)

# Add existing data points and Ricker curve
points(recruits~ssb, data = poperch, pch = 19, col = rgb(0,0,0,1/2))
lines(POPpredR~POPx, lwd = 2)

### Rock sole

plot(recruits~ssb, data = rsole, 
     xlim = RSxlmts, ylim = RSylmts, 
     col = "white", 
     ylab = "Recruits (in millions)", xlab = "SSB in (thousand) tonnes", cex.lab = 1.4,
     main = "Rock sole", cex.main = 1.75,
     yaxt = "n", xaxt = "n")

# Add axis in thousands and millions
axis(1, at = pretty(rsole$ssb), labels = label_number(scale = 1e-3)(pretty(rsole$ssb)))
axis(2, at = pretty(rsole$recruits), labels = label_number(scale = 1e-6)(pretty(rsole$recruits)))

# Add 95% confidence interval for predictions onto plot
polygon(c(RSx, rev(RSx)), c(RSLCI,rev(RSUCI)), col = palette.colors(7)[6], border = NA)

# Add existing data points and Ricker curve
points(recruits~ssb, data = rsole, pch = 19, col = rgb(0,0,0,1/2))
lines(RSpredR~RSx, lwd = 2)

### Walleye pollock

plot(recruits~ssb, data = wpollock, 
     xlim = WPxlmts, ylim = WPylmts, 
     col = "white", 
     ylab = "Recruits (in millions)", xlab = "SSB in (thousand) tonnes", cex.lab = 1.4,
     main = "Walleye pollock", cex.main = 1.75,
     yaxt = "n", xaxt = "n")

# Add axis in thousands and millions
axis(1, at = pretty(wpollock$ssb), labels = label_number(scale = 1e-3)(pretty(wpollock$ssb)))
axis(2, at = pretty(wpollock$recruits), labels = label_number(scale = 1e-6)(pretty(wpollock$recruits)))

# Add 95% confidence interval for predictions onto plot
polygon(c(WPx, rev(WPx)), c(WPLCI,rev(WPUCI)), col = palette.colors(7)[7], border = NA)

# Add existing data points and Ricker curve
points(recruits~ssb, data = wpollock, pch = 19, col = rgb(0,0,0,1/2))
lines(WPpredR~WPx, lwd = 2)

