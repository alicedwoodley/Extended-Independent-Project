##### Walleye pollock single species and multispecies Ricker models #####

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

# Extract walleye pollock data
wpollock <- GoA.data[GoA.data$species == "Walleye pollock",]
wpollock <- wpollock[,-2] # removes species name

##### Single species model #####

# Calculate starting values for parameter estimation
WPstart <- srStarts(recruits~ssb, data = wpollock, type = "Ricker")
WPstart # "Warning - b negative likely poor starting value" - nls still does not improve on it!

# Fit the model using nls
WPmodel <- nls(log.recruits~log(ricker(ssb,a,b)), data = wpollock, start = WPstart)

##### What happens to single species model if the crazy 2012 point is excluded? #####

wpollock_test <- wpollock[-which(wpollock$year=="2013"),]

# Calculate starting values for parameter estimation
WPstart_test <- srStarts(recruits~ssb, data = wpollock_test, type = "Ricker")
WPstart_test # "Warning - b negative likely poor starting value" - nls still does not improve on it!

# Fit the model using nls
WPmodel_test <- nls(log.recruits~log(ricker(ssb,a,b)), data = wpollock_test, start = WPstart_test)

# Calculate 95% confidence intervals for parameter estimates using the bootstrap method
WPbootR_test <- nlsBoot(WPmodel_test) # 50 or more warnings!
cbind(estimates = coef(WPmodel_test), confint(WPbootR_test))

# Produce values of S to predict new values of R
WPx_test <- seq(min(wpollock_test$ssb), max(wpollock_test$ssb), length.out = 199)
WPpredR_test <- ricker(WPx_test, a = coef(WPmodel_test))
WPLCI_test <- WPUCI_test <- numeric(length(WPx_test))
for(i in 1:length(WPx_test)) { # stores a 95% confidence interval for each predicted value of R
  tmp <- apply(WPbootR_test$coefboot, MARGIN = 1, FUN = ricker, S = WPx_test[i])
  WPLCI_test[i] <- quantile(tmp, 0.025)
  WPUCI_test[i] <- quantile(tmp, 0.975)
}

# Create axis limits for plot
WPylmts_test <- range(c(WPpredR_test, WPLCI_test, WPUCI_test, wpollock_test$recruits))
WPxlmts_test <- range(c(WPx_test, wpollock_test$ssb))

# Plot:

par(mfrow = c(1,1))
plot(recruits~ssb, data = wpollock_test, 
     xlim = WPxlmts_test, ylim = WPylmts_test, 
     col = "white", 
     ylab = "Recruits (in millions)", xlab = "SSB in (thousand) tonnes",
     main = "Walleye pollock (with unusual year removed)",
     yaxt = "n", xaxt = "n")

# Add axis in thousands and millions
axis(1, at = pretty(wpollock_test$ssb), labels = label_number(scale = 1e-3)(pretty(wpollock_test$ssb)))
axis(2, at = pretty(wpollock_test$recruits), labels = label_number(scale = 1e-6)(pretty(wpollock_test$recruits)))

# Add 95% confidence interval for predictions onto plot
polygon(c(WPx_test, rev(WPx_test)), c(WPLCI_test,rev(WPUCI_test)), col = palette.colors(7)[7], border = NA)

# Add existing data points and Ricker curve
points(recruits~ssb, data = wpollock_test, pch = 19, col = rgb(0,0,0,1/2))
lines(WPpredR_test~WPx_test, lwd = 2)


##### Multispecies model #####

# Add other species to wpollock
wpollock$afssb <- GoA.data[GoA.data$species == "Arrowtooth flounder", 3]
wpollock$fsssb <- GoA.data[GoA.data$species == "Flathead sole", 3]
wpollock$pcssb <- GoA.data[GoA.data$species == "Pacific cod", 3]
wpollock$popssb <- GoA.data[GoA.data$species == "Pacific ocean perch", 3]
wpollock$rsssb <- GoA.data[GoA.data$species == "Rock sole", 3]

WPricker1 <- function(S, X1, X2, X3, X4, X5, a, b = NULL, c = NULL, d = NULL, f = NULL, g = NULL, h  = NULL){ # S = wp stock, X1 = af, X2 = fs, X3 = pc, X4 = pop, X5 = rs
  if(length(a)>1){
    h <- a[7]
    g <- a[6]
    f <- a[5]
    d <- a[4]
    c <- a[3]
    b <- a[2]
    a <- a[1]
  }
  a*S*exp(-b*S+c*X1+d*X2+f*X3+g*X4+h*X5)
}

# Find starting values through linearising and using lm()
linear <- lm(log(recruits/ssb)~ ssb + afssb + fsssb + pcssb + popssb + rsssb, data = wpollock)
WPstart2 <- coef(linear) # extract coefficients
WPstart2 <- list(a = exp(WPstart2[[1]]), b = WPstart2[[2]], c = WPstart2[[3]], d = WPstart2[[4]], f = WPstart2[[5]], g = WPstart2[[6]], h = WPstart2[[7]]) # prep to put in model

# Fit model using nls
WPmultimodel <- nls(log.recruits~log(WPricker1(ssb, afssb, fsssb, pcssb, popssb, rsssb, a,b,c,d,f,g,h)), data = wpollock, start = WPstart2) 
summary(WPmultimodel)

# Calculate 95% confidence intervals for parameter estimates using the bootstrap method
WPbootR2 <- nlsBoot(WPmultimodel) # 50 or more warnings again - expected when initial parameter estimates are negative
nrow(WPbootR2$coefboot) # 623 successful iterations
cbind(estimates = coef(WPmultimodel), confint(WPbootR2))

# Produce values of S to predict new values of R
WPx2 <- seq(min(wpollock$ssb), max(wpollock$ssb), length.out = 38)
WPpredR2 <- WPricker1(WPx2, wpollock$afssb, wpollock$fsssb, wpollock$pcssb, wpollock$popssb, wpollock$rsssb, a = coef(WPmultimodel))
WPLCI2 <- WPUCI2 <- numeric(length(WPx2))
for(i in 1:length(WPx2)) { # stores a 95% confidence interval for each predicted value of R
  tmp <- apply(WPbootR2$coefboot, MARGIN = 1, function(a_tmp) 
    WPricker1(S = WPx2[i], X1 = wpollock$afssb[i], X2 = wpollock$fsssb[i], X3 = wpollock$pcssb[i], X4 = wpollock$popssb[i], X5 = wpollock$rsssb[i], a = a_tmp)
  )
  WPLCI2[i] <- quantile(tmp, 0.025)
  WPUCI2[i] <- quantile(tmp, 0.975)
}

# Create axis limits for plot
WPylmts2 <- range(c(WPpredR2, WPLCI2, WPUCI2, wpollock$recruits))
WPxlmts2 <- range(c(WPx2, wpollock$ssb))

### Plot:

plot(recruits~ssb, data = wpollock, 
     xlim = WPxlmts2, ylim = WPylmts2, 
     col = "white", 
     ylab = "Recruits (in millions)", xlab = "SSB in (thousand) tonnes",
     main = "Walleye pollock multispecies Ricker model",
     yaxt = "n", xaxt = "n")

# Add axis in thousands and millions
axis(1, at = pretty(wpollock$ssb), labels = label_number(scale = 1e-3)(pretty(wpollock$ssb)))
axis(2, at = pretty(WPUCI2), labels = label_number(scale = 1e-6)(pretty(WPUCI2)))

# Add 95% confidence intervals for predictions onto plot
polygon(c(WPx2, rev(WPx2)), c(WPLCI2,rev(WPUCI2)), col = palette.colors(7)[7], border = NA)

# Add existing data points and Ricker curve
points(recruits~ssb, data = wpollock, pch = 19, col = rgb(0,0,0,1/2))
lines(WPpredR2~WPx2, lwd = 2)

##### Comparison with single species model #####

cbind("Single species" = AIC(WPmodel), "Multispecies" = AIC(WPmultimodel))

# Verify residuals are normal so assumptions in AIC hold

par(mfrow = c(1,2))

qqnorm(resid(WPmodel), main = "Q-Q plot for walleye pollock Ricker model", cex.main = 0.9)
qqline(resid(WPmodel), col = palette.colors(7)[7], lwd = 1.5)

qqnorm(resid(WPmultimodel), main = "Q-Q plot for walleye pollock multispecies model", cex.main = 0.9)
qqline(resid(WPmultimodel), col = palette.colors(7)[7], lwd = 1.5)
