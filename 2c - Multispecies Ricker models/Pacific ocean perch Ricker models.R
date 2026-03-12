##### Pacific ocean perch single and multispecies Ricker models #####

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

##### Single species model for comparison #####

# Extract Pacific ocean perch data
poperch <- GoA.data[GoA.data$species == "Pacific ocean perch",]
poperch <- poperch[,-2] # removes species name

# Calculate starting values for parameter estimation
POPstart <- srStarts(recruits~ssb, data = poperch, type = "Ricker")
POPstart

# Fit model using nls
POPmodel <- nls(log.recruits~log(ricker(ssb,a,b)), data = poperch, start = POPstart) 

##### Multispecies model #####

# Add walleye pollock SSB data to poperch
poperch$wpssb <- GoA.data[GoA.data$species == "Walleye pollock", 3]

# Create model function
POPricker <- function(S, X, a, b = NULL, c = NULL){ # S = POP SSB, X = walleye pollock SSB
  if(length(a)>1){
    c <- a[3]
    b <- a[2]
    a <- a[1]
  }
  a*S*exp(-b*S+c*X) 
}

# Find starting values through linearising and using lm()
linear <- lm(log(recruits/ssb)~ ssb + wpssb, data = poperch)
POPstart2 <- coef(linear) # extract coefficients
POPstart2 <- list(a = exp(POPstart2[[1]]), b = POPstart2[[2]], c = POPstart2[[3]]) # prep to put in model

# Fit model using nls
POPmultimodel <- nls(log.recruits~log(POPricker(ssb, wpssb, a,b,c)), data = poperch, start = POPstart2) 
summary(POPmultimodel)

# Calculate 95% confidence intervals for parameter estimates using the bootstrap method
POPbootR2 <- nlsBoot(POPmultimodel) # Warning - fit did not converge 10 times during bootstrapping
cbind(estimates = coef(POPmultimodel), confint(POPbootR2))

# Produce values of S to predict new values of R
POPx2 <- seq(min(poperch$ssb), max(poperch$ssb), length.out = 38)
POPpredR2 <- POPricker(POPx2, poperch$wpssb, a = coef(POPmultimodel))
POPLCI2 <- POPUCI2 <- numeric(length(POPx2))
for(i in 1:length(POPx2)) { # stores a 95% confidence interval for each predicted value of R
  tmp <- apply(POPbootR2$coefboot, MARGIN = 1, function(a_tmp) 
    POPricker(S = POPx2[i], X = poperch$wpssb[i], a = a_tmp)
  )
  POPLCI2[i] <- quantile(tmp, 0.025)
  POPUCI2[i] <- quantile(tmp, 0.975)
}

# Create axis limits for plot
POPylmts2 <- range(c(POPpredR2, POPLCI2, POPUCI2, poperch$recruits))
POPxlmts2 <- range(c(POPx2, poperch$ssb))

### Plot:

plot(recruits~ssb, data = poperch, 
     xlim = POPxlmts2, ylim = POPylmts2, 
     col = "white", 
     ylab = "Recruits (in millions)", xlab = "SSB in (thousand) tonnes",
     main = "Pacific ocean perch multispecies Ricker model",
     yaxt = "n", xaxt = "n")

# Add axis in thousands and millions
axis(1, at = pretty(poperch$ssb), labels = label_number(scale = 1e-3)(pretty(poperch$ssb)))
axis(2, at = pretty(poperch$recruits), labels = label_number(scale = 1e-6)(pretty(poperch$recruits)))

# Add 95% confidence intervals for predictions onto plot
polygon(c(POPx2, rev(POPx2)), c(POPLCI2,rev(POPUCI2)), col = palette.colors(7)[5], border = NA)

# Add existing data points and Ricker curve
points(recruits~ssb, data = poperch, pch = 19, col = rgb(0,0,0,1/2))
lines(POPpredR2~POPx2, lwd = 2)

##### Comparison with single species model #####

cbind("Single species" = AIC(POPmodel), "Multispecies" = AIC(POPmultimodel))
