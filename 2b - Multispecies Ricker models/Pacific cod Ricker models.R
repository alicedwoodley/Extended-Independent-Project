##### Pacific cod Ricker models #####

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
pcod <- GoA.data[GoA.data$species == "Pacific cod",]
pcod <- pcod[,-2] # removes species name

# Calculate starting values for parameter estimation
PCstart <- srStarts(recruits~ssb, data = pcod, type = "Ricker")
PCstart

# Fit model using nls
PCmodel <- nls(log.recruits~log(ricker(ssb,a,b)), data = pcod, start = PCstart) 

##### Multispecies model #####

# Add other species SSB data to pcod
pcod$afssb <- GoA.data[GoA.data$species == "Arrowtooth flounder", 3]
pcod$fsssb <- GoA.data[GoA.data$species == "Flathead sole", 3]
pcod$rsssb <- GoA.data[GoA.data$species == "Rock sole", 3]
pcod$wpssb <- GoA.data[GoA.data$species == "Walleye pollock", 3]

# Create model function
PCricker <- function(S, X1, X2, X3, X4, a, b = NULL, c = NULL, d = NULL, f = NULL, g = NULL){ # S = pc , X1 = af, X2 = fs, X3 = rs, X4 = wp
  if(length(a)>1){
    g <- a[6]
    f <- a[5]
    d <- a[4]
    c <- a[3]
    b <- a[2]
    a <- a[1]
  }
  a*S*exp(-b*S+c*X1+d*X2+f*X3+g*X4) 
}

# Find starting values through linearising and using lm()
linear <- lm(log(recruits/ssb)~ ssb + afssb + fsssb + rsssb + wpssb, data = pcod)
PCstart2 <- coef(linear) # extract coefficients
PCstart2 <- list(a = exp(PCstart2[[1]]), b = PCstart2[[2]], c = PCstart2[[3]], d = PCstart2[[4]], f = PCstart2[[5]], g = PCstart2[[6]]) # prep to put in model

# Fit model using nls
PCmultimodel <- nls(log.recruits~log(PCricker(ssb, afssb, fsssb, rsssb, wpssb, a,b,c,d,f,g)), data = pcod, start = PCstart2) 
summary(PCmultimodel)

# Calculate 95% confidence intervals for parameter estimates using the bootstrap method
PCbootR2 <- nlsBoot(PCmultimodel) # 50 or more warnings
nrow(PCbootR2$coefboot) # 934 successful iterations
cbind(estimates = coef(PCmultimodel), confint(PCbootR2))

# Produce values of S to predict new values of R
PCx2 <- seq(min(pcod$ssb), max(pcod$ssb), length.out = 38)
PCpredR2 <- PCricker(PCx2, pcod$afssb, pcod$fsssb, pcod$rsssb, pcod$wpssb, a = coef(PCmultimodel))
PCLCI2 <- PCUCI2 <- numeric(length(PCx2))
for(i in 1:length(PCx2)) { # stores a 95% confidence interval for each predicted value of R
  tmp <- apply(PCbootR2$coefboot, MARGIN = 1, function(a_tmp) 
    PCricker(S = PCx2[i], X1 = pcod$afssb[i], X2 = pcod$fsssb[i], X3 = pcod$rsssb[i], X4 = pcod$wpssb[i], a = a_tmp)
  )
  PCLCI2[i] <- quantile(tmp, 0.025)
  PCUCI2[i] <- quantile(tmp, 0.975)
}

# Create axis limits for plot
PCylmts2 <- range(c(PCpredR2, PCLCI2, PCUCI2, pcod$recruits))
PCxlmts2 <- range(c(PCx2, pcod$ssb))

### Plot:

plot(recruits~ssb, data = pcod, 
     xlim = PCxlmts2, ylim = PCylmts2, 
     col = "white", 
     ylab = "Recruits (in millions)", xlab = "SSB in (thousand) tonnes",
     main = "Pacific cod multispecies Ricker model",
     yaxt = "n", xaxt = "n")

# Add axis in thousands and millions
axis(1, at = pretty(pcod$ssb), labels = label_number(scale = 1e-3)(pretty(pcod$ssb)))
axis(2, at = pretty(PCUCI2), labels = label_number(scale = 1e-6)(pretty(PCUCI2)))

# Add 95% confidence intervals for predictions onto plot
polygon(c(PCx2, rev(PCx2)), c(PCLCI2,rev(PCUCI2)), col = palette.colors(7)[4], border = NA)

# Add existing data points and Ricker curve
points(recruits~ssb, data = pcod, pch = 19, col = rgb(0,0,0,1/2))
lines(PCpredR2~PCx2, lwd = 2)

##### Comparison with single species model #####

cbind("Single species" = AIC(PCmodel), "Multispecies" = AIC(PCmultimodel))

# Verify residuals are normal so assumptions in AIC hold

par(mfrow = c(1,2))

qqnorm(resid(PCmodel), main = "Q-Q plot for Pacific cod Ricker model", cex.main = 0.9)
qqline(resid(PCmodel), col = palette.colors(7)[4], lwd = 1.5)

qqnorm(resid(PCmultimodel), main = "Q-Q plot for Pacific cod multispecies model", cex.main = 0.9)
qqline(resid(PCmultimodel), col = palette.colors(7)[4], lwd = 1.5)
