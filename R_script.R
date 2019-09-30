# Script used in the research paper "Modelling weight growth of Santa
# Ines sheep using Bayesian approach". Any doubts or comments should be
# addressed to thiagotsalles@gmail.com

#This script requires OpenBUGS to be already installed on the computer

# Load libraries
library(R2OpenBUGS)
library(coda)
library(nortest)

# Set working directory
setwd("C:/my_dir")

# Sample (Teixeira Neto et al., 2016)
x <- c(0, 30, 60, 90, 120, 150, 180, 210)
y <- c(3.4, 10.1, 17.6, 21.3, 23.8, 27, 31.7, 33.6)

# Prior information about weight
x.sarm <- c(0, 28, 56, 84, 112, 140, 168, 196) # SARMENTO et al. (2006)
y.sarm <- c(3.5, 8.3, 12.7, 16.4, 19.5, 21.4, 22.2, 23.2)

x.malhad <- c(0, 60, 120, 180, 240, 300, 365) # MALHADO et al. (2008)
y.malhad <- c(3.42, 11.98, 16.21, 21.01, 24.04, 24.65, 31.09)

x.car <- c(0, 60, 90, 120, 180, 365) # CARNEIRO et al. (2010)
y.car <- c(3.72, 11.61, 15.89, 18.73, 28.53, 37.05)

x.o <- c(0, 28, 56, 112, 350) # Ó et al. (2012)
y.o <- c(3.38, 6.3, 9.4, 13.5, 25.79)

x.teix <- c(0, 60, 120, 180, 240, 300, 360) # TEIXEIRA et al. (2012)
y.teix <- c(3.00, 11.73, 17.06, 23.35, 24.85, 26.97, 27.78)

x.santos <- c(0, 56, 112, 168, 224, 280, 336) # SANTOS et al. (2014)
y.santos <- c(3.37, 11.6, 17.62, 21.51, 24.01, 25.72, 26.98)

# Items required to fit the model with OpenBUGS software (Bayesian
# approach). They are the sample size, the names of the variables, the
# names of the parameters and the initial values for the parameters.
n <- length(x)
data <- list("x", "y", "n")
theta <- c("alpha", "beta", "gamma", "tau2")
inits.gomp <- function() {list(alpha=30.0, beta=0.9, gamma=0.01,
                               tau2=1.0)}

# Function for calculating BIC
BIC_func = function(y, res, df=4) {
  n = length(y)
  w = rep(1, n)
  ll = 0.5 *
    (sum(log(w)) - n * (log(2 * pi) + 1 - log(n) + log(sum(w * res^2))))
  bic <- -2 * ll + log(n) * df
}

# Function for calculating RMSE
RMSE = function(res) {
  n = length(res)
  sq_err = 0
  for (i in 1:n) {
    sq_err = sq_err + res[i] ^2
  }
  rmse = sqrt(sq_err / n)
}

#### ==================== FREQUENTIST APPROACH ==================== ####
# Model fit (nonlinear least sqares)
gomp <- nls(y ~ (alpha * exp(-beta * exp(-gamma * x))),
            start=c(alpha=30.0, beta=0.9, gamma=0.01))

# Results table
results.freq <- matrix(c(coefficients(gomp), confint(gomp)),
                       nrow=3, ncol=3)

colnames(results.freq) <- c("Estimate",	"CI lower limit",
                            "CI upper limit")

rownames(results.freq) <- c("alpha", "beta", "gamma")

# Adjusted R²
r2 <- 1 - (summary(gomp)[[3]] ^ 2 * (length(y) - 3) /
             sum((y - mean(y)) ^ 2))

r2.adj <- 1 - ((length(y) - 1) / (length(y) - (3 + 1))) * (1 - r2)

# Residual normality tests
shapiro.test(summary(gomp)$resid)
lillie.test(summary(gomp)$resid)

# Table with estimated weights and confidence intervals
est.y <- fitted(gomp)
IC.upp <- (est.y + qt(0.975, summary(gomp)$df[2]) *
             summary(gomp)$sigma)

IC.low <- (est.y - qt(0.975, summary(gomp)$df[2]) *
             summary(gomp)$sigma)

est.freq <- matrix(c(est.y, IC.low, IC.upp), nrow=8, ncol=3)
colnames(est.freq) <- c("Estimate",	"CI lower limit",	"CI upper limit")
rownames(est.freq) <- paste(x)

# RMSE and BIC
res.f <- summary(gomp)$resid
rmse.f <- RMSE(res.f)
bic.f <- BIC_func(y, res.f)

#### ====================== PRIOR INFORMATION ===================== ####
# Model fit for obtaining prior information about the parameters
gomp.sarm <- nls(y.sarm ~ (
  alpha * exp(-beta * exp(-gamma * x.sarm))),
  start=c(alpha=30.0, beta=0.9, gamma=0.01))

gomp.malhad <- nls(y.malhad ~ (
  alpha * exp(-beta * exp(-gamma * x.malhad))),
  start=c(alpha=30.0, beta=0.9, gamma=0.01))

gomp.car <- nls(y.car ~ (
  alpha * exp(-beta * exp(-gamma * x.car))),
  start=c(alpha=30.0, beta=0.9, gamma=0.01))

gomp.o <- nls(y.o ~ (
  alpha * exp(-beta * exp(-gamma * x.o))),
  start=c(alpha=30.0, beta=0.9, gamma=0.01))

gomp.teix <- nls(y.teix ~ (
  alpha * exp(-beta * exp(-gamma * x.teix))),
  start=c(alpha=30.0, beta=0.9, gamma=0.01))

gomp.santos <- nls(y.santos ~ (
  alpha * exp(-beta * exp(-gamma * x.santos))),
  start=c(alpha=30.0, beta=0.9, gamma=0.01))

# Results table
prior.params <- matrix(
  c(coefficients(gomp.sarm), summary(gomp.sarm)[[3]],
    coefficients(gomp.malhad), summary(gomp.malhad)[[3]],
    coefficients(gomp.car), summary(gomp.car)[[3]],
    coefficients(gomp.o), summary(gomp.o)[[3]],
    coefficients(gomp.teix), summary(gomp.teix)[[3]],
    coefficients(gomp.santos), summary(gomp.santos)[[3]]),
  nrow=4, ncol=6)

colnames(prior.params) <- c("SARMENTO et al. (2006)",
                            "MALHADO et al. (2008)",
                            "CARNEIRO et al. (2010)",
                            "Ó et al. (2012)",
                            "TEIXEIRA et al. (2012)",
                            "SANTOS et al. (2014)")

rownames(prior.params) <- c("alpha", "beta", "gamma", "sigma")


#### ========= BAYESIAN APPROACH - NONINFORMATIVE PRIORS ========== ####
# Creating a text file to run the "bugs" function in OpenBUGS software
sink("gomp_0.txt")
cat("
    model{
    
    # PRIOR DISTRIBUTIONS
    alpha ~ dnorm(0.0, 0.0001)
    beta ~ dnorm(0.0, 0.0001)
    gamma ~ dunif(0.0, 1.0000)
    tau2 ~ dgamma(0.1, 0.1)
    
    # MODEL
    for(i in 1 : n) {
    y[i] ~ dnorm(mu[i], tau2)
    mu[i] <- alpha * exp(-beta * exp(-gamma * x[i]))
    }
    
    # ESTIMATES FOR WEIGHT
    for(i in 1 : n) {
    yp[i] <- alpha * exp(-beta * exp(-gamma * x[i]))
    }
    
    }", fill=TRUE)
sink()

# Model fit. OpenBUGS will open to run the "bugs" function.
# OpenBUGS has to be closed after it finishes running its algorithm so
# the results can return to R.
gomp.0 = bugs(data, inits.gomp,
              model.file="gomp_0.txt",
              parameters=c(theta, "yp"),
              n.chains=1, n.iter=25000, n.burnin=4500, n.thin=5,
              codaPkg=FALSE, debug=T)

## Convergence tests ##
# Save bugs object as MCMC object
gomp.0.mcmc <- as.mcmc(gomp.0$sims.matrix[, ])
# Tests
geweke.diag(gomp.0.mcmc)
raftery.diag(gomp.0.mcmc)
heidel.diag(gomp.0.mcmc)

## Results table, with parameters and weight estimates ##
hpd.interval.0 <- HPDinterval(gomp.0.mcmc)
results.bayes.0 <- matrix(c(gomp.0$summary[, 1], hpd.interval.0),
                          nrow=14, ncol=3)

colnames(results.bayes.0) <- c("Estimate", "HPD lower limit",
                               "HPD upper limit")

rownames(results.bayes.0) <- c(rownames(hpd.interval.0))

# RMSE and BIC
res.0 <- unname(y - results.bayes.0[, 1][5:12])
rmse.0 <- RMSE(res.0)
bic.0 <- BIC_func(y, res.0)

#### ========== BAYESIAN APPROACH - INFORMATIVE PRIORS ============ ####
# Creating a text file to run the "bugs" function in OpenBUGS software
sink("gomp_1.txt")
cat("
    model{
    
    # PRIOR DISTRIBUTIONS
    alpha ~ dnorm(29.6760, 0.0360)
    beta ~ dnorm(1.9966, 26.2314)
    gamma ~ dunif(0.0081, 0.0193)
    tau2 ~ dgamma(19.9023, 21.7393)
    
    # MODEL
    for(i in 1 : n) {
    y[i] ~ dnorm(mu[i], tau2)
    mu[i] <- alpha * exp(-beta * exp(-gamma * x[i]))
    }
    
    # ESTIMATES FOR WEIGHT
    for(i in 1 : n) {
    yp[i] <- alpha * exp(-beta * exp(-gamma * x[i]))
    }
    
    }", fill=TRUE)
sink()

# Model fit. OpenBUGS will open to run the "bugs" function.
# OpenBUGS has to be closed after it finishes running its algorithm so
# the results can return to R.
gomp.1 = bugs(data, inits.gomp,
              model.file="gomp_1.txt",
              parameters=c(theta, "yp"),
              n.chains=1, n.iter=25000, n.burnin=4500, n.thin=5,
              codaPkg=FALSE, debug=T)

## Convergence tests ##
# Save bugs object as MCMC object
gomp.1.mcmc <- as.mcmc(gomp.1$sims.matrix[, ])
# Tests
geweke.diag(gomp.1.mcmc)
raftery.diag(gomp.1.mcmc)
heidel.diag(gomp.1.mcmc)

## Results table ##
hpd.interval.1 <- HPDinterval(gomp.1.mcmc)
results.bayes.1 <- matrix(c(gomp.1$summary[, 1], hpd.interval.1),
                          nrow=14, ncol=3)

colnames(results.bayes.1) <- c("Estimate", "HPD lower limit",
                               "HPD upper limit")

rownames(results.bayes.1) <- c(rownames(hpd.interval.1))

# RMSE and BIC
res.1 <- unname(y - results.bayes.1[, 1][5:12])
rmse.1 <- RMSE(res.1)
bic.1 <- BIC_func(y, res.1)


#### ==== BAYESIAN APPROACH - INF. PRIORS (GREATER DISPERSION) ==== ####
# Creating a text file to run the "bugs" function in OpenBUGS software
sink("gomp_2.txt")
cat("
    model{
    
    # PRIOR DISTRIBUTIONS
    alpha ~ dnorm(29.6760, 0.0090)
    beta ~ dnorm(1.9966, 6.5578)
    gamma ~ dunif(0.0000, 0.0361)
    tau2 ~ dgamma(6.3176, 6.9007)
    
    # MODEL
    for(i in 1 : n) {
    y[i] ~ dnorm(mu[i], tau2)
    mu[i] <- alpha * exp(-beta * exp(-gamma * x[i]))
    }
    
    # ESTIMATES FOR WEIGHT
    for(i in 1 : n) {
    yp[i] <- alpha * exp(-beta * exp(-gamma * x[i]))
    }
    
    }", fill=TRUE)

sink()

# Model fit. OpenBUGS will open to run the "bugs" function.
# OpenBUGS has to be closed after it finishes running its algorithm so
# the results can return to R.
gomp.2 = bugs(data, inits.gomp,
              model.file="gomp_2.txt",
              parameters=c(theta, "yp"),
              n.chains=1, n.iter=25000, n.burnin=4500, n.thin=5,
              codaPkg=FALSE, debug=T)

## Convergence tests ##
# Save bugs object as MCMC object
gomp.2.mcmc <- as.mcmc(gomp.2$sims.matrix[, ])
# Tests
geweke.diag(gomp.2.mcmc)
raftery.diag(gomp.2.mcmc)
heidel.diag(gomp.2.mcmc)

## Results table ##
hpd.interval.2 <- HPDinterval(gomp.2.mcmc)
results.bayes.2 <- matrix(c(gomp.2$summary[, 1], hpd.interval.2),
                          nrow=14, ncol=3)

colnames(results.bayes.2) <- c("Estimate", "HPD lower limit",
                               "HPD upper limit")

rownames(results.bayes.2) <- c(rownames(hpd.interval.2))

# RMSE and BIC
res.2 <- unname(y - results.bayes.2[, 1][5:12])
rmse.2 <- RMSE(res.2)
bic.2 <- BIC_func(y, res.2)

