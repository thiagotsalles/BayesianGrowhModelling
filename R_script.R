
# Load libraries
library(R2OpenBUGS)
library(coda)
library(nortest)

# Sample
sampl <- read.csv("sample.csv", header=T, sep=";")
x <- sampl$Age
y <- sampl$Weight

# Prior information about weight
prior_info <- read.csv("prior_info.csv", header=T, sep=";")

# Items required to fit the model with OpenBUGS software (Bayesian
# approach). They are the sample size, the names of the variables, the
# names of the parameters and the initial values for the parameters.
n <- length(x)
data <- list("x", "y", "n")
theta <- c("alpha", "beta", "gamma", "tau2")
inits.gomp <- function() {list(alpha=30.0, beta=0.9, gamma=0.01,
                               tau2=1.0)}

# Function for calculating BIC
# res = residuals and df = degrees of freedom
BIC_func = function(y, res, df=4) {
  n = length(y)
  w = rep(1, n)
  ll = 0.5 *
    (sum(log(w)) - n * (log(2 * pi) + 1 - log(n) + log(sum(w * res^2))))
  bic = -2 * ll + log(n) * df
}

# Function for calculating RMSE
# res = residuals
RMSE = function(res) {
  n = length(res)
  sq_err = 0
  for (i in 1:n) {
    sq_err = sq_err + res[i] ^2
  }
  rmse = sqrt(sq_err / n)
}

# Function for calculating adjusted R²
# k = numer of regressors, rse = residual standard error
r2Adj = function(y, rse, k=3) {
  r2 = 1 - (rse ^ 2 * (length(y) - 3) / sum((y - mean(y)) ^ 2))
  r2.adj = 1 - ((length(y) - 1) / (length(y) - (k + 1))) * (1 - r2)
}

#### ==================== FREQUENTIST APPROACH ==================== ####
# Model fit (Gompertz, nonlinear least sqares)
gomp <- nls(y ~ (alpha * exp(-beta * exp(-gamma * x))),
            start=c(alpha=30.0, beta=0.9, gamma=0.01))

# Results table
results.freq <- data.frame("Estimate"=coefficients(gomp),
                           "CI_lower_limit"=confint(gomp)[, 1],
                           "CI_upper_limit"=confint(gomp)[, 2])

# Adjusted R²
r2.adj = r2Adj(y, summary(gomp)[[3]])

# Normality tests for residuals
shapiro.test(summary(gomp)$resid)
lillie.test(summary(gomp)$resid)

# Table with estimated weights and confidence intervals
est.y <- as.vector(fitted(gomp))
est.freq <- data.frame(
  "Age"=x, "Weight"=est.y,
  "CI_lower_limit"=(est.y - qt(0.975, summary(gomp)$df[2]) *
                      summary(gomp)$sigma),
  "CI_upper_limit"=(est.y + qt(0.975, summary(gomp)$df[2]) *
                      summary(gomp)$sigma))

# RMSE and BIC
res.freq <- summary(gomp)$resid
rmse.freq <- RMSE(res.freq)
bic.freq <- BIC_func(y, res.freq)

#### ====================== PRIOR INFORMATION ===================== ####
# Model fit for obtaining prior information about the parameters
# and its results table
prior.params = data.frame(row.names=c("alpha", "beta",
                                      "gamma", "sigma"))

for (i in levels(prior_info$Source)) {
  x.pr = prior_info[prior_info$Source == i, ]$Age
  y.pr = prior_info[prior_info$Source == i, ]$Weight
  gomp.pr = nls(y.pr ~ (alpha * exp(-beta * exp(-gamma * x.pr))),
                start=c(alpha=30.0, beta=0.9, gamma=0.01))
  params = c(coefficients(gomp.pr), summary(gomp.pr)[[3]])
  prior.params[i] <- params
}


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

