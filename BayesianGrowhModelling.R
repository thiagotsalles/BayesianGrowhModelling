# This script requires JAGS to be installed on the computer

library(rjags)
library(R2jags) # for jags
library(coda) # for as.mcmc
library(nortest) # for lillie.test
library(ggthemr) # for ggthemr
library(ggplot2) # for plots
library(cowplot) # for get_legend
ggthemr("fresh", layout="scientific") # Set theme for plots

# Sample
sampl <- read.csv("sample.csv", header=T, sep=";")
x <- sampl$Age
y <- sampl$Weight

# Prior information about weight
prior.info <- read.csv("prior_info.csv", header=T, sep=";")

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


#### ======================= PLOTS FOR DATA ======================= ####
# Function for formatting axes
formatAxes = function(ggpobj) {
  ggpobj +
    scale_y_continuous(
      breaks=seq(0, 40, 10), minor_breaks=seq(0, 40, 5), limits=c(0, 42)
    ) +
    scale_x_continuous(
      breaks=seq(0, 400, 100), minor_breaks=seq(0, 400, 50),
      limits=c(0, 420)
    ) +
    theme(plot.title = element_text(size=12),
          axis.text=element_text(size=12))
}

# Plot for sample
sample.plot <- ggplot(sampl, aes(x=Age, y=Weight)) +
  geom_point(size=2, color=swatch()[3]) + 
  labs(title="Sample", x="Age (days)", y="Weight (kg)")

# Plot for prior information
prior.plot <- ggplot(prior.info, aes(x=Age, y=Weight, color=Source)) +
  geom_point(aes(color="Sample"), alpha="0") + # dummy
  geom_point(size=2) +
  labs(title="Prior information", x="Age (days)", y=element_blank()) +
  scale_color_manual(breaks=c("Sample", levels(prior.info$Source)),
                     values=c(swatch()[2], swatch()[9], swatch()[4],
                              swatch()[5], swatch()[7], swatch()[3]))

# Legend for plots
legend.data <- get_legend(prior.plot +
                          theme(legend.title=element_blank(),
                                legend.direction="vertical",
                                legend.box.margin=margin(-20, 0, 0, 0)))

# Joining sample and prior info
samp.pri.plot <- plot_grid(
  formatAxes(sample.plot + theme(legend.position="none")),
  formatAxes(prior.plot + theme(legend.position="none")))

# Complete plot
data.plot <- plot_grid(samp.pri.plot, legend.data, ncol=2,
                       rel_widths=c(1, 0.2))


#### ==================== FREQUENTIST APPROACH ==================== ####
# Model fit (Gompertz, nonlinear least sqares)
gomp <- nls(y ~ (alpha * exp(-beta * exp(-gamma * x))),
            start=c(alpha=30.0, beta=0.9, gamma=0.01))

# Results table
results.freq <- data.frame("Estimate"=coefficients(gomp),
                           "CI_lower_limit"=confint(gomp)[, 1],
                           "CI_upper_limit"=confint(gomp)[, 2])

# Normality tests for residuals
shapiro.test(summary(gomp)$resid)
lillie.test(summary(gomp)$resid)

# Table with estimated weights and confidence intervals
est.y <- as.vector(fitted(gomp))
est.freq <- data.frame(
  "Age"=x, "Freq_est"=est.y,
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
prior.params <- data.frame(row.names=c("alpha", "beta",
                                       "gamma", "sigma"))

for (i in levels(prior.info$Source)) {
  x.pr = prior.info[prior.info$Source == i, ]$Age
  y.pr = prior.info[prior.info$Source == i, ]$Weight
  gomp.pr = nls(y.pr ~ (alpha * exp(-beta * exp(-gamma * x.pr))),
                start=c(alpha=30.0, beta=0.9, gamma=0.01))
  params = c(coefficients(gomp.pr), summary(gomp.pr)[[3]])
  prior.params[i] <- params
}


#### ========= BAYESIAN APPROACH - NONINFORMATIVE PRIORS ========== ####
# Creating a function to run the jags function
gomp.0.F <- function() {
  # Prior distributions
  alpha ~ dnorm(0.0, 0.0001)
  beta ~ dnorm(0.0, 0.0001)
  gamma ~ dunif(0.0, 1.0000)
  tau2 ~ dgamma(0.1, 0.1)
  
  # Model
  for(i in 1 : n) {
    y[i] ~ dnorm(mu[i], tau2)
    mu[i] <- alpha * exp(-beta * exp(-gamma * x[i]))
  }
  
  # Estimates for weight
  for(i in 1 : n) {
    yp[i] <- alpha * exp(-beta * exp(-gamma * x[i]))
  }
}

# Model fit
gomp.0 <- jags(data, inits.gomp,
               parameters=c(theta, "yp"),
               model.file=gomp.0.F,
               n.chains=1, n.iter=25000, n.burnin=4500, n.thin=5)

## Convergence tests ##
# Save bugs object as MCMC object
gomp.0.mcmc <- as.mcmc(gomp.0$BUGSoutput$sims.matrix[, ])
# Tests
geweke.diag(gomp.0.mcmc)
raftery.diag(gomp.0.mcmc)
heidel.diag(gomp.0.mcmc)

## Results table, with parameters and weight estimates ##
hpd.interval.0 <- HPDinterval(gomp.0.mcmc)
results.bayes.0 <- data.frame(gomp.0$BUGSoutput$summary[, 1],
                              hpd.interval.0)[-3, ]
colnames(results.bayes.0) <- c("Estimate", "HPD lower limit",
                               "HPD upper limit")

# RMSE and BIC
res.0 <- unname(y - results.bayes.0[, 1][5:12])
rmse.0 <- RMSE(res.0)
bic.0 <- BIC_func(y, res.0)

#### ========== BAYESIAN APPROACH - INFORMATIVE PRIORS ============ ####
# Creating a function to run the jags function
gomp.1.F <- function() {
  # Prior distributions
  alpha ~ dnorm(29.6760, 0.0360)
  beta ~ dnorm(1.9966, 26.2314)
  gamma ~ dunif(0.0081, 0.0193)
  tau2 ~ dgamma(19.9023, 21.7393)
  
  # Model
  for(i in 1 : n) {
    y[i] ~ dnorm(mu[i], tau2)
    mu[i] <- alpha * exp(-beta * exp(-gamma * x[i]))
  }
  
  # Estimates for weight
  for(i in 1 : n) {
    yp[i] <- alpha * exp(-beta * exp(-gamma * x[i]))
  }
}

# Model fit
gomp.1 <- jags(data, inits.gomp,
               parameters=c(theta, "yp"),
               model.file=gomp.1.F,
               n.chains=1, n.iter=25000, n.burnin=4500, n.thin=5)

## Convergence tests ##
# Save bugs object as MCMC object
gomp.1.mcmc <- as.mcmc(gomp.1$BUGSoutput$sims.matrix[, ])
# Tests
geweke.diag(gomp.1.mcmc)
raftery.diag(gomp.1.mcmc)
heidel.diag(gomp.1.mcmc)

## Results table ##
hpd.interval.1 <- HPDinterval(gomp.1.mcmc)
results.bayes.1 <- data.frame(gomp.1$BUGSoutput$summary[, 1],
                              hpd.interval.1)[-3, ]
colnames(results.bayes.1) <- c("Estimate", "HPD lower limit",
                               "HPD upper limit")

# RMSE and BIC
res.1 <- unname(y - results.bayes.1[, 1][5:12])
rmse.1 <- RMSE(res.1)
bic.1 <- BIC_func(y, res.1)


#### ======= BAYESIAN APPROACH - INF. PRIORS (6x VARIANCE) ======== ####
# Creating a function to run the jags function
gomp.2.F <- function() {
  # Prior distributions
  alpha ~ dnorm(29.6760, 0.0060)
  beta ~ dnorm(1.9966, 4.3719)
  gamma ~ dunif(0.0000, 0.0274)
  tau2 ~ dgamma(0.1, 0.1)

  # Model
  for(i in 1 : n) {
    y[i] ~ dnorm(mu[i], tau2)
    mu[i] <- alpha * exp(-beta * exp(-gamma * x[i]))
  }
  
  # Estimates for weight
  for(i in 1 : n) {
    yp[i] <- alpha * exp(-beta * exp(-gamma * x[i]))
  }
}

# Model fit
gomp.2 <- jags(data, inits.gomp,
               parameters=c(theta, "yp"),
               model.file=gomp.2.F,
               n.chains=1, n.iter=25000, n.burnin=4500, n.thin=5)

## Convergence tests ##
# Save bugs object as MCMC object
gomp.2.mcmc <- as.mcmc(gomp.2$BUGSoutput$sims.matrix[, ])
# Tests
geweke.diag(gomp.2.mcmc)
raftery.diag(gomp.2.mcmc)
heidel.diag(gomp.2.mcmc)

## Results table ##
hpd.interval.2 <- HPDinterval(gomp.2.mcmc)
results.bayes.2 <- data.frame(gomp.2$BUGSoutput$summary[, 1],
                              hpd.interval.2)[-3, ]
colnames(results.bayes.2) <- c("Estimate", "HPD lower limit",
                               "HPD upper limit")

# RMSE and BIC
res.2 <- unname(y - results.bayes.2[, 1][5:12])
rmse.2 <- RMSE(res.2)
bic.2 <- BIC_func(y, res.2)


#### ===================== PLOTS FOR RESULTS ====================== ####
# Function for making plots for the results
makePlot = function(data) {
  ggplot(data) +
    geom_line(aes(x=Age, y=Freq_est, color="Frequentist fit"),
              stat="smooth", method=loess, size=0.75, linetype=1) +
    geom_line(aes(x=x, y=CI_lower_limit, color="Confidence interval"),
              stat="smooth", method=loess, size=0.75, linetype=1) +
    geom_line(aes(x=x, y=CI_upper_limit, color="Confidence interval"),
              stat="smooth", method=loess, size=0.75, linetype=1) +
    geom_line(aes(x=Age, y=Estimate, color="Bayesian fit"),
              stat="smooth", method=loess, size=0.75) +
    geom_line(aes(x=x, y=HPD.lower.limit, color="95% HPD"),
              stat="smooth", method=loess, size=0.75, linetype=2) +
    geom_line(aes(x=x, y=HPD.upper.limit, color="95% HPD"),
              stat="smooth", method=loess, size=0.75, linetype=2) +
    geom_point(aes(x=Age, y=Weight, color="Sample"), size=2
    ) +
    labs(x="Age (days)", y="Weight (kg)") +
    theme(plot.title = element_text(size=12),
          axis.text=element_text(size=12)) +
    scale_y_continuous(
      breaks=seq(0, 40, 10), minor_breaks=seq(0, 40, 5),
      limits=c(0, 42)
    ) +
    scale_x_continuous(
      breaks=seq(0, 210, 105), minor_breaks=seq(0, 210, 26.25),
      limits=c(0, 220)
    ) +
    scale_color_manual(name=NULL,
      breaks=c("Sample", "Frequentist fit", "Confidence interval",
               "Bayesian fit", "95% HPD"),
      values=c(swatch()[4], swatch()[4], swatch()[2],
               swatch()[2], swatch()[3]))
}

# Joined data of sample and results
all.results.0 <- data.frame(sampl[-1], est.freq[-1],
                            results.bayes.0[5:12, ])

all.results.1 <- data.frame(sampl[-1], est.freq[-1],
                            results.bayes.1[5:12, ])

all.results.2 <- data.frame(sampl[-1], est.freq[-1],
                            results.bayes.2[5:12, ])

# Plots for results containing frequentist and the 3 Bayeian approaches 
results.0.plot <- makePlot(all.results.0)
results.1.plot <- makePlot(all.results.1)
results.2.plot <- makePlot(all.results.2) +
  guides(colour=guide_legend(  # Format legend
    override.aes=list(shape=c(19, rep(32, 4)),
                      linetype=c("blank", rep("solid", 3), "dashed"))))

# Legend for plots
legend.results <- get_legend(results.2.plot + theme(
  legend.title=element_blank(),
  legend.direction="horizontal",
  legend.box.margin=margin(0, 0, 0, 0)))

# Joining the three plots
res.012.plot <- plot_grid(nrow=1,
  results.0.plot +
    theme(legend.position="none") + labs(title="Noninformative"),
  results.1.plot +
    theme(legend.position="none") +
    labs(title="Informative", y=element_blank()),
  results.2.plot +
    theme(legend.position="none") +
    labs(title="Inf. with greater variance", y=element_blank())
  )

# Complete plot
results.plot <- plot_grid(res.012.plot, legend.results, nrow=2,
                          rel_heights=c(1, 0.2))


