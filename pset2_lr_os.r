##################################
# Topics in Labor Markets, PSet 2
# Lina Ramirez and Olivia Stiegman
# 5/1/24
##################################
# install.packages('readstata13')
# install.packages('cowplot')

# Clear environment
rm(list=ls())
# setwd("/Users/linaramirez/Library/CloudStorage/Dropbox/phd/spring-2024/labor3/pset2")
# setwd("/Users/oliviastiegman/Documents/UChicago/2nd Year/Spring Quarter/Topics in Labor Markets/pset2")

# Load packages
library(gtools)
library(data.table)
library(ggplot2)
library(Matrix)
library(reshape)
library(ggplot2)
library(readstata13)
library(cowplot)
library(SparseM)

# For testing EM algorithm
model.mixture.new <- function(nk) {
  model <- list()

  model$A <- array(10 * (1 + 1.5 * runif(3 * nk)), c(3, nk))
  model$S <- array(1, c(3, nk))
  model$pk <- rdirichlet(1, rep(1, nk))
  model$nk <- nk

  return(model)
}

# Generate simulated data
model.mixture.simulate <- function(model, N, sd.scale) {
  Y1 <- array(0, sum(N))
  Y2 <- array(0, sum(N))
  Y3 <- array(0, sum(N))
  K  <- array(0, sum(N))

  A <- model$A
  S <- model$S
  pk <- model$pk
  nk <- model$nk

  # draw K
  K <- sample.int(nk, N, TRUE, pk)

  # draw Y1, Y2, Y3
  Y1 <- A[1, K] + S[1, K] * rnorm(N) * sd.scale
  Y2 <- A[2, K] + S[2, K] * rnorm(N) * sd.scale
  Y3 <- A[3, K] + S[3, K] * rnorm(N) * sd.scale

  data.sim <- data.table(k = K, y1 = Y1, y2 = Y2, y3 = Y3)

  return(data.sim)
}

# Mixture model
model <- model.mixture.new(3)
df <- model.mixture.simulate(model, 10000, sd.scale = 0.5)
df_melt <- melt(df, id = "k")
plot <- ggplot(df_melt, aes(x = value, group = k, fill = factor(k))) +
  geom_density() +
  facet_grid(~variable) +
  theme_bw()
print(plot)

########################################################################################
#QUESTION 1 AND 2

# Log of a normal pdf with arbitrary mean and variance
lognormpdf <- function(Y, mu, sigma)  {
  -0.5 * ((Y - mu) / sigma)^2   - 0.5 * log(2.0 * pi) - log(sigma)
}

# Avoid underflow
logsumexp <- function(v) {
  vm <- max(v)
  log(sum(exp(v - vm))) + vm
}

# Expectation function
expectation <- function(N, nk, pk, Y1, Y2, Y3, A, S) {
  # For the posteriors
  omega <- array(0, c(N, nk))
  lpm <- array(0, c(N, nk))

  # Initialize the likelihood
  lik <- 0

  # Loop over each latent variable type
  for (i in 1:N) {
    # Log of the probability of each latent variable type
    lomega <- log(pk)

    # Log of the probability of each observation given a mean and variance
    lnorm1 <- lognormpdf(Y1[i], A[1, ], S[1, ])
    lnorm2 <- lognormpdf(Y2[i], A[2, ], S[2, ])
    lnorm3 <- lognormpdf(Y3[i], A[3, ], S[3, ])

    # Calculate log of the likelihood for that latent variable type
    # (also the numerator of the posterior)
    lall <- lomega + lnorm2 + lnorm1 + lnorm3
    lpm[i, ] <- lall

    # Add to the likelihood
    lik <- lik + logsumexp(lall)

    # Calculate the posterior
    omega[i, ] <- exp(lall - logsumexp(lall))
  }

  return(list(omega = omega, lpm = lpm, lik = lik))
}

# Maximization function
maximization <- function(N, nk, omega, Y1, Y2, Y3, A, S) {
  Y <- as.matrix(c(Y1, Y2, Y3))
  DY <- as.matrix(kronecker(Y ,rep(1,nk)))
  d <- as.matrix(kronecker(rep(1, N), diag(nk)))
  Dkj <- as.matrix(kronecker(diag(3), d))

  rw <- as.matrix(c(omega, omega, omega))
  fit <- lm(DY ~ Dkj - 1, weights = rw)
  A <- coef(fit)
  A <- matrix(A, nrow = 3, ncol = nk)

  resids <- DY - predict(fit)
  resids_scale <- resids^2 / rw

  fit_v <- lm(resids_scale ~ Dkj - 1, weights = rw)
  S <- sqrt(coef(fit_v))
  S <- matrix(S, nrow = 3, ncol = nk)

  return(list(A = A, S = S))
}

# EM algorithm
em <- function(N, nk, df, tol) {

  # Initialize variables
  round <- 1
  lik_old <- 0
  check <- Inf

  # Initialize parameters
  pk <- rep(1 / nk, nk)
  A <- matrix(5, 3, nk)
  S <- matrix(1, 3, nk)

  # Initialize variables for convergence and Q, H functions
  omega_old <- array(pk, c(N, nk))
  lpm_old <- array(0, c(N, nk))

  # Arrange data
  Y1 <- df$y1
  Y2 <- df$y2
  Y3 <- df$y3

  # Perform EM algorithm until convergence
  while (check > tol) {
    print(paste("Round:", round))

    # Expectation step
    results <- expectation(N, nk, pk, Y1, Y2, Y3, A, S)

    # Update probabilities of types
    pk <- colSums(results$omega) / N

    # Update variables for convergence and Q, H functions
    lik_new <- results$lik
    omega_new <- results$omega
    lpm_new <- results$lpm

    # Maximization step
    results <- maximization(N, nk, omega_new, Y1, Y2, Y3, A, S)

    # Update means and variances
    A <- results$A
    S <- results$S

    # Calculate Q and H functions
    Q1 <- sum(((omega_old) * lpm_old))
    Q2 <- sum(((omega_old) * lpm_new))
    H1 <- - sum((omega_old) * log(omega_old))
    H2 <- - sum((omega_old) * log(omega_new))

    # Calculate convergence criterion
    check <- abs(lik_new - lik_old)

    print(paste("Check:", check))

    # Update variables for next iteration
    round <- round + 1
    lik_old <- lik_new
    omega_old <- omega_new
    lpm_old <- lpm_new
  }

  return(list(pk = pk, A = A, S = S))

}

em_results <- em(10000, 3, df, tol=1e-6)

##### QUESTION 3: With simulated data ####

# Plot the simulated data

ggplot(df_melt,aes(x=value,group=k,fill=factor(k))) + geom_density() + facet_grid(~variable) + theme_bw()

# Save the plot as PDF
ggsave("./output/histograms.pdf")

### Estimated vs. True means 


# Extract the required columns from em_results$A and model$A
As <- cbind(unlist(em_results$A), unlist(model$A))

# Create a dataframe with Estimated and True values for each y variable
df_As <- data.frame(
  Estimated = c(As[, 1], As[, 2], As[, 3]),
  True = c(As[, 4], As[, 5], As[, 6]),
  y = rep(c('y1', 'y2', 'y3'), each = nrow(As) / 3)
)

# Create scatter plot using ggplot
ggplot(df_As, aes(x = Estimated, y = True, color=y)) +
  geom_point() +  # Add points for scatter plot
  labs(title = "Estimated vs True",  # Set title
       x = "Estimated",              # Set x-axis label
       y = "True") +                 # Set y-axis label
  geom_text(aes(label = "Mean"), x = Inf, y = Inf, vjust = -1) +  # Add label "Mean"
  theme_minimal()                   # Set plot theme to minimal

ggsave("./output/means.pdf")

### Estimated vs. True Variances


# Extract the required columns from em_results$S and model$S
Ss <- cbind(unlist(em_results$S), unlist(model$S))

# Create a dataframe with Estimated and True values for each y variable
df_Ss <- data.frame(
  Estimated = c(Ss[, 1], Ss[, 2], Ss[, 3]),
  True = c(Ss[, 4], Ss[, 5], Ss[, 6]),
  y = rep(c('y1', 'y2', 'y3'), each = nrow(As) / 3)
)

# Create scatter plot using ggplot
ggplot(df_Ss, aes(x = Estimated, y = True, color=y)) +
  geom_point() +  # Add points for scatter plot
  labs(title = "Estimated vs True",  # Set title
       x = "Estimated",              # Set x-axis label
       y = "True") +                 # Set y-axis label
  geom_text(aes(label = "Mean"), x = Inf, y = Inf, vjust = -1) +  # Add label "Mean"
  theme_minimal()                   # Set plot theme to minimal

ggsave("./output/variance.pdf")


### Estimated vs. True Pks


# Extract the required columns from em_results$S and model$S
em_results[[1]] <- matrix(em_results[[1]], nrow = 1)
Pks <- cbind(unlist(em_results$pk), unlist(model$pk))


# Create a dataframe with Estimated and True values for each y variable
df_Pks <- data.frame(
  Estimated = c(Pks[, 1], Pks[, 2], Pks[, 3]),
  True = c(Pks[, 4], Pks[, 5], Pks[, 6]),
  y = rep(c('y1', 'y2', 'y3'), each = nrow(As) / 3)
)

# Create scatter plot using ggplot
ggplot(df_Pks, aes(x = Estimated, y = True, color=y)) +
  geom_point() +  # Add points for scatter plot
  labs(title = "Estimated vs True",  # Set title
       x = "Estimated",              # Set x-axis label
       y = "True") +                 # Set y-axis label
  geom_text(aes(label = "Mean"), x = Inf, y = Inf, vjust = -1) +  # Add label "Mean"
  theme_minimal()                   # Set plot theme to minimal

ggsave("./output/pk.pdf")


#We need to report the values of the likelihood, the H function and the Q function but I don't know how :(


##### QUESTION 4: Estimating on PSID Data #### 
psid_df <- data.table(read.dta13("./AER_2012_1549_data/output/data4estimation.dta"))

#computing the wage residuals 
fit <- lm(log(log_y) ~ year + marit + state_st, psid_df, na.action = na.exclude)
psid_df[, log_yr := residuals(fit)]

# extract lags
setkey(psid_df,person,year)
psid_df[, log_yr_l1 := psid_df[J(person, year - 2), log_yr]]
psid_df[, log_yr_l2 := psid_df[J(person, year - 4), log_yr]]

# compute difference from start
fdata <- psid_df[!is.na(log_yr * log_yr_l1 * log_yr_l2)][, list(y1 = log_yr_l2, y2 = log_yr_l1, y3 = log_yr)]

K <- sample.int(3, 4941, TRUE, model$pk)
K_df <- data.frame(k = K)

fdata <- cbind(fdata, K_df)

#Running the EM algorithm
psid_em3_results <- em(4941, 3, fdata, tol = 1e-3)
psid_em4_results <- em(4941, 4, fdata, tol = 1e-3)
psid_em5_results <- em(4941, 5, fdata, tol = 1e-3)

