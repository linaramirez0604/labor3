##################################
# Topics in Labor Markets, PSet 2
# Lina Ramirez and Olivia Stiegman
# 5/1/24
##################################

# Clear environment
rm(list=ls())

# Load packages
library(gtools)
library(data.table)
library(ggplot2)
library(Matrix)
library(reshape)
library(readstata13)
library(SparseM)

# For testing EM algorithm
model.mixture.new <-function(nk) {
  
  model = list()
  # model for Y1,Y2,Y3|k 
  model$A     = array(3*(1 + 0.8*runif(3*nk)),c(3,nk))
  model$S     = array(1,c(3,nk))
  model$pk    = rdirichlet(1,rep(1,nk))
  model$nk    = nk
  return(model)
}

# Generate simulated data
model.mixture.simulate <-function(model,N,sd.scale) {
  
  Y1 = array(0,sum(N)) 
  Y2 = array(0,sum(N)) 
  Y3 = array(0,sum(N)) 
  K  = array(0,sum(N)) 
  
  A   = model$A
  S   = model$S
  pk  = model$pk
  nk  = model$nk

  # draw K
  K = sample.int(nk,N,TRUE,pk)

  # draw Y1, Y2, Y3
  Y1  = A[1,K] + S[1,K] * rnorm(N) *sd.scale
  Y2  = A[2,K] + S[2,K] * rnorm(N) *sd.scale
  Y3  = A[3,K] + S[3,K] * rnorm(N) *sd.scale

  data.sim = data.table(k=K,y1=Y1,y2=Y2,y3=Y3)
  return(data.sim)  
}

# Mixture model
model <- model.mixture.new(3)
df <- model.mixture.simulate(model,10000,sd.scale=0.5) # simulating with lower sd to see separation
df_melt <- melt(df,id="k")
plot <- ggplot(df_melt,aes(x=value,group=k,fill=factor(k))) + geom_density() + facet_grid(~variable) + theme_bw()
print(plot)

########################################################################################

# Log of a normal pdf with arbitrary mean and variance
lognormpdf <- function(Y,mu,sigma)  {
  -0.5 * (  (Y-mu) / sigma )^2   - 0.5 * log(2.0*pi) - log(sigma)  
}

# Avoid underflow
logsumexp <- function(v) {
  vm = max(v)
  log(sum(exp(v-vm))) + vm
}

# Expectation function
expectation <- function(N, nk, pk, Y1, Y2, Y3, A, S) {
    # For the posteriors
    omega <- array(0,c(N,nk))
    lpm <- array(0,c(N,nk))

    # Initialize the likelihood
    lik <- 0

    # Loop over each latent variable type
    for (i in 1:N) {
        # Log of the probability of each latent variable type
        lomega <- log(pk)

        # Log of the probability of each observation given a mean and variance
        lnorm1 <- lognormpdf(Y1[i], A[1,], S[1,])
        lnorm2 <- lognormpdf(Y2[i], A[2,], S[2,])
        lnorm3 <- lognormpdf(Y3[i], A[3,], S[3,])
        
        # Calculate log of the likelihood for that latent variable type (also the numerator of the posterior)
        lall <- lomega + lnorm2 + lnorm1 +lnorm3
        lpm[i,] <- lall

        # Add to the likelihood
        lik <- lik + logsumexp(lall)
    }   

    # Loop over each latent variable type
    for (i in 1:N) {
        # Calculate the posterior
        omega[i,] <- exp(lpm[i,] - logsumexp(lpm))
    }

    return(list(omega=omega, lpm=lpm, lik=lik))
}

# Maximization function
maximization <- function(N, nk, omega, Y1, Y2, Y3, A, S) {
    DY1 <- as.matrix(kronecker(Y1 ,rep(1,nk)))
    DY2 <- as.matrix(kronecker(Y2 ,rep(1,nk)))
    DY3 <- as.matrix(kronecker(Y3 ,rep(1,nk)))
    Dkj1 <- as.matrix.csr(kronecker(rep(1,N),diag(nk)))
    Dkj2 <- as.matrix.csr(kronecker(rep(1,N),diag(nk)))  
    Dkj3 <- as.matrix.csr(kronecker(rep(1,N),diag(nk))) 

    rw <- c(t(omega))
    fit <- slm.wfit(Dkj1, cbind(DY1, DY2, DY3), rw)
    A <- coef(fit)
    fit_v <- slm.wfit(Dkj1, resid(fit)^2/rw, rw)
    S <- sqrt(coef(fit_v))

    return(list(A=A, S=S))
}

# EM algorithm
em <- function(N, nk, df, tol=1e-6) {

    # Initialize variables
    round <- 1
    lik_old <- 0
    check <- Inf

    # Initialize parameters
    pk <- rep(1/nk, nk)
    A <- matrix(0, nk, nk)
    S <- matrix(1, nk, nk)

    # Initialize variables for convergence and Q, H functions
    omega_old <- pk
    lpm_old <- array(0,c(N,nk))

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
        
        print(paste("Q1:", Q1))
        print(paste("Q2:", Q2))
        print(paste("H1:", H1))
        print(paste("H2:", H2))
        
        # Calculate convergence criterion
        check <- abs(lik_new - lik_old)

        print(paste("Check:", check))
        print(paste("pk:", pk))
        print(paste("A:", A))
        print(paste("S:", S))

        # Update variables for next iteration
        round <- round + 1
        lik_old <- lik_new
        omega_old <- omega_new
        lpm_old <- lpm_new
    }

    return(list(pk=pk, A=A, S=S))

}

em_results <- em(10000, 3, df, tol=1e-6)

