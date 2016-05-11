# Author:      Georgia Papadogeorgou
# Date:        4/20/2016
# Description: Fucntions for fitting a mean zero Gaussian Process with exponential
#              correlation function and unknown hyperparameters.

library(mvtnorm)
library(MASS)
library(truncnorm)
library(msm)
library(actuar)


ExponCovar <- function(x, y = NULL, range, var) {
  # Function that calculates the exponential covariance function for one or two
  # sets of observation.
  #
  # Args:
  #  x:     Set of covariates for the first group.
  #  y:     The second set of observations. Defaults to NULL. When set to NULL,
  #         the covariance of x on itself will be calculated.
  #  range: The range of the exponential covariance function.
  #  var:   The marginal variance of the process.
  
  if (class(x) == 'numeric') {
    x <- matrix(x, nrow = length(x), ncol = 1)
  }
  
  if (is.null(y)) {
    y <- x
  }
  
  if (class(y) == 'numeric') {
    y <- matrix(y, nrow = length(y), ncol = 1)
  }
  
  if (ncol(x) != ncol(y)) {
    stop('The two location sets should be of the same column size.')
  }
  
  D <- 0
  for (ii in 1:ncol(x)) {
    covar <- matrix(x[, ii], nrow = nrow(x), ncol = 1) %*% rep(1, nrow(y)) -
      matrix(rep(1, nrow(x)), nrow = nrow(x), ncol = 1) %*% y[, ii]
    covar <- abs(covar)
    D <- D + covar
  }
  
  covar <- exp( - D / range) * var
  return(covar)
}


ExpCovarMat <- function(x, y, range, var, nugget) {
  # Function that calculates the covariance matrix for two sets of observations,
  # where the first set of observations can potentially include a nugget effect
  # for the variance.
  #
  # Args:
  #  x:      Set of covariates for the first group.
  #  y:      The second set of observations.
  #  range:  The range of the exponential covariance function.
  #  var:    The marginal variance of the process.
  #  nugget: The nugget for the fist set of observations.

  K <- ExponCovar(x = x, y = NULL, range = range, var = var)
  if (identical(x, y)) {
    Kstar <- K
    K2star <- K
  } else {
    Kstar <- ExponCovar(x = y, y = x, range = range, var = var)
    K2star <- ExponCovar(x = y, y = NULL, range = range, var = var)
  }
  
  VarX <- K + nugget * diag(dim(K)[1])
  
  Sigma <- cbind(VarX, t(Kstar))
  Sigma <- rbind(Sigma, cbind(Kstar, K2star))
  return(Sigma)
}


CalcLogPost <- function(x, mean = 0, Sigma, range, var, nugget, alphaf, betaf,
                        alphan, betan, alphal, betal) {
  # Calculating the log posterior for a mean 0 gaussian process with exponential
  # covariance function and inverse gamma priors on the variance, the nugget and
  # the range parameter.
  #
  # Args:
  #  x:      Vector of the values of the process at all points
  #  mean:   Mean of the multivariate normal for Y and fstar. Defaults to 0.
  #  Sigma:  Covariance matrix of all observed and imputed values
  #  range:  Range of the exponential covariance function
  #  var:    Variance of the process
  #  nugget: Nugget variance at observed points
  #  alphaf: The alpha parameter of the IG prior on the marginal variance
  #  betaf:  The beta parameter of the IG prior on the marginal variance
  #  alphan: The alpha parameter of the IG prior on the nugget
  #  betan:  The beta paramater of the IG prior on the nugget
  #  alphal: The alpha parameter of the IG prior on the range
  #  betal:  The beta parameter of the IG prior on the range
  #
  # Returns the posterior evaluated at the specific values of data and
  # parameters.
  
  post <- dmvnorm(x, sigma = Sigma, log = TRUE)
  post <- post + actuar::dinvgamma(var, shape = alphaf, scale = betaf, log = TRUE)
  post <- post + actuar::dinvgamma(nugget, shape = alphan, scale = betan, log = TRUE)
  post <- post + actuar::dinvgamma(range, shape = alphal, scale = betal, log = TRUE)
  
  return(post)
}

