# Author: Georgia Papadogeorgou
# Date: 5/2/2016
# Description: Functions that will be used to update the range, variance and nugget
#              parameter for Gaussian Process with Exponential Covariance functions.
#              It also includes functions that updates the betas of the linear mean
#              function.

PredictFstar <- function(betas, Desstar, Y, Kstar, K2star, VarX_inv) {
  # Function that takes in the observed values, prediction points and covariance
  # matrices, and returns the predictions.
  #
  # Args:
  #  betas:    Vector of the betas in the linear predictor for the mean function.
  #  Desstar:  The design matrix of the locations where we want to predict. The number
  #            of columns for Desstar should be equal to the length of the betas.
  #  Y:        The observed outcome values at the observed locations.
  #  Kstar:    The covariance matrix based of the predicting locations (rows) against
  #            the observed locations (columns).
  #  K2star:   The covariance matrix of the predicting locations.
  #  VarX_inv: The inverse of the covariance matrix of the observed locations.
  #
  # Returns:
  #  Vector of fstar values for the prediciting locations.

  b <- matrix(betas, nrow = length(betas), ncol = 1)
  mu_star <- Desstar %*% b + Kstar %*% VarX_inv %*% (Y - Des %*% b)
  sigma_star <- K2star - Kstar %*% VarX_inv %*% t(Kstar)
  fstar <- rmvnorm(1, mean = mu_star, sigma = sigma_star)
  return(fstar)
}


UpdateVariance <- function(var_prev, prop_f, outcome, out_mean, out_Sigma,
                           range, nugget,
                           alphaf, betaf, alphan, betan, alphal, betal) {
  # Function that proposed new value for the marginal variance of an exponential
  # covariance function and accepts or rejects it, returning the new varinace
  # value and an indicator of whether we accepted.
  #
  # Args:
  #  var_prev:   Previous value of the variance around which the truncated normal
  #              will be centered.
  #  prop_f:     Variance of the truncated normal from which we will propose.
  #  outcome:    Vector of outcome values at osberved and predicting locations.
  #  out_mean:   The mean of the mutlivariate normal from which the data arise.
  #  out_Sigma:  The outcome covariance matrix.
  #  range:      The current value of the range parameter.
  #  nugget:     The current value of the nugget parameter.
  #  alpha,beta: The values of alpha and beta for the inverse gamma prior on the
  #              marginal variance (f), the nugget (n), and range (l).
  #
  # Returns a list of two elements. The element value corresponds to the new
  # variance value, and the element accept is an indicator of acceptance of proposed
  # move.
  
  var_prop <- rtruncnorm(1, a = 0, mean = var_prev, sd = sqrt(prop_f))
  
  logAR <- CalcLogPost(x = outcome, mean = out_mean, Sigma = out_Sigma,
                       range = range, var = var_prop, nugget = nugget,
                       alphaf = alphaf, betaf = betaf, alphan = alphan,
                       betan = betan, alphal = alphal, betal = betal) -
    CalcLogPost(x = outcome, mean = out_mean, Sigma = out_Sigma,
                range = range, var = var_prev, nugget = nugget,
                alphaf = alphaf, betaf = betaf, alphan = alphan,
                betan = betan, alphal = alphal, betal = betal) +
    dtnorm(var_prev, lower = 0, mean = var_prop, sd = sqrt(prop_f),
           log = TRUE) -
    dtnorm(var_prop, lower = 0, mean = var_prev, sd = sqrt(prop_f),
           log = TRUE)
  
  r <- NULL
  r$value <- var_prev
  r$accept <- FALSE
  if (runif(1) < exp(logAR)) {
    r$value <- var_prop
    r$accept <- TRUE
  }
  return(r)
}






UpdateRange <- function(range_prev, prop_l, outcome, out_mean, out_Sigma, variance,
                        nugget, alphaf, betaf, alphan, betan, alphal, betal) {
  # Function that proposed new value for the range parameter of an exponential
  # covariance function and accepts or rejects it, returning the new ramge value
  # and an indicator of whether we accepted.
  #
  # Args:
  #  range_prev: Previous value of around around which the truncated normal proposal
  #              will be centered.
  #  prop_l:     Variance of the truncated normal from which we will propose.
  #  outcome:    Vector of outcome values at osberved and predicting locations.
  #  out_mean:   The mean of the mutlivariate normal from which the data arise.
  #  out_Sigma:  The outcome covariance matrix.
  #  range:      The current value of the range parameter.
  #  nugget:     The current value of the nugget parameter.
  #  alpha,beta: The values of alpha and beta for the inverse gamma prior on the
  #              marginal variance (f), the nugget (n), and range (l).
  #
  # Returns a list of two elements. The element value corresponds to the new range
  # value, and the element accept is an indicator of acceptance of proposed move.

  # Proposing a new value from a truncated normal centered at the current value.
  range_prop <- rtruncnorm(1, a = 0, mean = range_prev, sd = sqrt(prop_l))
  
  logAR <- CalcLogPost(x = outcome, mean = out_mean,
                       Sigma = out_Sigma, range = range_prop,
                       var = variance, nugget = nugget,
                       alphaf = alphaf, betaf = betaf,
                       alphan = alphan, betan = betan,
                       alphal = alphal, betal = betal) -
    CalcLogPost(x = outcome, mean = out_mean,
                Sigma = out_Sigma, range = range_prev,
                var = variance, nugget = nugget,
                alphaf = alphaf, betaf = betaf, alphan = alphan,
                betan = betan, alphal = alphal, betal = betal) +
    dtnorm(range_prev, lower = 0, mean = range_prop, sd = sqrt(prop_l),
           log = TRUE) -
    dtnorm(range_prop, lower = 0, mean = range_prev, sd = sqrt(prop_l),
           log = TRUE)
  
  r <- NULL
  r$value <- range_prev
  r$accept <- FALSE
  if (runif(1) < exp(logAR)) {
    r$value <- range_prop
    r$accept <- TRUE
  }
  return(r)
}



UpdateNugget <- function(nug_prev, prop_n, outcome, out_mean, out_Sigma, variance,
                         range, alphaf, betaf, alphan, betan, alphal, betal) {
  # Function that proposed new value for the nugget parameter of an exponential
  # covariance function and accepts or rejects it, returning the new ramge value
  # and an indicator of whether we accepted.
  #
  # Args:
  #  nug_prev:   Previous value of around around which the truncated normal proposal
  #              will be centered.
  #  prop_n:     Variance of the truncated normal from which we will propose.
  #  outcome:    Vector of outcome values at osberved and predicting locations.
  #  out_mean:   The mean of the mutlivariate normal from which the data arise.
  #  out_Sigma:  The outcome covariance matrix.
  #  range:      The current value of the range parameter.
  #  nugget:     The current value of the nugget parameter.
  #  alpha,beta: The values of alpha and beta for the inverse gamma prior on the
  #              marginal variance (f), the nugget (n), and range (l).
  #
  # Returns a list of two elements. The element value corresponds to the new nugget
  # value, and the element accept is an indicator of acceptance of proposed move.
  
  # Proposing a new value from a truncated normal centered at the current value.
  nug_prop <- rtruncnorm(1, a = 0, mean = nug_prev, sd = sqrt(prop_n))
  
  logAR <- CalcLogPost(x = outcome, mean = out_mean, Sigma = out_Sigma,
                       range = range, var = variance, nugget = nug_prop,
                       alphaf = alphaf, betaf = betaf, alphan = alphan,
                       betan = betan, alphal = alphal, betal = betal) -
    CalcLogPost(x = outcome, mean = out_mean, Sigma = out_Sigma,
                range = range, var = variance, nugget = nug_prev,
                alphaf = alphaf, betaf = betaf, alphan = alphan,
                betan = betan, alphal = alphal, betal = betal) +
    dtnorm(nug_prev, lower = 0, mean = nug_prop, sd = sqrt(prop_n), log = TRUE) -
    dtnorm(nug_prop, lower = 0, mean = nug_prev, sd = sqrt(prop_n), log = TRUE)
  
  r <- NULL
  r$value <- nug_prev
  r$accept <- FALSE
  if (runif(1) < exp(logAR)) {
    r$value <- nug_prop
    r$accept <- TRUE
  }
  return(r)
}



UpdateBetas <- function(full_Des, Sigma, outcome, uninf = 100 ^ 2) {
  # Function that updates the betas of the linear mean functions.
  #
  # Args:
  #  full_Des: Design matrix for observed and predicting locations.
  #  Sigma:    Covariance matrix over all locations.
  #  outcome:  Observed or generated outcome value at all locations.
  #  uninf:    Variance of the normal mean zero prior on betas.
  #
  # Returns vector of betas sampled from the posterior distribution.
  
  Sigma_inv <- solve(Sigma)
  betas_var <- solve(t(full_Des) %*% Sigma_inv %*% full_Des +
                       uninf ^ ( - 1) * diag(ncol(full_Des)))
  betas_mean <- betas_var %*% t(full_Des) %*% Sigma_inv %*% outcome
  betas <- rmvnorm(1, mean = betas_mean, sigma = betas_var)
  return(betas)
} 

