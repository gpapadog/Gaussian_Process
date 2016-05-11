# Author:      Georgia Papadogeorgou
# Date:        4/20/2016
# Description: Fitting a Gaussian Process with linear mean function on multiple covariates
#              and exponential correlation function on some of them, unknown hyperparameters.


# SET THE WORKING DIRECTORY
setwd('/Users/georgiapapadogeorgou/Documents/6.882/Project/')

source('Code/Exp_GP_functions.R')
source('Code/Update_Exp_GP_functions.R')
library(plot3D)


# Specify:
sample_size <- 30
num_conf <- 1
conf <- matrix(rep(1, num_conf * 2), num_conf, 2)
prop_f <- 30
prop_l <- 30
prop_n <- 80
alphaf <- betaf <- 1
alphan <- 1
betan <- 1
alphal <- 1
betal <- 1

# Creating linear confounders:
covs <- matrix(rnorm(sample_size * num_conf), nrow = sample_size, ncol = num_conf)
colnames(covs) <- paste0('X', 1:num_conf)
probs <- covs %*% conf[, 1]
probs <- exp(probs) / (1 + exp(probs))
Z <- rbinom(sample_size, 1, prob = probs)
Y <- Z + covs %*% conf[, 2] + rnorm(sample_size)

# Standardizing the covariates and creating a design matrix.
covs <- apply(covs, 2, function(x) (x - mean(x)) / sd(x))
Des <- cbind(Int = 1, Z, covs)

# The locations that the mean Gaussian Process should be predicted the same
# as the measured locations:
Xstar <- covs
Desstar <- Des
nstar <- nrow(Xstar)
full_Des <- rbind(Des, Desstar)


Nsamples <- 10000
range <- numeric(Nsamples) + 1
var <- numeric(Nsamples) + 1
nugget <- numeric(Nsamples) + 1
fstar <- array(NA, dim = c(Nsamples, nstar))
accept <- numeric(3)
names(accept) <- c('Variance', 'Range', 'Nugget')
beta <- matrix(NA, nrow = Nsamples, ncol = ncol(Des))
beta[1, ] <- 1
colnames(beta) <- c('intercept', 'slope', paste0('X', 1:num_conf))
uninf <- 100 ^ 2

# At every step we will update fstar given the rest of the parameters


for (ss in 2:Nsamples) {
  if (ss %% 1000 == 0) {
    print(ss)
    print(accept / ss)
    plot(nugget[1:(ss - 1)], type = 'l')
    plot(range[1:(ss - 1)], type = 'l')
    plot(var[1:(ss - 1)], type = 'l')
    
  }
  Sigma <- ExpCovarMat(x = covs, y = Xstar, range = range[ss - 1], var = var[ss - 1],
                       nugget = nugget[ss - 1])
  VarX <- Sigma[1:sample_size, 1:sample_size]
  Kstar <- Sigma[(sample_size + 1):dim(Sigma)[1], 1:sample_size]
  K2star <- Sigma[(sample_size + 1):dim(Sigma)[1], (sample_size + 1):dim(Sigma)[1]]
  VarX_inv <- solve(VarX)
  
  # Generating fstar
  fstar[ss, ] <- PredictFstar(betas = beta[ss - 1,], Desstar = Desstar, Y = Y, 
                              Kstar = Kstar, K2star = K2star, VarX_inv = VarX_inv)
  
  # Mean of Y and fstar:
  mvr_mean <- full_Des %*% matrix(beta[ss - 1, ], nrow = ncol(beta), ncol = 1)
  
  # Updating the marginal variance:
  upd_var <- UpdateVariance(var_prev = var[ss - 1], prop_f, outcome = c(Y, fstar[ss, ]),
                            out_mean = mvr_mean, out_Sigma = Sigma, range[ss - 1],
                            nugget[ss - 1], alphaf, betaf, alphan, betan, alphal, betal)
  var[ss] <- upd_var$value
  if (upd_var$accept) {
    accept[1] <- accept[1] + 1
  }
  # Updating Sigma
  Sigma <- ExpCovarMat(x = covs, y = Xstar, range = range[ss - 1], var = var[ss],
                       nugget = nugget[ss - 1])
  
  
  
  # Updating the range parameter.
  upd_range <- UpdateRange(range[ss - 1], prop_l, c(Y, fstar[ss, ]), out_mean = mvr_mean,
                           out_Sigma = Sigma, var[ss], nugget[ss - 1],
                           alphaf, betaf, alphan, betan, alphal, betal)
  range[ss] <- upd_range$value
  if (upd_range$accept) {
    accept[2] <- accept[2] + 1
  }
  Sigma <- ExpCovarMat(x = covs, y = Xstar, range = range[ss], var = var[ss],
                       nugget = nugget[ss - 1])
  
  
  # Updating the nugget:
  upd_nug <- UpdateNugget(nugget[ss - 1], prop_n, c(Y, fstar[ss, ]), mvr_mean, Sigma,
                          var[ss], range[ss], alphaf, betaf, alphan, betan, alphal, betal)
  nugget[ss] <- upd_nug$value
  if (upd_nug$accept) {
    accept[3] <- accept[3] + 1
  }
  Sigma <- ExpCovarMat(x = covs, y = Xstar, range = range[ss], var = var[ss],
                       nugget = nugget[ss])
  
  
  # Updating the betas:
  beta[ss, ] <- UpdateBetas(full_Des = full_Des, Sigma = Sigma,
                            outcome = c(Y, fstar[ss, ]), uninf = uninf)
    
  
}


plot(nugget, type = 'l')
plot(range, type = 'l')
plot(var, type = 'l')

scatter3D(covs[, 1], covs[, 2], Y, col = 'black', pch = 16)
points3D(covs[, 1], covs[, 2], apply(fstar, 2, mean, na.rm = TRUE), add = TRUE, col = 'red',
         pch = 16, cex = 0.7)

accept/Nsamples
