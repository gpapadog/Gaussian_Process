# Gaussian_Process

Code to support my 6.882 project on the effect of using Gaussian Processes for confounding adjustment in linear model settings. The files included are:

- Exp_GP_Causal.R: Main simulation code. Generating data with a binary treatment and continuous covariates. The Gaussian Process part is fit on the covariates. The mean of the outcome is linear in treatment and covariates, to encourage extrapolation of the predictions towards a linear model near the edges of the observations.
- Exp_GP_functions.R: Functions to calculate the exponential kernel and the log posterior distribution.
- Update_Exp_GP_functions.R: Functions that update all the parameters (predict the Gaussian Process, update the coefficients of the linear mean, update the hyperparameters of the covariance matrix such as the signal variance, the nugget variance and the range).
