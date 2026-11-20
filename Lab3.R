################################################################################
##          ECONOMETRICS 4b: EMPIRICAL MACROECONOMICS - HT2025                ##
################################################################################

################################################################################
##                          Code for Assignment 2                             ##
################################################################################

####################### Setting the work environment ###########################

# Clear the environment
rm(list = ls())

# Loading the necessary packages
packages <- c("readxl", "vars", "lpirfs", "dplyr", "ggplot2", "zoo", "hpfilter")
# Make sure to install any packages you haven't installed yet using 
# install.packages("name_of_the_package") in the console, then run:
lapply(packages, library, character.only=TRUE)

##################### Loading and preparing the data ###########################
raw_data <- read.csv(file = "Data3.csv", sep = ";") # ; is used to separate columns

# Take a look at the data:
View(raw_data)


# We will be using the function lp_lin() from the package lpirfs. This function
# computes linear local projections. The data needs to be a dataframe so no need 
# to transform to TS.

col <- 2:6
varnames <- c("unemp", "indprod", "MP", "realrate", "rstar")
data <- lapply(col, function(i) as.numeric(t(raw_data[,i])))
names(data) <- varnames
list2env(data, envir = .GlobalEnv) 

# Transform unemp to 1st diff and indprod to logdiff
unemp <- diff(unemp)
indprod <- diff(log(indprod))

# Construct dataframe
MP <- MP[2:466]
#because the other variables loose one obs when differencing

data_lp <- as.data.frame(cbind(unemp, indprod, MP))

###################### Compute local projections ###############################

results_lp <- lp_lin(data_lp,
                     lags_endog_lin = 4,
                     trend = 0,
                     shock_type = 1, # unit shock type
                     confint = 1.96,
                     hor = 8)

# Point estimates responses to MP shock on each variable at each horizon:
beta_lp <- results_lp$irf_lin_mean[, , 3]

# Plot impulse response functions: response of unemployment and industrial prod
print(plot_lin(results_lp)[[3]]) 
print(plot_lin(results_lp)[[6]]) 

# Understanding the magnitude of the effect: 

# indprod is in log diff: a value of 0.0025 corresponds to a 0.25 percentage points
# (or 25 basis points) change in industrial production following 
# a 1 unit shock (i.e. a 1 percentage point increase in MP).
# unemp is period to period change in the unemployment rate. A value of 0.05 
# corresponds to a 0.05 percentage point (or 5 bp) change in unemployment following
# a 1 unit shock.

############################## Compute VAR #####################################

# Construct the deviation of the real interest rate from its time-varying 
# equilibrium

rdev <- realrate-rstar 
rdev <- rdev[2:466]
rdiff = diff(realrate)
# Construct the dataset (in time series)
data_var <- ts(cbind(unemp, indprod, rdiff), start = c(1990, 1), frequency = 12)
#plot.ts(data_var)

# Here, I select the lag length automatically by AIC.
p <- as.numeric(VARselect(data_var, lag.max = 10)[[1]][1])

# Compute a VAR model:
VARmod <- VAR(data_var, p = 4, season = 12)

# Finally, compute the SVAR identified using Cholesky decomposition, short-run 
# restrictions

# Prepare the A and B matrices, imposing Kx(k-1)/2=3 restrictions
# The A matrix defines the contemporanous relationships among the endogenous 
# variables: lower triangular structure so that unemp is contemporaneously 
# affected by indprod and MP, indprod contemporaneously affected by MP but not 
# vice versa. 
Amat <- diag(3) # 3x3 identity matrix
diag(Amat) <- NA
Amat[2,1] <- NA
Amat[3,1] <- NA
Amat[3,2] <- NA
# B matrix defines the relationships between the structural shocks, such that
# they are orthogonal.
Bmat <- diag(3)
diag(Bmat) <- NA

SVARmod <- SVAR(VARmod, Amat = Amat, Bmat = Bmat, max.iter = 10000, hessian = TRUE)
# Here you should get a warning but you can ignore. The model is just identified 
# because we coded it that way.

# Plot the IRF: responses of unemployment and industrial production to monetary
# policy shocks
 

irf_indprod <- irf(SVARmod, response = "indprod", impulse = "rdiff", 
                 n.ahead = 8, ortho = TRUE,runs=1000)
plot(irf_indprod)

irf_unemp <- irf(SVARmod, response = "unemp", impulse = "rdiff", 
                   n.ahead = 8, ortho = TRUE,runs=1000)
plot(irf_unemp)