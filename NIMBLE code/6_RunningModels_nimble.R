library(nimble)

# Set seed
mySeed <- 0
set.seed(mySeed)

##### Running the models; 
# MCMC settings - run
# niter <- 2500
# nthin <- 2
# nburn <- 1000
# nchains <- 3    

# MCMC settings - test run
niter <- 2
nthin <- 1
nburn <- 0
nchains <- 3 

# Load data
rype.data <- readRDS("RypeData_forIM.rds")

# Make prior for survival based on telemetry study
shape_from_stats <- function(mu , sigma ){
  a <-(mu^2-mu^3-mu*sigma^2)/sigma^2
  b <- (mu-2*mu^2+mu^3-sigma^2+mu*sigma^2)/sigma^2
  shape_ps <- c(a,b)
  return(shape_ps)
}

S_priors <- shape_from_stats(0.425, 0.035)

## Plot - to see if this works; 
x <- seq(0, 1, 0.001)
w <- dbeta(x, S_priors[1], S_priors[2])
plot(x, w, type = "l")


# Assembling data and constants for all models
# NOTE: Data/constants not used in any particular model are ignored by nimble

nim.data <- list(R_obs = rype.data$R_obs, y = rype.data$y, 
                 zeros.dist = rype.data$zeros_dist, L = rype.data$L, 
                 N_line_year = rype.data$N_line_year, A = rype.data$A)

nim.constants <- list(N_years = rype.data$N_years, W = rype.data$W, scale1 = 1000,
                     N_obs = rype.data$N_obs, Year_obs = rype.data$Year_obs,
                     N_sites = rype.data$N_sites, 
                     R_obs_year = rype.data$R_obs_year, N_R_obs = rype.data$N_R_obs,
                     a = S_priors[1], b = S_priors[2])


################################################################################
################################################################################
### No prior on S - variable S between years (S[t])

# Setting parameters to monitor            
params <- c("esw", "R_year", "p", "S", "D")

# Function for setting initial values
inits2 <- function(){
  
  mu.D1 <- runif(1, 3, 4)
  N_exp <- matrix(NA, nrow = N_sites, ncol = N_years)
  for(j in 1:N_sites){
    N_exp[j, 1] <- rpois(1, mu.D1*L[j, 1]*(W/nim.constants$scale1)*2)      ## Expected number of birds
  }
  
  list(
  mu.dd = runif(1, 4, 5), 
  sigma.dd = runif(1, 0.05, 2),
  mu.D1 = mu.D1, 
  sigma.D = runif(1, 0.05, 2),
  mu.R = runif(1, -2, 2), 
  sigma.R = runif(1, 0.05, 2),
  eps.dd = rep(0, N_years), 
  eps.R = rep(0, N_years), 
  eps.D1 = rep(0, N_sites),
  S = runif(N_years, 0.4, 0.5),
  N_exp = N_exp
  )
}

# Sample initial values
initVals2 <- list(inits2(), inits2(), inits2())

# Run code file
source('3_Combined_M2_nimble.R')

# Test run
out_real2 <- nimbleMCMC(code = modM2.code, 
                        data = nim.data, constants = nim.constants,
                        inits = initVals2, monitors = params,
                        nchains = nchains, niter = niter, 
                        nburnin = nburn, thin = nthin, 
                        samplesAsCodaMCMC = TRUE, setSeed = mySeed)

################################################################################
################################################################################
## No prior on S - common S across years (S)

# Setting parameters to monitor            
params <- c("esw", "R_year", "p", "S", "D")

# Function for setting initial values
inits3 <- function(){
  
  mu.D1 <- runif(1, 3, 4)
  N_exp <- matrix(NA, nrow = N_sites, ncol = N_years)
  for(j in 1:N_sites){
    N_exp[j, 1] <- rpois(1, mu.D1*L[j, 1]*(W/nim.constants$scale1)*2)      ## Expected number of birds
  }
  
  list(
    mu.dd = runif(1, 4, 5), 
    sigma.dd = runif(1, 0.05, 2),
    mu.D1 = mu.D1, 
    sigma.D = runif(1, 0.05, 2),
    mu.R = runif(1, -2, 2), 
    sigma.R = runif(1, 0.05, 2),
    eps.dd = rep(0, N_years), 
    eps.R = rep(0, N_years), 
    eps.D1 = rep(0, N_sites),
    S = runif(1, 0.4, 0.5),
    N_exp = N_exp
  )
}

# Sample initial values
initVals3 <- list(inits3(), inits3(), inits3())

# Run code file
source('3_Combined_M3_nimble.R')

# Test run
out_real3 <- nimbleMCMC(code = modM3.code, 
                        data = nim.data, constants = nim.constants,
                        inits = initVals3, monitors = params,
                        nchains = nchains, niter = niter, 
                        nburnin = nburn, thin = nthin, 
                        samplesAsCodaMCMC = TRUE, setSeed = mySeed)


################################################################################
################################################################################
## ## Prior on S based on radiotelemetry study - extracted from Israelsen et al. 2020 [Ecol & Evo]

# Function for setting initial values
inits4 <- function(){
  
  mu.D1 <- runif(1, 3, 4)
  N_exp <- matrix(NA, nrow = N_sites, ncol = N_years)
  for(j in 1:N_sites){
    N_exp[j, 1] <- rpois(1, mu.D1*L[j, 1]*(W/nim.constants$scale1)*2)      ## Expected number of birds
  }
  
  list(
    mu.dd = runif(1, 4, 5), 
    sigma.dd = runif(1, 0.05, 2),
    mu.D1 = mu.D1, 
    sigma.D = runif(1, 0.05, 2),
    mu.R = runif(1, -2, 2), 
    sigma.R = runif(1, 0.05, 2),
    eps.dd = rep(0, N_years), 
    eps.R = rep(0, N_years), 
    eps.D1 = rep(0, N_sites),
    S = rbeta(1, S_priors[1], S_priors[2]),
    N_exp = N_exp
  )
}

# Sample initial values
initVals4 <- list(inits4(), inits4(), inits4())

# Run code file
source('3_Combined_M4_nimble.R')

# Test run
out_real4 <- nimbleMCMC(code = modM4.code, 
                        data = nim.data, constants = nim.constants,
                        inits = initVals4, monitors = params,
                        nchains = nchains, niter = niter, 
                        nburnin = nburn, thin = nthin, 
                        samplesAsCodaMCMC = TRUE, setSeed = mySeed)


##########################################################################################
##########################################################################################
### Running the integrated distance sampling model where we include a known-fate formulation 
### for the survival process. 

## Add fake data for S1 and S2
Survs1 <- cbind(c(16, 28, 29, 29, 36), c(6, 19, 16, 13, 26))
Survs2 <- cbind(c(47, 53, 54, 50, 52), c(34, 32, 35, 42, 33))

nim.data$Survs1 <- Survs1
nim.data$Survs2 <- Survs2

## Set parameters
params3b <- c("esw", "R_year", "p", "S", "D", "S1", "S2")


## Known fate model for S - common S across years (S)

# Function for setting initial values
inits3b <- function(){
  
  mu.D1 <- runif(1, 3, 4)
  N_exp <- matrix(NA, nrow = N_sites, ncol = N_years)
  for(j in 1:N_sites){
    N_exp[j, 1] <- rpois(1, mu.D1*L[j, 1]*(W/nim.constants$scale1)*2)      ## Expected number of birds
  }
  
  list(
    mu.dd = runif(1, 4, 5), 
    sigma.dd = runif(1, 0.05, 2),
    mu.D1 = mu.D1, 
    sigma.D = runif(1, 0.05, 2),
    mu.R = runif(1, -2, 2), 
    sigma.R = runif(1, 0.05, 2),
    eps.dd = rep(0, N_years), 
    eps.R = rep(0, N_years), 
    eps.D1 = rep(0, N_sites),
    S1 = runif(1, 0.6, 0.7), 
    S2 = runif(1, 0.6, 0.7),
    N_exp = N_exp
  )
}

# Sample initial values
initVals3b <- list(inits3b(), inits3b(), inits3b())

# Run code file
source('3_Combined_M3b_KnownFate_nimble.R')

# Test run
out_real3b <- nimbleMCMC(code = modM3b.code, 
                        data = nim.data, constants = nim.constants,
                        inits = initVals3b, monitors = params,
                        nchains = nchains, niter = niter, 
                        nburnin = nburn, thin = nthin, 
                        samplesAsCodaMCMC = TRUE, setSeed = mySeed)







