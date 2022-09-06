library(nimble)
library(magrittr)

#---------------#
# GENERAL SETUP #
#---------------#

## Set seed
mySeed <- 0
set.seed(mySeed)

## Load line transect data
rype.data <- readRDS("RypeData_forIM.rds")

## Load and format known fate CMR data
CMR_data <- tibble::as_tibble(read.csv("Demographic_data/CMR_Data.csv", header=T, sep=",")) 

Survs1 <- CMR_data %>% dplyr::filter(TimePeriod == 1) %>% dplyr::select(-YearPeriod, -TimePeriod) %>%
  as.matrix()
Survs2 <- CMR_data %>% dplyr::filter(TimePeriod == 2) %>% dplyr::select(-YearPeriod, -TimePeriod) %>%
  as.matrix()

## Assembling data and constants for all models
nim.data <- list(R_obs = rype.data$R_obs, y = rype.data$y, 
                 zeros.dist = rype.data$zeros_dist, L = rype.data$L, 
                 N_line_year = rype.data$N_line_year, 
                 N_a_line_year = rype.data$N_a_line_year, 
                 A = rype.data$A,
                 Survs1 = Survs1, Survs2 = Survs2)

nim.constants <- list(N_years = rype.data$N_years, W = rype.data$W, scale1 = 1000,
                     N_obs = rype.data$N_obs, Year_obs = rype.data$Year_obs,
                     N_sites = rype.data$N_sites, 
                     R_obs_year = rype.data$R_obs_year, N_R_obs = rype.data$N_R_obs,
                     N_ageC = 2)


#-------------#
# MODEL CODES #
#-------------#

## Model A - Original process model formulation
source('NIMBLE code/3_Combined_M3b_KnownFate_nimble.R')

## Model B - Altered process model without additional data
## Model C - Altered process model + explicit use of age-structure in data
source('NIMBLE code/3_Combined_M3b_KnownFate_nimble_Alt.R')


#-----------------------------#
# INITIALIZATION & MCMC SETUP #
#-----------------------------#

## Set parameters to monitor
params3b <- c("esw", "R_year", "p", "S", "D", "S1", "S2", 
              "Density", "N_exp", "mu.D1")

## Function for setting initial values (altered process model)
inits3b <- function(){
  
  mu.D1 <- runif(1, 3, 4)
  ratio.JA1 <- runif(1, 0.2, 0.6)
  
  N_exp <- array(NA, dim = c(nim.constants$N_ageC, nim.constants$N_sites, nim.constants$N_years))
  
  for(j in 1:nim.constants$N_sites){
    N_exp1 <- rpois(1, mu.D1*nim.data$L[j, 1]*(nim.constants$W/nim.constants$scale1)*2) ## Expected number of birds
    
    N_exp[1, j, 1] <- round(N_exp1*ratio.JA1)     
    N_exp[2, j, 1] <- N_exp1 - N_exp[1, j, 1]
  }
  
  list(
    mu.dd = runif(1, 4, 5), 
    sigma.dd = runif(1, 0.05, 2),
    mu.D1 = mu.D1, 
    sigma.D = runif(1, 0.05, 2),
    mu.R = runif(1, -2, 2), 
    sigma.R = runif(1, 0.05, 2),
    eps.dd = rep(0, nim.constants$N_years), 
    eps.R = rep(0, nim.constants$N_years), 
    eps.D1 = rep(0, nim.constants$N_sites),
    Mu.S1 = runif(1, 0.6, 0.7), 
    Mu.S2 = runif(1, 0.6, 0.7),
    N_exp = N_exp,
    ratio.JA1 = ratio.JA1
  )
}

## Sample initial values (altered process model)
initVals3b.alt <- list(inits3b(), inits3b(), inits3b())

## Summarize nodes for original process model
initVals3b.orig <- initVals3b.alt
for(i in 1:3){
  initVals3b.orig[[i]]$N_exp <- apply(initVals3b.alt[[i]]$N_exp, c(2,3), sum)
}

## Determine if this is a testrun or not
testRun <- TRUE
#testRun <- FALSE

## Set MCMC parameters
if(testRun){
  niter <- 2
  nthin <- 1
  nburn <- 0
  nchains <- 3  
}else{
  niter <- 25000
  nthin <- 5
  nburn <- 5000
  nchains <- 3  
}


#-----------#
# TEST RUNS #
#-----------#

## Test run - Model A
out_modA <- nimbleMCMC(code = modM3b.code.A, 
                        data = nim.data, constants = nim.constants,
                        inits = initVals3b.orig, monitors = params3b,
                        nchains = nchains, niter = niter, 
                        nburnin = nburn, thin = nthin, 
                        samplesAsCodaMCMC = TRUE, setSeed = mySeed)

saveRDS(out_modA, file = 'mod3b_versionA_realData.rds')


## Test run - Model B
out_modB <- nimbleMCMC(code = modM3b.code.B, 
                       data = nim.data, constants = nim.constants,
                       inits = initVals3b.alt, monitors = params3b,
                       nchains = nchains, niter = niter, 
                       nburnin = nburn, thin = nthin, 
                       samplesAsCodaMCMC = TRUE, setSeed = mySeed)

saveRDS(out_modB, file = 'mod3b_versionB_realData.rds')


## Test run - Model C
out_modC <- nimbleMCMC(code = modM3b.code.C, 
                       data = nim.data, constants = nim.constants,
                       inits = initVals3b.alt, monitors = params3b,
                       nchains = nchains, niter = niter, 
                       nburnin = nburn, thin = nthin, 
                       samplesAsCodaMCMC = TRUE, setSeed = mySeed)

saveRDS(out_modC, file = 'mod3b_versionC_realData.rds')


