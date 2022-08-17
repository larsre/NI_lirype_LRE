
cat("model{
####################################################
#### Distance sampling half normal detection function; 
## Data; 

# y <- vector with distances to transect line; 
# N_years <- number of years
# N_obs <- number of observations
# Year_obs <- vector with year for each observation


# priors for distance model; 

pi <- 3.141593

#random year effect for distance sampling model; 
## Priors for hyper-parameters
mu.dd ~ dunif(0,3)
tau.dd <- pow(sigma.dd, -2)
sigma.dd ~ dunif(0, 10)

for(t in 1:N_years){
  eps.dd[t] ~dnorm(0, tau.dd) 
}



########################################################
for(t in 1:N_years){
  sigma[t] <- mu.dd + eps.dd[t]
  sigma2[t] <- sigma[t]*sigma[t]
  
  # effective strip width
  esw[t] <- sqrt(pi * sigma2[t] / 2) 
  #f0[t] <- 1/esw[t] #assuming all detected on the line
}

########################################################   
for (i in 1:N_obs){ 
  # LIKELIHOOD
  # using zeros trick
  y[i] ~ dunif(0,W) 
  L.f0[i] <- exp(-y[i]*y[i] / (2*sigma2[Year_obs[i]])) * 1/esw[Year_obs[i]] #y are the distances
  nlogL.f0[i] <-  -log(L.f0[i])
  zeros.dist[i] ~ dpois(nlogL.f0[i])
}

## Hierarchicla node: esw[t]; 
##
###################################################
## Random effects model for R (i.e. )
## Data: 

# R_obs <- vector with number of chicks / obs [0 - 12]
# N_Years <- number of years in time series 
# N_obs <- number of observations 
# Year_obs <- vector with years for each observation. 

## Priors; 
for (t in 1:N_years){
  eps.R[t] ~ dnorm(0, tau.R)
}

mu.R  ~ dunif(-3,3)
tau.R <- pow(sigma.R, -2)
sigma.R ~ dunif(0, 10)

## Likelihood;
for (i in 1:N_obs){
  
  R_obs[i] ~ dpois(R_exp[i])
  # Linear predictor
  R_exp[i] <- exp(mu.R + eps.R[Year_obs[i]])   
}

### Derive parameter: annual predictions for R
for (t in 1:N_years){
  R_year[t] <- exp(mu.R + eps.R[t])
}

########################################################################
### MODEL FOR Density in year 1:
### Simple random effects model 
# Data; 

# N_sites <- number of sites
# N_line_year <- number of birds pr. line 
# L <- length of transect lines


## Priors; 

for(j in 1:N_sites){
  eps.D1[j] ~dnorm(0, tau.D1)
}

mu.D1 ~dunif(-30, 30)
tau.D1 <- pow(sigma.D, -2)
sigma.D ~ dunif(0, 20)

 p <- esw[1] / W

# State model
for (j in 1:N_sites){
  N_exp[j,1] ~ dpois(Density[j,1]*L[j,1]*W*2)      ## Expected number of birds
  Density[j,1] <- exp(mu.D1 + eps.D1[j])             ## random effects model for spatial variation in density for year 1
  ## Detection model year 1
  N_line_year[j,1] ~ dbin(p, N_exp[j,1])
}
} 
", fill=TRUE,file="obs_R_D1_model.bugs")

##########################################################
##########################################################

library(jagsUI)
library(R2jags)

# y <- vector with distances to transect line; 
# N_years <- number of years
# N_obs <- number of observations
# Year_obs <- vector with year for each observation

y <- abs(rnorm(100, mean=0, sd=10))
N_years <- 5
N_obs <- length(y)
Year_obs <- rep(1:5, each=20)
zeros_dist <- rep(0, length(y))
R_obs <- rpois(100, 4)
L <- matrix(ncol=1, nrow=10, rep(c(1000), each=10))
N_sites <- 10
N_line_year <- matrix(ncol=1, nrow=10, rpois(10, 100))

jags.dat2 <- list(R_obs=R_obs, y=y, N_years=N_years, 
                  N_obs=N_obs, Year_obs=Year_obs, W=20, 
                  zeros.dist=zeros_dist, L=L, 
                  N_sites=N_sites, N_line_year=N_line_year)

# MCMC settings
niter <- 50000
nthin <- 2
nburn <- 25000
nchains <- 3    

# Setting parameters to monitor            

#params <- c("b.df.0")

params <- c("esw", "sigma", "eps.R", "R_year", "Density", "N_exp", "p")

inits2 <- function() {list(mu.dd=runif(0,1))}

out2 <- jagsUI::jags(jags.dat2, inits=inits1, params, model.file="obs_R_D1_model.bugs",
                     n.chain=nchains, n.iter=niter, 
                     n.burnin=nburn, parallel = TRUE, DIC=FALSE, 
                     codaOnly=c("Deviance", "Density"))



