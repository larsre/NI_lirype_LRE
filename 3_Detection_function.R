
cat("model{
# priors for distance model; 

pi <- 3.141593

#random year effect for distance sampling model; 
## Priors for hyper-parameters
mu.dd ~ dunif(0,30)
tau.dd <- pow(sigma.dd, -2)
sigma.dd ~ dunif(0, 10)

for(t in 1:N_years){
  eps.dd[t] ~dnorm(0, tau.dd) 
}



########################################################
for(t in 1:N_years){
  sigma[t] <- exp(mu.dd + eps.dd[t])
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





}
", fill=TRUE,file="obs_model.bugs")


#####################################################
### 
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

jags.dat <- list(y=y*1000, N_years=N_years, N_obs=N_obs, Year_obs=Year_obs, 
                 W=200, zeros.dist=zeros_dist)

# MCMC settings
niter <- 100000
nthin <- 2
nburn <- 50000
nchains <- 3    

# Setting parameters to monitor            

#params <- c("b.df.0")

params <- c("esw", "sigma")

inits1 <- function() {list(mu.dd=runif(0,1))}

out1 <- jagsUI::jags(jags.dat, inits=inits1, params, model.file="obs_model.bugs",
                     n.chain=nchains, n.iter=niter, 
                     n.burnin=nburn, parallel = TRUE, DIC=FALSE, 
                     codaOnly=c("Deviance", "Density"))






