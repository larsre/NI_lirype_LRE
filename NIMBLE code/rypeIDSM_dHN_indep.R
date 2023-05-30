rypeIDSM <- nimbleCode({
  

  ###################################################
  ## Random effects model for R (i.e. )
  ## Data: 
  
  # R_obs <- vector with number of chicks / obs [0 - 12]
  # N_Years <- number of years in time series 
  # N_obs <- number of observations 
  # Year_obs <- vector with years for each observation. 
  
  ## Priors; 
  for (t in 1:N_years){
    epsT.R[t] ~ dnorm(0, sd = sigmaT.R)
  }
  
  Mu.R  ~ dunif(0, 10)
  sigmaT.R ~ dunif(0, 5)
  
  ## Constraints;
  R_year[1:N_years] <- exp(log(Mu.R) + epsT.R[1:N_years])
  
  ## Likelihood;
  for (i in 1:N_sumR_obs){
    
    sumR_obs[i] ~ dpois(R_year[sumR_obs_year[i]]*sumAd_obs[i])
  }
  
  
  #####################################################
  ## Model for survival; 
  
  ## Priors
  Mu.S ~ dunif(0, 1)
  Mu.S1 ~ dunif(0.25, 0.9)
  
  ## Constraints
  S[1:N_years] <- Mu.S
  
  S1[1:N_years] <- Mu.S1
  S2[1:N_years] <- S[1:N_years]/S1[1:N_years]
  
  ## Data likelihoods
  for (t in 1:N_years_RT){
    
    Survs1[t, 2] ~ dbinom(S1[year_Survs[t]], Survs1[t, 1])
    Survs2[t, 2] ~ dbinom(S2[year_Survs[t]], Survs2[t, 1])
    
  }
  
  
})
