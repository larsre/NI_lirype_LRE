

library(R2jags)

##### Running the models; 
# MCMC settings
niter <- 25000
nthin <- 2
nburn <- 10000
nchains <- 3    



################################################################################
################################################################################
### No prior on S - variable S between years (S[t])

# Setting parameters to monitor            
params <- c("esw", "R_year", "p", "S", "D")


jags.dat2 <- list(R_obs=R_obs, y=y, N_years=N_years, W=W, scale1=1000,
                  N_obs=N_obs, Year_obs=Year_obs, 
                  zeros.dist=zeros_dist, L=L, 
                  N_sites=N_sites, N_line_year=N_line_year, A=A, 
                  R_obs_year=R_obs_year, N_R_obs=N_R_obs)

inits2 <- function() {list(mu.dd=runif(1, 4,5), mu.D1=runif(1, 3,4), S=runif(N_years, 0.4, 0.5))}

out_real2 <- jagsUI::jags(jags.dat2, inits=inits2, params, model.file="Combined_M2.bugs",
                          n.chain=nchains, n.iter=niter, 
                          n.burnin=nburn, parallel = TRUE, DIC=FALSE, 
                          codaOnly=c("Deviance", "Density"))

################################################################################
################################################################################
## No prior on S - common S across years (S)

# Setting parameters to monitor            
params <- c("esw", "R_year", "p", "S", "D")

jags.dat3 <- list(R_obs=R_obs, y=y, N_years=N_years, W=W, scale1=1000,
                  N_obs=N_obs, Year_obs=Year_obs, 
                  zeros.dist=zeros_dist, L=L, 
                  N_sites=N_sites, N_line_year=N_line_year, A=A, 
                  R_obs_year=R_obs_year, N_R_obs=N_R_obs)

inits3 <- function() {list(mu.dd=runif(1, 4,5), mu.D1=runif(1, 2,3), S=runif(1, 0.4, 0.5))}

out_real3 <- jagsUI::jags(jags.dat3, inits=inits3, params, model.file="Combined_M3.bugs",
                          n.chain=nchains, n.iter=niter, 
                          n.burnin=nburn, parallel = TRUE, DIC=FALSE, 
                          codaOnly=c("Deviance", "Density"))


################################################################################
################################################################################
## ## Prior on S based on radiotelemetry study;
shape_from_stats <- function(mu , sigma ){
  a <-(mu^2-mu^3-mu*sigma^2)/sigma^2
  b <- (mu-2*mu^2+mu^3-sigma^2+mu*sigma^2)/sigma^2
  shape_ps <- c(a,b)
  return(shape_ps)
}

S_priors <- shape_from_stats(0.425, 0.035)

## Plott - to see if this works; 

x=seq(0,1,.001)
w=dbeta(x,S_priors[1],S_priors[2])
plot(x,w,typ="l")



jags.dat4 <- list(R_obs=R_obs, y=y, N_years=N_years, W=W, scale1=1000,
                  N_obs=N_obs, Year_obs=Year_obs, 
                  zeros.dist=zeros_dist, L=L, 
                  N_sites=N_sites, N_line_year=N_line_year, A=A, 
                  R_obs_year=R_obs_year, N_R_obs=N_R_obs, a=S_priors[1], b=S_priors[2])

inits4 <- function() {list(mu.dd=runif(1, 4,5), mu.D1=runif(1, 2,3))}

out_real4 <- jagsUI::jags(jags.dat4, inits=inits4, params, model.file="Combined_M4.bugs",
                          n.chain=nchains, n.iter=niter, 
                          n.burnin=nburn, parallel = TRUE, DIC=FALSE, 
                          codaOnly=c("Deviance", "Density"))




##########################################################################################
##########################################################################################
### 


