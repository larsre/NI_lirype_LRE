

##############################
##### Known-fatemodel for annual survival - 
##### based on the binomial
##### Divide the year in two: 1. September - 28. February AND 1. March - 31. August

cat("model{

S1 ~ dunif(0.25,1)
S2 ~ dunif(0.25,1)


for (t in 1:T){
  
  Survs1[t,2] ~ dbinom(S1, Survs1[t,1])
  Survs2[t,2] ~ dbinom(S2, Survs2[t,1])

}


### Derived parameter; 

S <- S1 * S2 

} 
", fill=TRUE,file="KF_Surv_M1.bugs")

#####################################################
#####################################################
library(jagsUI)

Survs1 <- cbind(c(16, 28, 29, 29, 36), c(6, 19, 16, 13, 26))
Survs2 <- cbind(c(47, 53, 54, 50, 52), c(34, 32, 35, 42, 33))




# MCMC settings
niter <- 5000
nthin <- 2
nburn <- 2500
nchains <- 3    

# Setting parameters to monitor            


params <- c("S1", "S2", "S")


### No prior on S - variable S between years (S[t])
jags.dat10 <- list(Survs1=Survs1, Survs2=Survs2, T=4)

inits2 <- function() {list(S1=0.5, S2=0.5)}

out_real_S <- jagsUI::jags(jags.dat10, inits=inits2, params, model.file="KF_Surv_M1.bugs",
                          n.chain=nchains, n.iter=niter, 
                          n.burnin=nburn, parallel = TRUE, DIC=FALSE, 
                          codaOnly=c("Deviance", "Density"))









