## Function to simulate year- and site-specific vital rates
simulateVRs <- function(Tmax, Jmax,
                        Mu.S, Mu.R, 
                        sigmaT.S, sigmaT.R, 
                        sigmaJ.S, sigmaJ.R){
  
  # Sample random year effects for all vital rates
  epsilonT.S <- rnorm(Tmax, mean = 0, sd = sigmaT.S)
  epsilonT.R <- rnorm(Tmax, mean = 0, sd = sigmaT.R)
  #epsilonT.sJ <- rnorm(Tmax, mean = 0, sd = sigmaT.sJ)
  
  # Sample random site effects for all vital rates
  epsilonJ.S <- rnorm(Jmax, mean = 0, sd = sigmaJ.S)
  epsilonJ.R <- rnorm(Jmax, mean = 0, sd = sigmaJ.R)
  #epsilonJ.sJ <- rnorm(Jmax, mean = 0, sd = sigmaJ.sJ)
  
  # Calculate year- and site-specific vital rates
  S <- R <- sJ <- matrix(NA, nrow = Jmax, ncol = Tmax)
  for(t in 1:Tmax){
    for(j in 1:Jmax){
      
      S[j,t] <- plogis(qlogis(Mu.S) + epsilonT.S[t] + epsilonJ.S[j])
      R[j,t] <- exp(log(Mu.R) + epsilonT.R[t] + epsilonJ.R[j])
      #sJ[j,t] <- plogis(qlogis(Mu.sJ) + epsilonT.sJ[t] + epsilonJ.sJ[j])
      
    }
  }
  
  # Arrange vital rates in a list
  VR.list <- list(S = S, R = R)
  
}