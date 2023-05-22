## Function for simulating temporal population dynamics in each site
simulatePopDyn <- function(Amax, Tmax, Jmax, VR.list, stochastic = TRUE, plot = FALSE){
  
  # Prepare arrays for storing population projections
  N <- array(NA, dim = c(Jmax, Amax, Tmax))
  #Chicks <- matrix(NA, nrow = Jmax, ncol = Tmax)
  
  # Initialize population projection
  for(j in 1:Jmax){
    N[j,2,1] <- round(runif(1, N1_juv_limits[1], N1_juv_limits[2])) # Number of adults
    N[j,1,1] <- rpois(1, N[j,2,1]*VR.list$R[1]) # Number of juveniles
  }
  
  # Print simulation info
  if(stochastic){
    message('Simulating population dynamics with demographic stochasticity...')
  }else{
    message('Simulating deterministic population dynamics...')
  }
  
  # Simulate survival & reproduction over time in each site
  for(j in 1:Jmax){
    for(t in 1:(Tmax-1)){
      
      if(stochastic){
        
        # - Survivors
        N[j,2,t+1] <- rbinom(1, size = sum(N[j,1:Amax,t]), prob = VR.list$S[j,t])
        
        # - Reproduction of survivors
        N[j,1,t+1] <- rpois(1, lambda = N[j,2,t+1]*(VR.list$R[j,t+1]))
        
        # - Young-of-the-year recruiting in the end of summer
        # N[j,1,t+1] <- rbinom(1, size = Chicks[j,t+1], prob = VR.list$sJ[j,t+1])
        
        
      }else{
        
        # - Survivors
        N[j,2,t+1] <- sum(N[j,1:Amax,t])*VR.list$S[j,t]
        
        # - Reproduction of survivors
        N[j,1,t+1] <- N[j,2,t+1]*(VR.list$R[j,t+1])
        
        # - Young-of-the-year recruiting in the end of summer
        #N[j,1,t+1] <- Chicks[j,t+1]*VR.list$sJ[j,t+1]
      }
    }
  }
  
  # Plot trajectories
  #par(mfrow = c(1,2))
  plot(colSums(N[,1,]), type = 'l', col = 'red', lty = 'dashed', 
       ylab = 'Number', xlab = 'Time', ylim = c(0, max(colSums(N[,1,]), na.rm = T)))
  lines(colSums(N[,2,]), col = 'blue')
  legend('topleft', legend = c('Recruits', 'Adults'), 
         col = c('red', 'blue'),
         lty = c('dashed', 'solid'), 
         cex = 0.8, y.intersp = 0.1, bty = 'n')
  
  matplot(apply(N, 1, colSums), type = 'l', ylab = 'Number per site', xlab = 'Time')
  
  # Arrange simulated numbers in a list & return
  SimData <- list(N = N)
  return(SimData)
}
