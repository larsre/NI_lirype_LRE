# NOTE: As of now, only one (= annual) survival interval is included in the data
#       simulation. However, extending this to conditional survival over two
#       seasonal periods (as in the real data) is straightforward. 

## Function for simulating known-fate telemetry data
simulateData_RT <- function(nind.avg.RT, Tmin.RT, Tmax.RT, Tmax, SurvProbs){
  
  # Make vectors for storing data
  n.rel.S1 <- n.surv.S1 <- rep(NA, Tmax)
  n.rel.S2 <- n.surv.S2 <- rep(NA, Tmax)
  
  for(t in Tmin.RT:Tmax.RT){
    
    # Sample number of individuals fitted with transmitters each year
    n.rel.S1[t] <- rpois(1, nind.avg.RT)
    n.rel.S2[t] <- rpois(1, nind.avg.RT)
    
    # Simulate survivors to the next year
    n.surv.S1[t] <- rbinom(1, size = n.rel.S1[t], prob = sqrt(SurvProbs[t]))
    n.surv.S2[t] <- rbinom(1, size = n.rel.S2[t], prob = sqrt(SurvProbs[t]))
  }
  
  # Assemble data in a list and return
  RT.data <- list(Survs1 = cbind(n.rel.S1, n.surv.S1),
                  Survs2 = cbind(n.rel.S2, n.surv.S2))
  return(RT.data)
}
