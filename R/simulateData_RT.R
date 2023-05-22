# NOTE: As of now, only one (= annual) survival interval is included in the data
#       simulation. However, extending this to conditional survival over two
#       seasonal periods (as in the real data) is straightforward. 

## Function for simulating known-fate telemetry data
simulateData_RT <- function(nind.avg.RT, Tmin.RT, Tmax.RT, Tmax, SurvProbs){
  
  # Set year range and number
  year_Survs <- c(Tmin.RT:Tmax.RT)
  N_years_RT <- length(year_Survs)
  
  # Make vectors for storing data
  n.rel.S1 <- n.surv.S1 <- rep(NA, N_years_RT)
  n.rel.S2 <- n.surv.S2 <- rep(NA, N_years_RT)
  
  for(t in 1:N_years_RT){
    
    # Sample number of individuals fitted with transmitters each year
    n.rel.S1[t] <- rpois(1, nind.avg.RT)
    n.rel.S2[t] <- rpois(1, nind.avg.RT)
    
    # Simulate survivors to the next year
    n.surv.S1[t] <- rbinom(1, size = n.rel.S1[t], prob = sqrt(SurvProbs[year_Survs[t]]))
    n.surv.S2[t] <- rbinom(1, size = n.rel.S2[t], prob = sqrt(SurvProbs[year_Survs[t]]))
  }
  
  # Assemble data in a list and return
  RT.data <- list(Survs1 = cbind(n.rel.S1, n.surv.S1),
                  Survs2 = cbind(n.rel.S2, n.surv.S2),
                  year_Survs = year_Survs,
                  N_years_RT = N_years_RT)
  return(RT.data)
}
