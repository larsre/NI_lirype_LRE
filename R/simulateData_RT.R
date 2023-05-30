#' Simulate known-fate telemetry data
#' This functions simulates a dataset of known-fate telemetry data for two 
#' seasons (S1 and S2). At present, it assumes that survival is equal for both
#' seasons and equal to the square-root of annual survival. It would be 
#' straightforward to extend this to allow for different survival in the two 
#' seasons. 
#' @param nind.avg.RT integer. Average number of individuals fitted with transmitters
#' every season.
#' @param Tmin.RT integer. Index of the first year for which to simulate telemetry
#' data.  
#' @param Tmax.RT integer. Index of the last year for which to simulate telemetry 
#' data.
#' @param SurvProbs vector of time-dependent survival probabilities. 
#'
#' @return a list containing numbers of individuals with transmitters and subset 
#' of survivors for the first (Survs1) and second (Survs2) seasons, as well as
#' corresponding year indices.  
#' @export
#'
#' @examples

simulateData_RT <- function(nind.avg.RT, Tmin.RT, Tmax.RT, SurvProbs){
  
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
