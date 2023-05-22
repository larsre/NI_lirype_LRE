## Function to extract reproductive data from distance sampling data
extractSimData_Rep <- function(Jmax, Tmax, DS.count){
  
  # Extract necessary data from DS counts
  sumR_obs <- c(DS.count[1,,])
  sumAd_obs <- c(DS.count[2,,])
  sumR_obs_year <- rep(1:Tmax, each = Jmax)
  
  # Discard entries with 0 adults present
  drop.idx <- which(sumAd_obs == 0)
  sumR_obs <- sumR_obs[-drop.idx]
  sumAd_obs <- sumAd_obs[-drop.idx]
  sumR_obs_year <- sumR_obs_year[-drop.idx]
  
  # Count observations
  N_sumR_obs <- length(sumR_obs)
  
  # Collate and return data
  Rep.data <- list(sumR_obs = sumR_obs, sumAd_obs = sumAd_obs,
                   sumR_obs_year = sumR_obs_year, N_sumR_obs = N_sumR_obs)
  return(Rep.data)
  
}
