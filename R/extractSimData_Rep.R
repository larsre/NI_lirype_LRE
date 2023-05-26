#' Extract reproduction data from simulated distance sampling data
#'
#' @param Jmax integer. The total number of sites/transects surveyed. 
#' @param Tmax integer. The total number of years with distance sampling surveys.
#' @param DS.count array containing counts of adult and juvenile birds from 
#' simulated distance sampling surveys. Dimensions: 1 = age class (1 = juveniles,
#' 2 = adults), 2 = site/transect, 3 = year. 
#' @param R_perF logical. If TRUE, treats recruitment rate as juvenile per adult female.
#' If FALSE, treats recruitment rate as juvenile per adult (sum of both sexes).
#' @param R_parent_drop0 logical. If TRUE, removes observations of juveniles without adults
#' from recruitment data. If FALSE, sets 1 as the number of adults/adults females when none
#' are observed. 
#' 
#' @return a list containing data on juveniles per adult from distance sampling
#' surveys.
#' @export
#'
#' @examples

extractSimData_Rep <- function(Jmax, Tmax, DS.count, R_perF, R_parent_drop0){
  
  # Extract necessary data from DS counts
  sumR_obs <- c(DS.count[1,,])
  
  if(R_perF){
    sumAd_obs <- c(round(DS.count[2,,]*0.5))
  }else{
    sumAd_obs <- c(DS.count[2,,])
  }
  
  sumR_obs_year <- rep(1:Tmax, each = Jmax)
  
  # Discard entries with 0 adults present
  if(R_parent_drop0){
    drop.idx <- which(sumAd_obs == 0)
    sumR_obs <- sumR_obs[-drop.idx]
    sumAd_obs <- sumAd_obs[-drop.idx]
    sumR_obs_year <- sumR_obs_year[-drop.idx]
  }

  # Count observations
  N_sumR_obs <- length(sumR_obs)
  
  # Collate and return data
  Rep.data <- list(sumR_obs = sumR_obs, sumAd_obs = sumAd_obs,
                   sumR_obs_year = sumR_obs_year, N_sumR_obs = N_sumR_obs)
  return(Rep.data)
  
}
