#' Simulate hierarchical distance sampling data
#'
#' @param Jmax integer. Total number of sites to simulate for.
#' @param Tmax integer. Total number of years to simulate for.
#' @param G.age 
#' @param Mu.dd numeric. Average detection decline rate.
#' @param sigmaT.dd numeric. Standard deviation of random year effects on detection.
#' @param sigmaJ.dd numeric. Standard deviation of random site effects on detection.
#' @param W numeric. Truncation distance.
#' @param min.Tlength numeric. Minimum transect length.
#' @param max.Tlength numeric. Maximum transect length.
#' @param discard0 logical. If TRUE, retain only data from observations (drop
#' non-observations). 
#'
#' @return a list containing simulated observation distances ("d") with corresponding
#' years ("d_year") and sites ("d_site"), age-structured count observations ("DS.count")
#' and transect lengths ("L").
#' @export
#'
#' @examples

simulateData_HDS <- function(Jmax, Tmax, G.age,
                             Mu.dd, sigmaT.dd, sigmaJ.dd,
                             W, min.Tlength, max.Tlength,
                             discard0){
  
  # Set site-specific transect lengths
  # NOTE: These are not currently influencing the data simulation in any way
  L <- matrix(runif(Jmax*Tmax, min.Tlength, max.Tlength), nrow = Jmax, ncol = Tmax)
  
  # Sample random year and site effects on sigma parameter
  epsilonT <- rnorm(Tmax, mean = 0, sd = sigmaT.dd)
  epsilonJ <- rnorm(Jmax, mean = 0, sd = sigmaJ.dd)
  
  # Calculate year- and site-specific sigma parameter
  Sigma <- exp(log(Mu.dd) + outer(epsilonJ, epsilonT, FUN = '+'))
  
  # Copy group data
  data <- G.age
  
  # Sample distance from transect line for each group (assuming uniform distribution)
  data$distance <- round(runif(nrow(data), 0, W))
  
  # Calculate group detection probabilities p based on distances d & simulate observation process
  data$p <- NA
  data$obs <- NA
  for(i in 1:nrow(data)){
    data$p[i] <- exp(-data$distance[i] * data$distance[i]/(2 * (Sigma[data$site[i], data$year[i]]^2)))
    data$obs[i] <- rbinom(n = 1, size = 1, prob = data$p[i])
  }
  
  # Retain only data from observations if "discard0"
  if(discard0){
    data <- subset(data, obs == 1)
  }
  
  # Summarise number of individuals observed per year-site combination
  #DS.count <- array(NA, nrow = Jmax, ncol = Tmax)
  DS.count <- array(NA, c(2, Jmax, Tmax))
  for(t in 1:Tmax){
    for(j in 1:Jmax){
      data.sub <- subset(data, year == t & site == j)
      DS.count[1, j, t] <- sum(data.sub$no_juv)
      DS.count[2, j, t] <- sum(data.sub$no_ad)
    }
  }
  
  # Collate and return data
  DS.data <- list(d = data$distance, d_year = data$year, d_site = data$site,
                  DS.count = DS.count, L = L)
  return(DS.data)
}
