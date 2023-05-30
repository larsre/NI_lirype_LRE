#' Simulate group aggregation of a population
#'
#' @param Jmax integer. Total number of sites to simulate for. 
#' @param Tmax integer. Total number of years to simulate for.
#' @param N.age array of simulated population sizes with dimension [j,a,t], 
#' where j = site, a = age class, t = year.
#' @param avg_Gsize numeric. Average group size. 
#' @param discard0 logical. If TRUE, removes 0 groups (units that were not 
#' asssigned any individuals). 
#'
#' @return
#' @export a dataframe containing information on composition of simulated groups
#' within each year and site. 
#'
#' @examples

simulateGroups <- function(Jmax, Tmax, N.age, avg_Gsize, discard0 = TRUE){
  
  # Set up group data frame and matrix
  G.age <- data.frame()
  G.no <- matrix(0, nrow = Jmax, ncol = Tmax)
  
  for(j in 1:Jmax){
    for(t in 1:Tmax){
      
      if(sum(N.age[j, 1:2, t]) > 0){
        # Set number of groups to simulate
        G.noTot <- sum(N.age[j, 1:2, t])/avg_Gsize
        G.noTot <- ifelse(G.noTot < 1, 1, round(G.noTot))
        
        # Distribute juveniles among groups
        G.juv <- rmultinom(1, size = N.age[j, 1, t], prob = rep(1/G.noTot, G.noTot))
        
        # Distribute adults among groups
        G.ad <- rmultinom(1, size = N.age[j, 2, t], prob = rep(1/G.noTot, G.noTot))
        
        # Remove any potential 0 observations
        if(discard0){
          obs0 <- which(G.juv + G.ad == 0)
          if(length(obs0) > 0){
            G.juv <- G.juv[-obs0]
            G.ad <- G.ad[-obs0]
          }
        }
        
        
        # Collate in data frame
        G.data <- data.frame(site = j,
                             year = t,
                             no_juv = G.juv,
                             no_ad = G.ad)
        G.age <- rbind(G.age, G.data)
        
        # Set number of groups for site-year
        G.no[j, t] <- nrow(G.data)
      }
    }
  }    
  
  return(G.age)
}