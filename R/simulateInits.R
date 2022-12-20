#' Simulate values for initializing integrated model
#'
#' @param nim.data list of input objects representing data
#' @param nim.constants list of input objects representing constants
#'
#' @return A list containing one complete set of initial values for the model.
#' @export
#'
#' @examples

simulateInits <- function(nim.data, nim.constants){
  
  # Limits and constants #
  #----------------------#
  
  N_areas <- nim.constants$N_areas
  N_ageC <- nim.constants$N_ageC
  N_years <- nim.constants$N_years
  N_sites <- nim.constants$N_sites
  
  L <- nim.data$L
  W <- nim.constants$W
  pi <- 3.141593
  A <- nim.data$A
  
  
  # Vital rates #
  #-------------#
  
  ## Area-specific survival parameters
  h.Mu.S1 <- runif(1, 0.6, 0.7) 
  h.Mu.S2 <- runif(1, 0.6, 0.7)
  h.sigma.S1 <- runif(1, 0, 0.1)
  h.sigma.S2 <- runif(1, 0, 0.1)
  
  mu.S1 <- rnorm(N_areas, plogis(h.Mu.S1), sd = h.sigma.S1)
  mu.S2 <- rnorm(N_areas, plogis(h.Mu.S2), sd = h.sigma.S2)
  
  Mu.S1 <- Mu.S2 <- rep(NA, N_areas)
  S1 <- S2 <- matrix(NA, nrow = N_areas, ncol = N_years)
  
  for(x in 1:N_areas){
    Mu.S1[x] <- plogis(mu.S1[x])
    Mu.S2[x] <- plogis(mu.S2[x])
    
    S1[x, 1:N_years] <- Mu.S1[x]
    S2[x, 1:N_years] <- Mu.S2[x]
  }
  S <- S1*S2
  
  ## Area-specific reproductive parameters
  h.mu.R  <- runif(1, 0, 0.5)
  h.sigma.R <- runif(1, 0, 0.05)
  
  mu.R <- rnorm(N_areas, h.mu.R, sd =  h.sigma.R)
  sigmaT.R <- runif(1, 0, 2)
  
  epsT.R <- rep(0, N_years)
  #epsT.T <- rnorm(N_year, 0, sigmaT.R)
  
  R_year <- matrix(NA, nrow = N_areas, ncol = N_years)
  
  for(x in 1:N_areas){
    R_year[x, 1:N_years] <- exp(mu.R[x] + epsT.R[1:N_years])
  }
  
  
  # Detection parameters #
  #----------------------#
  
  ## Area-specific detection parameters
  h.mu.dd <- runif(1, 2, 5)
  h.sigma.dd <- runif(1, 0, 1)
  
  #mu.dd <- rnorm(N_areas, h.mu.dd, sd = h.sigma.dd)
  mu.dd <- rep(h.mu.dd, N_areas)
  sigmaT.dd <- runif(1, 1, 10)
  
  epsT.dd <- rep(0, N_years)
  #epsT.dd <- rnorm(N_years, 0, sd = sigmaT.dd)
  
  sigma <- esw <- p <- matrix(NA, nrow = N_areas, ncol = N_years)
  
  for(x in 1:N_areas){
    for(t in 1:N_years){
      sigma[x, t] <- exp(mu.dd[x] + epsT.dd[t])
      esw[x, t] <- sqrt(pi * sigma[x, t]^2 / 2) 
      p[x, t] <- min(esw[x, t], W) / W
    }
  }
  
  
  # Population model #
  #------------------#
  
  ## Initial densities / population sizes
  mu.D1 <- rep(NA, N_areas)
  sigma.D <- runif(N_areas, 0.05, 2)
  ratio.JA1 <- runif(N_areas, 0.2, 0.6)
  
  N_exp <- Density <- array(0, dim = c(N_areas, N_ageC, max(nim.constants$N_sites), N_years))
  
  for(x in 1:N_areas){
    
    D_x_sum <- apply(nim.data$N_a_line_year[x,,,], c(2,3), sum) / (L[x,,]*W*2)
    D_data <- D_x_sum[which(!is.na(D_x_sum) & D_x_sum > 0)]
    mu.D1[x] <- runif(1, quantile(D_data, 0.25), quantile(D_data, 0.75))  
    
    for(j in 1:N_sites[x]){
      
      Density[x, 1, j, 1] <- mu.D1[x]*ratio.JA1[x]
      Density[x, 2, j, 1] <- mu.D1[x]*(1-ratio.JA1[x])
      
      N_exp[x, 1, j, 1] <- extraDistr::rtpois(1, lambda = Density[x, 1, j, 1]*L[x, j, 1]*W*2, a = 1)
      N_exp[x, 2, j, 1] <- extraDistr::rtpois(1, lambda = Density[x, 2, j, 1]*L[x, j, 1]*W*2, a = 1)
    }
  }

  ## Population projection over time
  for(x in 1:N_areas){
    for(j in 1:N_sites[x]){
      for(t in 2:N_years){
        
        Density[x, 2, j, t] <- sum(Density[x, 1:N_ageC, j, t-1])*S[x, t-1] # Adults
        Density[x, 1, j, t] <- Density[x, 2, j, t]*R_year[x, t]/2 # Juveniles
        
        N_exp[x, 1:N_ageC, j, t] <- Density[x, 1:N_ageC, j, t]*L[x, j, t]*W*2
      }
    }
  }
  
  ## Area-specific population size and density
  N_tot_exp <- matrix(NA, nrow = N_areas, ncol = N_years)
  
  for(x in 1:N_areas){
    for (t in 1:N_years){
      N_tot_exp[x, t] <- sum(N_exp[x, 1, 1:N_sites[x], t] + N_exp[x, 2, 1:N_sites[x], t])    ## Summing up expected number of birds in covered area; 
    }
  }
  
  # Assembly #
  #----------#
  
  list(
    #b = runif(1, 1, 50), 
    
    mu.D1 = mu.D1, 
    sigma.D = sigma.D,
    eps.D1 = matrix(0, nrow = nim.constants$N_areas, ncol = max(nim.constants$N_sites)),
    ratio.JA1 = ratio.JA1,
    
    mu.R = mu.R,
    h.mu.R = h.mu.R, h.sigma.R = h.sigma.R,
    sigmaT.R = sigmaT.R,
    epsT.R = epsT.R, 
    R_year = R_year,
    
    mu.dd = mu.dd,
    h.mu.dd = h.mu.dd, h.sigma.dd = h.sigma.dd,
    sigmaT.dd = sigmaT.dd,
    epsT.dd = epsT.dd,
    sigma = sigma, sigma2 = sigma^2,
    esw = esw,
    p = p,
    
    h.Mu.S1 = h.Mu.S1, h.Mu.S2 = h.Mu.S2,
    h.sigma.S1 = h.sigma.S1, h.sigma.S2 = h.sigma.S2,
    mu.S1 = mu.S1, mu.S2 = mu.S2,
    S1 = S1, S2 = S2, S = S,

    Density = Density,
    N_exp = N_exp,
    N_tot_exp = N_tot_exp
  )
}
