simulateInits <- function(nim.data, nim.constants){
  
  ## Initial densities / population sizes
  mu.D1 <- runif(nim.constants$N_areas, 3, 4)
  ratio.JA1 <- runif(nim.constants$N_areas, 0.2, 0.6)
  
  N_exp <- array(NA, dim = c(nim.constants$N_areas, nim.constants$N_ageC, max(nim.constants$N_sites), nim.constants$N_years))
  
  for(x in 1:nim.constants$N_areas){
    for(j in 1:nim.constants$N_sites[x]){
      #N_exp1 <- rpois(1, mu.D1[x]*nim.data$L[x, j, 1]*(nim.constants$W/nim.constants$scale1)*2) ## Expected number of birds
      N_exp1 <- extraDistr::rtpois(1, lambda = mu.D1[x]*nim.data$L[x, j, 1]*(nim.constants$W/nim.constants$scale1)*2, a = 2) ## Expected number of birds
      
      N_exp[x, 1, j, 1] <- round(N_exp1*ratio.JA1[x])     
      N_exp[x, 2, j, 1] <- N_exp1 - N_exp[x, 1, j, 1]
    }
  }

  ## Area-specific survival parameters
  h.Mu.S1 <- runif(1, 0.6, 0.7) 
  h.Mu.S2 <- runif(1, 0.6, 0.7)
  h.sigma.S1 <- runif(1, 0, 0.1)
  h.sigma.S2 <- runif(1, 0, 0.1)
  
  mu.S1 <- rnorm(nim.constants$N_areas, plogis(h.Mu.S1), sd = h.sigma.S1)
  mu.S2 <- rnorm(nim.constants$N_areas, plogis(h.Mu.S2), sd = h.sigma.S2)
  
  ## Area-specific reproductive parameters
  h.mu.R  <- runif(1, -5, 5)
  h.sigma.R <- runif(1, 0, 2)
  
  mu.R <- rnorm(nim.constants$N_areas,h.mu.R, sd =  h.sigma.R)
  
  ## Area-specific detection parameters
  h.mu.dd <- runif(1, 2, 5)
  h.sigma.dd <- runif(1, 0, 1)

  #mu.dd <- rnorm(nim.constants$N_areas, h.mu.dd, sd = h.sigma.dd)
  mu.dd <- rep(h.mu.dd, 3)
  
  ## Assemble all in list
  list(
    b = runif(1, 1, 50), 
    mu.D1 = mu.D1, 
    sigma.D = runif(nim.constants$N_areas, 0.05, 2),
    mu.R = mu.R,
    h.mu.R = h.mu.R, h.sigma.R = h.sigma.R,
    mu.dd = mu.dd,
    h.mu.dd = h.mu.dd, h.sigma.dd = h.sigma.dd,
    sigmaT.dd = runif(1, 0.05, 2),
    sigmaT.R = runif(1, 0.05, 2),
    epsT.dd = rep(0, nim.constants$N_years), 
    epsT.R = rep(0, nim.constants$N_years), 
    eps.D1 = matrix(0, nrow = nim.constants$N_areas, ncol = max(nim.constants$N_sites)),
    h.Mu.S1 = h.Mu.S1, h.Mu.S2 = h.Mu.S2,
    h.sigma.S1 = h.sigma.S1, h.sigma.S2 = h.sigma.S2,
    mu.S1 = mu.S1, mu.S2 = mu.S2,
    N_exp = N_exp,
    ratio.JA1 = ratio.JA1
  )
}
