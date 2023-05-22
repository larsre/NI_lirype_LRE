

#' Title
#'
#' @param type The type of parameter to compare. Either "demographic", "density" or "detection" 
#' @param simData List object. Output from data simulation model
#' @param mcmcData mcmc.list object. output from NIMLBE model run.  
#' @return Plot comparing posterior distributions with parameters used to simulate data 
#' @export
#'
#' @examples



plotSimCheck <- function (type = "demographic", simData = d1, mcmcData = d2) {

  require(coda)
  require(tidyverse)
  require(tidybayes)
  require(ggforce)
  require(see)
  require(cowplot)
  
################################################################################    
if(type == "demographic") {
  ## Prepare tidy data
  R_year <- mcmcData %>% spread_draws(R_year[year])
  Mu_R <- mcmcData %>% spread_draws(mu.R) %>% mutate(lab_code = "mu.R")
  sigmaT_R <- mcmcData %>% spread_draws(sigma.R) %>% mutate(lab_code = "sigma.R")
  
  Mu_S1 <- mcmcData %>% spread_draws(Mu.S1) %>% mutate(Surv = "S1") %>% rename(S = Mu.S1) %>% select(S, Surv)
  Mu_S2 <- mcmcData %>% spread_draws(Mu.S2) %>% mutate(Surv = "S2") %>% rename(S = Mu.S2) %>% select(S, Surv)
  Mu_S <-  tibble(S = Mu_S1$S * Mu_S2$S, Surv = "S") %>% bind_rows(., Mu_S1, Mu_S2)
  
   # surival probabilities
  p1 <- ggplot(data = Mu_S, aes(x = Surv, y = S)) +
    geom_violinhalf(fill = "lightgreen", alpha = 0.7)  +
    ylim(0, 1) +
    geom_point(aes(x = "S", y = simData$SimParams$Mu.S), 
               col = "darkblue", size = 3) +
    theme_cowplot() +
    theme(panel.grid.major.y = element_line(color = "#8ccde3",
                                          size = 0.25,
                                          linetype = 2), 
          text = element_text(size = 10)) +
    ylab("Survival probability") +
    xlab("Survival component")
  
    
  
  # mean recruitment (mu.R)
  p2 <- ggplot(data = Mu_R, aes(x = lab_code, y = exp(mu.R))) +
    geom_violinhalf(fill = "lightgreen", alpha = 0.7)  +
    ylim(0, 4) +
    geom_point(aes(x = "mu.R", y = simData$SimParams$Mu.R), 
               col = "darkblue", size = 3) +
    theme_cowplot() +
    theme(panel.grid.major.y = element_line(color = "#8ccde3",
                                            size = 0.25,
                                            linetype = 2), 
          text = element_text(size = 10)) +
    ylab("Recruitment rate") +
    xlab("")
  
  # tempral variation in R (sigma.R)
  p3 <- ggplot(data = sigmaT_R, aes(x = lab_code, y = sigma.R)) +
    geom_violinhalf(fill = "lightgreen", alpha = 0.7)  +
    ylim(0, max(sigmaT_R$sigma.R)+1) +
    geom_point(aes(x = "sigma.R", y = simData$SimParams$sigmaT.R), 
               col = "darkblue", size = 3) +
    theme_cowplot() +
    theme(panel.grid.major.y = element_line(color = "#8ccde3",
                                            size = 0.25,
                                            linetype = 2), 
          text = element_text(size = 10)) +
    ylab("sigma R") +
    xlab("")
  
    
  ## Annual recruitment rates in simulated data
  Na_temp <- apply(simData$N.data, c(2,3), sum) 
  Na_temp2 <- as_tibble(t(Na_temp)) %>%
    mutate(R = V1 / V2, year = seq(1:simData$SimParams$Tmax))
  
  p4 <- ggplot(data = R_year, aes(x = as.factor(year), y = R_year)) +
    geom_violinhalf(fill = "lightgreen", alpha = 0.7)  +
    ylim(0, 5) +
    geom_point(data=Na_temp2, aes(x = as.factor(year), y = R), 
               col = "darkblue", size = 3) +
    theme_cowplot() +
    theme(panel.grid.major.y = element_line(color = "#8ccde3",
                                            size = 0.25,
                                            linetype = 2), 
          text = element_text(size = 10)) +
    ylab("R") +
    xlab("Year")

    ## Putting the plots together in a multipanel plot
    p_a <- plot_grid(p1, p2, p3, nrow = 1)
    p_out <- plot_grid(p_a, p4, nrow = 2, label_size = 10)
  
    }

################################################################################
  if(type == "density"){
    
    # Prepare density (for each time step) for simulated data
    N_temp <- apply(simData$N.data, 3, sum) 
    A_temp <- apply(simData$DS.data$L, 2, sum) * simData$SimParams$W*2 / (1000 *1000)
    D_temp <- tibble(year=seq(1:simData$SimParams$Tmax), Mean.D = N_temp / A_temp) 
    
    # Summarize mcmc data - annual density. This could (should) be done with 
    # spread_draws()
    
    Density_year <- mcmcData %>% spread_draws(D[year]) %>% mutate(density = D*1000) 
     
    p_out <- ggplot(data = Density_year, aes(x = as.factor(year), y = density)) +
      geom_violinhalf(fill = "lightgreen", alpha = 0.7)  +
      ylim(0, 160) +
      geom_point(data=D_temp, aes(x = as.factor(year), y = Mean.D), 
                 col = "darkblue", size = 3) +
      theme_cowplot() +
      theme(panel.grid.major.y = element_line(color = "#8ccde3",
                                              size = 0.25,
                                              linetype = 2), 
            text = element_text(size = 10)) +
      ylab("Density") +
      xlab("Year")
    

  }  

if(type == "detection"){

  # mean half normal scale parameter (mu.dd)
  Mu_dd <- mcmcData %>% spread_draws(mu.dd) %>% mutate(lab_code = "mu.dd")
  p1 <- ggplot(data = Mu_dd, aes(x = lab_code, y = exp(mu.dd))) +
    geom_violinhalf(fill = "lightgreen", alpha = 0.7)  +
    #ylim(0, 4) +
    geom_point(aes(x = "mu.dd", y = simData$SimParams$Mu.dd ), 
               col = "darkblue", size = 3) +
    theme_cowplot() +
    theme(panel.grid.major.y = element_line(color = "#8ccde3",
                                            size = 0.25,
                                            linetype = 2), 
          text = element_text(size = 10)) +
    ylab("HN scale para") +
    xlab("")
  
  ## temporal variation in scale parameter (sigma.dd) - i.e. sd of random var.
  sigmaT_dd <- mcmcData %>% spread_draws(sigma.dd) %>% mutate(lab_code = "sigma.dd")
  
  p2 <- ggplot(data = sigmaT_dd, aes(x = lab_code, y = sigma.dd)) +
    geom_violinhalf(fill = "lightgreen", alpha = 0.7)  +
    #ylim(0, 4) +
    geom_point(aes(x = "sigma.dd", y = simData$SimParams$sigmaT.dd), 
               col = "darkblue", size = 3) +
    theme_cowplot() +
    theme(panel.grid.major.y = element_line(color = "#8ccde3",
                                            size = 0.25,
                                            linetype = 2), 
          text = element_text(size = 10)) +
    ylab("sigmaT.dd") +
    xlab("")
  
  ## esw across years
  esw_year <- mcmcData %>% spread_draws(esw[year])
  
  
  p4 <- ggplot(data = esw_year, aes(x = as.factor(year), y = esw)) +
    geom_violinhalf(fill = "lightgreen", alpha = 0.7)  +
   # ylim(0, 160) +
    theme_cowplot() +
    theme(panel.grid.major.y = element_line(color = "#8ccde3",
                                            size = 0.25,
                                            linetype = 2), 
          text = element_text(size = 10)) +
    ylab("Esw") +
    xlab("Year")
  
  ## esw across years
  p_year <- mcmcData %>% spread_draws(p[year])
  
  p5 <- ggplot(data = p_year, aes(x = as.factor(year), y = p)) +
    geom_violinhalf(fill = "lightgreen", alpha = 0.7)  +
    # ylim(0, 160) +
    theme_cowplot() +
    theme(panel.grid.major.y = element_line(color = "#8ccde3",
                                            size = 0.25,
                                            linetype = 2), 
          text = element_text(size = 10)) +
    ylab("Detection probability") +
    xlab("Year")
  
  ########################
  p3 <- plot_grid(p1, p2, nrow = 1)
  p_out <- plot_grid(p3, p4, p5, nrow=3)
  
  
}  
  
  p_out  
}