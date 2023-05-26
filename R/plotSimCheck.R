#' Visualize comparison of model estimates with true parameter values from data simulation
#'
#' @param SimData a list containing all data necessary for running the model, as well
#' as all simulation parameters and true quantitites. Output of assembleSimData().
#' @param mcmc.out an mcmc list containing posterior samples from a model run. 
#' @param VitalRates logical. If TRUE (default), plots traces and posterior densities for vital rate parameters.
#' @param DetectParams logical. If TRUE (default), plots traces and posterior densities for detection parameters.
#' @param PopSizes logical. If TRUE (default), plots traces and posterior densities for population sizes.
#' @param Densities logical. If TRUE (default), plots traces and posterior densities for population densities.
#' @return a vector of pdf plot names. The plots can be found in Plots/TimeSeries.
#' 
#' @export 
#'
#' @examples



plotSimCheck <- function(SimData, mcmc.out, VitalRates = TRUE, DetectParams = TRUE, PopSizes = TRUE, Densities = TRUE) {
  
  require(coda)
  require(tidyverse)
  require(tidybayes)
  require(ggforce)
  require(see)
  require(cowplot)
  
  
  # Setup #
  #-------#
  
  ## Plotting directory
  ifelse(!dir.exists("Plots/SimCheck"), dir.create("Plots/SimCheck"), FALSE) ## Check if folder exists, if not create folder
  
  ## List of plot names (required for targets integration)
  plot.paths <- c()
  
  
  # Plotting vital rates #
  #----------------------#
  
  if(VitalRates) {
    
    ## Prepare tidy data
    R_year <- mcmc.out %>% spread_draws(R_year[year])
    Mu_R <- mcmc.out %>% spread_draws(Mu.R) %>% mutate(lab_code = "Mu.R")
    sigmaT_R <- mcmc.out %>% spread_draws(sigmaT.R) %>% mutate(lab_code = "sigmaT.R")
    
    Mu_S1 <- mcmc.out %>% spread_draws(Mu.S1) %>% mutate(Surv = "S1") %>% rename(S = Mu.S1) %>% select(S, Surv)
    Mu_S <- mcmc.out %>% spread_draws(Mu.S) %>% mutate(Surv = "S") %>% rename(S = Mu.S) %>% select(S, Surv)
    Mu_S_data <-  tibble(S = Mu_S$S/Mu_S1$S, Surv = "S2") %>% bind_rows(., Mu_S1, Mu_S)
    
    ## Plot average survival probabilities (Mu.S)
    p1 <- ggplot(data = Mu_S_data, aes(x = Surv, y = S)) +
      geom_violinhalf(fill = "lightgreen", alpha = 0.7)  +
      ylim(0, 1) +
      geom_point(aes(x = "S", y = SimData$SimParams$Mu.S), 
                 col = "darkblue", size = 3) +
      geom_point(aes(x = "S1", y = sqrt(SimData$SimParams$Mu.S)), 
                 col = "darkblue", size = 3) + 
      geom_point(aes(x = "S2", y = sqrt(SimData$SimParams$Mu.S)), 
                 col = "darkblue", size = 3) + 
      theme_cowplot() +
      theme(panel.grid.major.y = element_line(color = "#8ccde3",  size = 0.25, linetype = 2), 
            text = element_text(size = 10)) +
      ylab("Survival probability") +
      xlab("Survival component")
    
    
    ## Plot average recruitment rate (Mu.R)
    p2 <- ggplot(data = Mu_R, aes(x = lab_code, y = Mu.R)) +
      geom_violinhalf(fill = "lightgreen", alpha = 0.7)  +
      ylim(0, 4) +
      geom_point(aes(x = "Mu.R", y = SimData$SimParams$Mu.R), 
                 col = "darkblue", size = 3) +
      theme_cowplot() +
      theme(panel.grid.major.y = element_line(color = "#8ccde3", size = 0.25, linetype = 2), 
            text = element_text(size = 10)) +
      ylab("Recruitment rate") +
      xlab("")
    
    ## Plot temporal variation in recruitment rate (sigma.R)
    p3 <- ggplot(data = sigmaT_R, aes(x = lab_code, y = sigmaT.R)) +
      geom_violinhalf(fill = "lightgreen", alpha = 0.7)  +
      ylim(0, max(sigmaT_R$sigmaT.R)+1) +
      geom_point(aes(x = "sigmaT.R", y = SimData$SimParams$sigmaT.R), 
                 col = "darkblue", size = 3) +
      theme_cowplot() +
      theme(panel.grid.major.y = element_line(color = "#8ccde3", size = 0.25, linetype = 2), 
            text = element_text(size = 10)) +
      ylab("Year RE SD") +
      xlab("")
    p3
    
    ## Plot annual recruitment rates (R_year)
    Na_temp <- apply(SimData$N.data, c(2,3), sum) 
    R_simData <- as_tibble(t(Na_temp)) %>%
      dplyr::mutate(realizedR = V1 / V2, year = seq(1:SimData$SimParams$Tmax)) %>%
      dplyr::select(year, realizedR) %>%
      dplyr::mutate(predictedR = colMeans(AllSimData$VR.list$R))
    
    p4 <- ggplot(data = R_year, aes(x = as.factor(year), y = R_year)) +
      geom_violinhalf(fill = "lightgreen", alpha = 0.7)  +
      ylim(0, 5) +
      geom_point(data = R_simData, aes(x = as.factor(year), y = predictedR), 
                 col = "darkblue", size = 3) +
      geom_point(data = R_simData, aes(x = as.factor(year), y = realizedR), 
                 col = "purple", size = 3, alpha = 0.5) +
      theme_cowplot() +
      theme(panel.grid.major.y = element_line(color = "#8ccde3", size = 0.25, linetype = 2), 
            text = element_text(size = 10)) +
      ylab("R") +
      xlab("Year")
    
    ## Putting the plots together in a multipanel plot
    p_a <- plot_grid(p1, p2, p3, nrow = 1)
    p_out <- plot_grid(p_a, p4, nrow = 2, label_size = 10)
    
    ## Plot to pdf
    pdf("Plots/SimCheck/SimCheck_VRs.pdf", width = 12, height = 7.5)
    print(p_out)
    dev.off()
    
    plot.paths <- c(plot.paths, "SimCheck_VRs.pdf")
  }
  
  
  # Plotting detection parameters #
  #-------------------------------#
  
  if(DetectParams){
    
    ## Mean half normal scale parameter (mu.dd)
    Mu_dd <- mcmc.out %>% spread_draws(mu.dd) %>% mutate(lab_code = "mu.dd")
    
    p1 <- ggplot(data = Mu_dd, aes(x = lab_code, y = exp(mu.dd))) +
      geom_violinhalf(fill = "lightgreen", alpha = 0.7)  +
      #ylim(0, 4) +
      geom_point(aes(x = "mu.dd", y = SimData$SimParams$Mu.dd), 
                 col = "darkblue", size = 3) +
      theme_cowplot() +
      theme(panel.grid.major.y = element_line(color = "#8ccde3", size = 0.25, linetype = 2), 
            text = element_text(size = 10)) +
      ylab("HN scale parameter") +
      xlab("")
    
    ## Temporal variation in scale parameter (sigmaT.dd) - i.e. sd of random var.
    sigmaT_dd <- mcmc.out %>% spread_draws(sigmaT.dd) %>% mutate(lab_code = "sigmaT.dd")
    
    p2 <- ggplot(data = sigmaT_dd, aes(x = lab_code, y = sigmaT.dd)) +
      geom_violinhalf(fill = "lightgreen", alpha = 0.7)  +
      #ylim(0, 4) +
      geom_point(aes(x = "sigmaT.dd", y = SimData$SimParams$sigmaT.dd), 
                 col = "darkblue", size = 3) +
      theme_cowplot() +
      theme(panel.grid.major.y = element_line(color = "#8ccde3", size = 0.25, linetype = 2), 
            text = element_text(size = 10)) +
      ylab("sigmaT.dd") +
      xlab("")
    
    ## Effective strip width across years
    esw_year <- mcmc.out %>% spread_draws(esw[year])
    
    p4 <- ggplot(data = esw_year, aes(x = as.factor(year), y = esw)) +
      geom_violinhalf(fill = "lightgreen", alpha = 0.7)  +
      # ylim(0, 160) +
      theme_cowplot() +
      theme(panel.grid.major.y = element_line(color = "#8ccde3", size = 0.25, linetype = 2), 
            text = element_text(size = 10)) +
      ylab("Eff. strip width") +
      xlab("Year")
    #TODO: Add theoretical true values calculated from simulated data
    
    ## Detection probability across years
    p_year <- mcmc.out %>% spread_draws(p[year])
    
    p5 <- ggplot(data = p_year, aes(x = as.factor(year), y = p)) +
      geom_violinhalf(fill = "lightgreen", alpha = 0.7)  +
      # ylim(0, 160) +
      theme_cowplot() +
      theme(panel.grid.major.y = element_line(color = "#8ccde3", size = 0.25, linetype = 2), 
            text = element_text(size = 10)) +
      ylab("Detection probability") +
      xlab("Year")
    #TODO: Add theoretical true values calculated from simulated data

    ## Assemble plots
    p3 <- plot_grid(p1, p2, nrow = 1)
    p_out <- plot_grid(p3, p4, p5, nrow = 3)
    
    ## Plot to pdf
    pdf("Plots/SimCheck/SimCheck_Detects.pdf", width = 10, height = 7.5)
    print(p_out)
    dev.off()
    
    plot.paths <- c(plot.paths, "SimCheck_Detects.pdf")
    
  } 
  
  if(PopSizes){
    message("Plotting ov SimChecks for population sizes not yet supported.") 
    #TODO: Add code for population sizes
  }
  
  
  if(Densities){
    
    # Prepare density (for each time step) for simulated data
    N_temp <- apply(SimData$N.data, 3, sum) 
    A_temp <- apply(SimData$DS.data$L, 2, sum) * SimData$SimParams$W*2 / (1000 *1000)
    D_temp <- tibble(year=seq(1:SimData$SimParams$Tmax), Mean.D = N_temp / A_temp) 
    
    # Summarize mcmc data - annual density. This could (should) be done with 
    # spread_draws()
    
    Density_year <- mcmc.out %>% spread_draws(N_tot_exp[year]) %>% 
      dplyr::mutate(density = (N_tot_exp/A_temp)) 
    
    p_out <- ggplot(data = Density_year, aes(x = as.factor(year), y = density)) +
      geom_violinhalf(fill = "lightgreen", alpha = 0.7)  +
      geom_point(data=D_temp, aes(x = as.factor(year), y = Mean.D), 
                 col = "darkblue", size = 3) +
      theme_cowplot() +
      theme(panel.grid.major.y = element_line(color = "#8ccde3", size = 0.25, linetype = 2), 
            text = element_text(size = 10)) +
      ylab("Density") +
      xlab("Year")
    p_out
    
    ## Plot to pdf
    pdf("Plots/SimCheck/SimCheck_Density.pdf", width = 12, height = 5)
    print(p_out)
    dev.off()
    
    plot.paths <- c(plot.paths, "SimCheck_Density.pdf")
    
  }  
  
  return(plot.paths)
}