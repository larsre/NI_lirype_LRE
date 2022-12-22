plotMCMCTraces <- function(mcmc.out, VitalRates = TRUE, DetectParams = TRUE, PopSizes = TRUE, Densities = TRUE){
  
  ## Make parameter lists
  mVR_params <- c("mu.R", "h.mu.R", "h.sigma.R", "sigmaT.R",
                  "mu.D1", "sigma.D", "ratio.JA1",
                  "Mu.S1", "Mu.S2", "h.Mu.S1", "h.Mu.S2", "h.sigma.S1", "h.sigma.S2")
  
  tVR_params <- c("R_year")
  
  mDet_params <- c("mu.dd", "sigmaT.dd")
  
  tDet_params <- c("esw", "p", "sigma")
  
  N_params <- c("N_exp")
  
  Ntot_params <- c("N_tot_exp")
  
  Dens_params <- c("Density")
  
  
  ## Make plots and print to pdf
  
  # Vital rates
  if(VitalRates){
    pdf("Plots/MCMCTraces/MCMCtrace_mVR.pdf", width = 8, height = 6)
    MCMCvis::MCMCtrace(mcmc.out,
                       params = mVR_params,
                       ind = TRUE,
                       pdf = FALSE)
    dev.off()
    
    pdf("Plots/MCMCTraces/MCMCtrace_tVR.pdf", width = 8, height = 6)
    MCMCvis::MCMCtrace(mcmc.out,
                       params = tVR_params,
                       ind = TRUE,
                       pdf = FALSE)
    dev.off()
  }
  
  # Detection parameters
  if(DetectParams){
    pdf("Plots/MCMCTraces/MCMCtrace_mDet.pdf", width = 8, height = 6)
    MCMCvis::MCMCtrace(mcmc.out,
                       params = mDet_params,
                       ind = TRUE,
                       pdf = FALSE)
    dev.off()
    
    pdf("Plots/MCMCTraces/MCMCtrace_tDet.pdf", width = 8, height = 6)
    MCMCvis::MCMCtrace(mcmc.out,
                       params = tDet_params,
                       ind = TRUE,
                       pdf = FALSE)
    dev.off()
  }
  
  # Population sizes
  if(PopSizes){
    pdf("Plots/MCMCTraces/MCMCtrace_Ntot.pdf", width = 8, height = 6)
    MCMCvis::MCMCtrace(mcmc.out,
                       params = Ntot_params,
                       ind = TRUE,
                       pdf = FALSE)
    dev.off()
    
    pdf("Plots/MCMCTraces/MCMCtrace_N.pdf", width = 8, height = 6)
    MCMCvis::MCMCtrace(mcmc.out,
                       params = N_params,
                       ind = TRUE,
                       pdf = FALSE)
    dev.off()
  }
  
  # Population densities
  if(Densities){
    pdf("Plots/MCMCTraces/MCMCtrace_Dens.pdf", width = 8, height = 6)
    MCMCvis::MCMCtrace(mcmc.out,
                       params = Dens_params,
                       ind = TRUE,
                       pdf = FALSE)
    dev.off()
  }
  
}