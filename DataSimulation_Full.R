

assembleSimData <- function(Amax, Tmax, Jmax,
                            avg_Gsize, 
                            Mu.S, sigmaT.S, sigmaJ.S = 0,
                            Mu.R, sigmaT.R, sigmaJ.R = 0,
                            Mu.dd, sigmaT.dd, sigmaJ.dd,
                            W, min.Tlength, max.Tlength,
                            nind.avg.RT, 
                            Tmin.RT, Tmax.RT,
                            seed = NA, 
                            stochasticSim = TRUE,
                            plotPopSim = FALSE,
                            save = TRUE){
  
  #######################
  # PROCESS SIMULATIONS #
  #######################
  
  ## Set seed
  if(!is.na(seed)){
    set.seed(seed)
  }

  
  ## Simulate vital rates
  VR.list <- simulateVRs(Tmax = Tmax, Jmax = Jmax,
                         Mu.S = Mu.S, Mu.R = Mu.R, 
                         sigmaT.S = sigmaT.S, sigmaT.R = sigmaT.R,
                         sigmaJ.S = sigmaJ.S, sigmaJ.R = sigmaJ.R)
  
  
  ## Simulate population trajectories for all sites
  SimData <- simulatePopDyn(Amax = Amax, Tmax = Tmax, Jmax = Jmax, 
                            VR.list = VR.list, 
                            stochastic = stochasticSim, plot = plotPopSim)
  
  
  ## Simulate group aggregation in population
  G.age <- simulateGroups(Jmax = Jmax, Tmax = Tmax, 
                          N.age = SimData$N, 
                          avg_Gsize = avg_Gsize, discard0 = TRUE)
  
  
  #################################
  # OBSERVATIONAL DATA SIMULATION #
  #################################
  
  ## Simulate hierarchical distance sampling data
  DS.data <- simulateData_HDS(Jmax, Tmax, 
                              G.age = G.age, 
                              Mu.dd, sigmaT.dd, sigmaJ.dd,
                              W, discard0 = TRUE, 
                              min.Tlength = min.Tlength, max.Tlength = max.Tlength)
  
  
  ## Extract reproductive data
  Rep.data <- extractSimData_Rep(Jmax = Jmax, Tmax = Tmax, 
                                 DS.count = DS.data$DS.count)
  
  
  ## Simulate known-fate telemetry data
  RT.data <- simulateData_RT(nind.avg.RT, 
                             Tmin.RT, Tmax.RT, Tmax, 
                             SurvProbs = VR.list$S[1,])
  
  
  ### Simulate nest survey data
  #NS.data <- simulateData_Nests(nind.avg.NS, 
  #                              Tmin.NS, Tmax.NS, 
  #                              RepRates = VR.list$R[1,])
  
  
  #################
  # DATA ASSEMBLY #
  #################
  
  ## Collect all simulated parameter values into a list
  SimParams <- list(
    Seed = seed,
    Amax = Amax, Tmax = Tmax, Jmax = Jmax,
    avg_Gsize = avg_Gsize,
    Mu.S = Mu.S, Mu.R = Mu.R, 
    sigmaT.S = sigmaT.S, sigmaT.R = sigmaT.R, 
    sigmaJ.S = sigmaJ.S, sigmaJ.R = sigmaJ.R, 
    
    min.Tlength = min.Tlength,
    max.Tlength = max.Tlength,
    W = W,
    Mu.dd = Mu.dd, sigmaT.dd = sigmaT.dd, sigmaJ.dd = sigmaJ.dd,
    
    Tmin.RT = Tmin.RT, Tmax.RT = Tmax.RT, nind.avg.RT = nind.avg.RT
    #Tmin.NS = Tmin.NS, Tmax.NS = Tmax.NS, nind.avg.NS = nind.avg.NS
  )
  
  ## Collate and save all data
  AllSimData <- list(
    SimParams = SimParams,
    VR.list = VR.list,
    N.data = SimData$N,
    group.data = G.age,
    DS.data = DS.data,
    Rep.data = Rep.data,
    RT.data = RT.data
  )
  
  if(save){
    saveRDS(AllSimData, "SimData_Full.rds")
  }
  
  return(AllSimData)
}



