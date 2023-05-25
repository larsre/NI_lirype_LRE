

#' Simulate complete dataset for testing the integrated distance-sampling model
#'
#' @param Amax integer. Number of age classes to consider.
#' @param Tmax integer. Total number of years to simulate for.
#' @param Jmax integer. Total number of sites/transects to simulate for.
#' @param avg_Gsize numeric. Average group size. 
#' @param Mu.S numeric. Average annual survival probability.
#' @param sigmaT.S numeric. Standard deviation of random year effects on survival.
#' @param sigmaJ.S numeric. Standard deviation of random site effects on survival. Default = 0.
#' @param Mu.R numeric. Average number of recruits/adult. 
#' @param sigmaT.R numeric. Standard deviation of random year effects on reproduction. 
#' @param sigmaJ.R numeric. Standard deviation of random site effects on reproduction. Default = 0.
#' @param Mu.dd numeric. Average detection decline rate.
#' @param sigmaT.dd numeric. Standard deviation of random year effects on detection.
#' @param sigmaJ.dd numeric. Standard deviation of random site effects on detection. 
#' @param W numeric. Truncation distance.
#' @param min.Tlength numeric. Minimum transect length.
#' @param max.Tlength numeric. Maximum transect length.
#' @param nind.avg.RT integer. Average number of individuals fitted with transmitters
#' every season.
#' @param Tmin.RT integer. Index of the first year for which to simulate telemetry
#' data.  
#' @param Tmax.RT integer. Index of the last year for which to simulate telemetry 
#' data.
#' @param R_perF logical. If TRUE, treats recruitment rate as juvenile per adult female.
#' If FALSE, treats recruitment rate as juvenile per adult (sum of both sexes).
#' @param R_parent_drop0 logical. If TRUE, removes observations of juveniles without adults
#' from recruitment data. If FALSE, sets 1 as the number of adults/adults females when none
#' are observed. 
#' @param seed integer. Seed for data simulation. 
#' @param stochasticSim logical. If TRUE (default), population dynamics are
#' simulated including demographic stochasticitiy. If FALSE, deterministic 
#' dynamics are simulated instead. 
#' @param plotPopSim logical. If TRUE, plots simulated site-specific population size
#' over time. Default = FALSE. 
#' @param save logical. If TRUE (default), saves dataset as "SimData_Full.rds". 
#'
#' @return a list containing all data necessary for running the model, as well
#' as all simulation parameters and true quantitites. 
#' @export
#'
#' @examples

assembleSimData <- function(Amax, Tmax, Jmax,
                            avg_Gsize, 
                            Mu.S, sigmaT.S, sigmaJ.S = 0,
                            Mu.R, sigmaT.R, sigmaJ.R = 0,
                            Mu.dd, sigmaT.dd, sigmaJ.dd,
                            W, min.Tlength, max.Tlength,
                            nind.avg.RT, 
                            Tmin.RT, Tmax.RT,
                            R_perF, R_parent_drop0,
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
                            R_perF = R_perF,
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
                                 DS.count = DS.data$DS.count,
                                 R_perF = R_perF,
                                 R_parent_drop0 = R_parent_drop0)
  
  
  ## Simulate known-fate telemetry data
  RT.data <- simulateData_RT(nind.avg.RT, 
                             Tmin.RT, Tmax.RT, 
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



