

#########################
# VITAL RATE SIMULATION #
#########################

## Simulate vital rates
VR.list <- simulateVRs(Tmax, Jmax,
                       Mu.S, Mu.R, 
                       sigmaT.S, sigmaT.R,
                       sigmaJ.S, sigmaJ.R)


#########################
# POPULATION SIMULATION #
#########################

## Simulate population trajectories for all sites
SimData <- simulatePopDyn(Amax, Tmax, Jmax, VR.list, N1, stochastic = TRUE, plot = TRUE)


################################
# GROUP AGGREGATION SIMULATION #
################################

## Simulate group aggregation of population
G.age <- simulateGroups(Jmax, Tmax, N.age = SimData$N, avg_Gsize, discard0 = TRUE)

#################################
# LINE TRANSECT DATA SIMULATION #
#################################

## Simulate hierarchical distance sampling data
DS.data <- simulateData_HDS(Jmax, Tmax, G.age = G.age, 
                            Mu.dd, sigmaT.dd, sigmaJ.dd,
                            W, discard0 = TRUE, 
                            min.Tlength, max.Tlength)


## Extract reproductive data
Rep.data <- extractSimData_Rep(Jmax, Tmax, DS.count = DS.data$DS.count)


########################################
# KNOWN-FATE TELEMETRY DATA SIMULATION #
########################################

## Simulate known-fate telemetry data
RT.data <- simulateData_RT(nind.avg.RT, Tmin.RT, Tmax.RT, Tmax, SurvProbs = VR.list$S[1,])


###############################
# NEST SURVEY DATA SIMULATION #
###############################

# ## Simulate nest survey data
# NS.data <- simulateData_Nests(nind.avg.NS, Tmin.NS, Tmax.NS, RepRates = VR.list$R[1,])


#########################
# FULL DATASET ASSEMBLY #
#########################

## Collect all simulated parameter values into a list
SimParams <- list(
  Seed = mySeed,
  Amax = Amax, Tmax = Tmax, Jmax = Jmax,
  avg_Gsize = avg_Gsize,
  N1 = N1,
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

saveRDS(AllSimData, "SimData_Full.rds")
