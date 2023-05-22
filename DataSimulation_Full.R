
mySeed <- 0
set.seed(mySeed)


#########################
# SIMULATION PARAMETERS #
#########################


# General simulation parameters #
#-------------------------------#

Amax <- 2 # Number of age classes
Tmax <- 15 # Number of years
#Tmax <- 30 # Number of years
Jmax <- 50 # Number of sites/transect lines


# Vital rate parameters #
#-----------------------#

## Annual survival
Mu.S <- 0.35 # Average annual survival probability
sigmaT.S <- 0 # SD of random year variation in survival
sigmaJ.S <- 0 # SD of random site variation in survival

## Reproduction
Mu.R <- 2 # Average number of chicks in August
sigmaT.R <- 0.5 # SD of random year variation in number of chicks
sigmaJ.R <- 0 # SD of random site variation in number of chicks

## Juvenile summer survival
#Mu.sJ <- 0.2 # Average summer survival of chicks
#sigmaT.sJ <- 0 # SD of random year variation in chick survival
#sigmaJ.sJ <- 0 # SD of random site variation in survival


# Population parameters #
#-----------------------#

# Initial population numbers per site
N1_juv_limits <- c(3, 8)

# Average group size
avg_Gsize <- 5.6


# Data & observation parameters #
#-------------------------------#

## Line-transect distance sampling
min.Tlength <- 1000 # Minimum transect length
max.Tlength <- 1000  # Maximum transect length

W <- 200 # Truncation distance (max. distance at which observation is possible)
  
#pi <- 3.141593 # Approximation of pi

Mu.dd <- 75 # Average width parameter for half-normal detection function
sigmaT.dd <- 0 # SD of random year variation in detection probability
sigmaJ.dd <- 0 # SD of random line variation in detection probability

## Known-fate radio-telemetry
Tmin.RT <- 5 # First year for which radio-telemetry data has been collected
Tmax.RT <- 10 # Last year for which radio-telemetry data has been collected

# Average number of individuals fitted with transmitters each year
nind.avg.RT <- 30


# ## Nest survey
# Tmin.NS <- 1 # First year for which nest survey data has been collected
# Tmax.NS <- 15 # Last year for which nest survey data has been collected
# 
# # Average number of nests monitored each year
# nind.avg.NS <- 40


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
