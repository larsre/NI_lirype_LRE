
# GENERAL SETUP #
#---------------#

## Source all functions in "R" folder
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir('R')


## Set switches

# Re-simulate data
resimulate <- TRUE

# Recruitment per adult or per adult female
R_perF <- FALSE

# Drop observations of juveniles with no adults present
R_parent_drop0 <- TRUE



# SET SIMUALATION PARAMETERS #
#----------------------------#

# General simulation parameters
#---

mySeed <- 0

Amax <- 2 # Number of age classes
Tmax <- 15 # Number of years
Jmax <- 50 # Number of sites/transect lines


# Vital rate parameters
#---

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


# Population parameters
#---

# Initial population numbers per site
N1_juv_limits <- c(3, 8)

# Average group size
avg_Gsize <- 5.6


# Data & observation parameters 
#---

## Line-transect distance sampling
min.Tlength <- 1000 # Minimum transect length
max.Tlength <- 1000  # Maximum transect length

W <- 200 # Truncation distance (max. distance at which observation is possible)

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



# SIMULATE DATA #
#---------------#

if(resimulate){
  AllSimData <- assembleSimData(Amax = Amax, Tmax = Tmax, Jmax = Jmax,
                                avg_Gsize = avg_Gsize, 
                                Mu.S = Mu.S, sigmaT.S = sigmaT.S, sigmaJ.S = sigmaJ.S,
                                Mu.R = Mu.R, sigmaT.R = sigmaT.R, sigmaJ.R = sigmaJ.R,
                                Mu.dd = Mu.dd, sigmaT.dd = sigmaT.dd, sigmaJ.dd = sigmaJ.dd,
                                W = W, min.Tlength = min.Tlength, max.Tlength = max.Tlength,
                                nind.avg.RT = nind.avg.RT, 
                                Tmin.RT = Tmin.RT, Tmax.RT = Tmax.RT,
                                seed = mySeed, 
                                stochasticSim = TRUE,
                                plotPopSim = TRUE,
                                save = TRUE)
}else{
  AllSimData <- readRDS("SimData_Full.rds")
}


# PREPARE INPUT DATA FOR INTEGRATED MODEL #
#-----------------------------------------#

## Reformat data into vector/array list for analysis with Nimble
input_data <- list(
  nim.data = list(
    sumR_obs = AllSimData$Rep.data$sumR_obs,
    sumAd_obs = AllSimData$Rep.data$sumAd_obs,
    y = AllSimData$DS.data$d,
    L = AllSimData$DS.data$L,
    N_a_line_year = AllSimData$DS.data$DS.count,
    Survs1 = AllSimData$RT.data$Survs1,
    Survs2 = AllSimData$RT.data$Survs2
  ),
  
  nim.constants = list(
    N_years = AllSimData$SimParams$Tmax,
    year_Survs = AllSimData$RT.data$year_Survs,
    N_years_RT = AllSimData$RT.data$N_years_RT,
    W = AllSimData$SimParams$W,
    scale1 = scale1,
    A = colSums(AllSimData$DS.data$L)*(AllSimData$SimParams$W/scale1)*2,
    N_obs = length(AllSimData$DS.data$d),
    Year_obs = AllSimData$DS.data$d_year,
    N_sites = AllSimData$SimParams$Jmax,
    sumR_obs_year = AllSimData$Rep.data$sumR_obs_year,
    N_sumR_obs = AllSimData$Rep.data$N_sumR_obs,
    N_ageC = AllSimData$SimParams$Amax
  )
)


# MODEL SETUP #
#-------------#
  
model_setup <- setupModel(modelCode.path = "NIMBLE Code/RypeIDSM_dHN.R",
                          customDist = TRUE,
                          nim.data = input_data$nim.data,
                          nim.constants = input_data$nim.constants,
                          testRun = FALSE, initVals.seed = 0)

# MODEL (TEST) RUN #
#------------------#
t.start <- Sys.time()
IDSM.out <- nimbleMCMC(code = model_setup$modelCode, 
                       data = input_data$nim.data, 
                       constants = input_data$nim.constants,
                       inits = model_setup$initVals, 
                       monitors = model_setup$modelParams,
                       nchains = model_setup$mcmcParams$nchains, 
                       niter = model_setup$mcmcParams$niter, 
                       nburnin = model_setup$mcmcParams$nburn, 
                       thin = model_setup$mcmcParams$nthin, 
                       samplesAsCodaMCMC = TRUE, 
                       setSeed = 0)
Sys.time() - t.start

#saveRDS(IDSM.out, file = 'rypeIDSM_dHN_simData_t15.rds')
#saveRDS(IDSM.out, file = 'rypeIDSM_dHN_simData_t30.rds')
saveRDS(IDSM.out, file = 'rypeIDSM_dHN_simData_t15_MuR.rds')
#saveRDS(IDSM.out, file = 'rypeIDSM_dHN_simData_t30_MuR.rds')

