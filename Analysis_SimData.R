
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

# Random effects shared across areas
shareRE <- TRUE

# Time variation in survival
survVarT <- FALSE

# Rodent covariate on reproduction
fitRodentCov <- FALSE


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

## Add dummy dimensions to observational data
N_a_line_year <- addDummyDimension(AllSimData$DS.data$DS.count)
L <- addDummyDimension(AllSimData$DS.data$L)
y <- addDummyDimension(AllSimData$DS.data$d)
Year_obs <- addDummyDimension(AllSimData$DS.data$d_year)
sumR_obs <- addDummyDimension(AllSimData$Rep.data$sumR_obs)
sumAd_obs <- addDummyDimension(AllSimData$Rep.data$sumAd_obs)
sumR_obs_year <- addDummyDimension(AllSimData$Rep.data$sumR_obs_year)
N_sites <- c(AllSimData$SimParams$Jmax, NA)
N_obs <- c(length(AllSimData$DS.data$d), NA)
N_sumR_obs <- c(AllSimData$Rep.data$N_sumR_obs, NA)

## Reformat data into vector/array list for analysis with Nimble
input_data <- list(
  nim.data = list(
    sumR_obs = sumR_obs,
    sumAd_obs = sumAd_obs,
    y = y,
    L = L,
    N_a_line_year = N_a_line_year,
    Survs1 = AllSimData$RT.data$Survs1,
    Survs2 = AllSimData$RT.data$Survs2
  ),
  
  nim.constants = list(
    N_areas = 1,
    SurvAreaIdx = 1,
    N_years = AllSimData$SimParams$Tmax,
    year_Survs = AllSimData$RT.data$year_Survs,
    N_years_RT = AllSimData$RT.data$N_years_RT,
    W = AllSimData$SimParams$W,
    N_obs = N_obs,
    Year_obs = Year_obs,
    N_sites = N_sites,
    sumR_obs_year = sumR_obs_year,
    N_sumR_obs = N_sumR_obs,
    N_ageC = AllSimData$SimParams$Amax
  )
)


# MODEL SETUP #
#-------------#
  
model_setup <- setupModel(modelCode.path = "NIMBLE Code/RypeIDSM_multiArea_dHN.R",
                          customDist = TRUE,
                          shareRE = shareRE, survVarT = survVarT, fitRodentCov = fitRodentCov,
                          nim.data = input_data$nim.data,
                          nim.constants = input_data$nim.constants,
                          testRun = FALSE, nchains = 3,
                          initVals.seed = mySeed)

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

saveRDS(IDSM.out, file = paste0("rypeIDSM_dHN_simData_t", Tmax, ".rds"))



# MCMC TRACE PLOTS #
#------------------#

plotMCMCTraces(mcmc.out = IDSM.out,
               fitRodentCov = fitRodentCov)



