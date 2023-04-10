library(tidyverse)

# SETUP #
#-------#

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

# (Re-)downloading data
# downloadData <- FALSE
downloadData <- TRUE

# Recruitment per adult or per adult female
R_perF <- TRUE

# Drop observations of juveniles with no adults present
R_parent_drop0 <- FALSE

# DOWNLOAD/FETCH DATA #
#---------------------#

if(downloadData){
  #Rype_arkiv <- downloadLN(datasets = "Fjellstyrene", versions = 1.6, save = TRUE)
  Rype_arkiv <- downloadLN(datasets = c("Fjellstyrene", "Statskog", "FeFo"), versions = c(1.7, 1.8, 1.12), save = TRUE)
}else{
  stop("downloadData = FALSE not supported yet. There is an issue with encoding when using LivingNorwayR::initializeDwCArchive() that needs to be resolved first.")
  #Rype_arkiv <- initializeDwCArchive("data/Rype_arkiv.zip")
}


# WRANGLE LINE TRANSECT DATA #
#----------------------------#

## Set localities/areas and time period of interest
localities <- listLocations()
areas <- listAreas()
minYear <- 2007
maxYear <- 2021

## List duplicate transects to remove
duplTransects <- listDuplTransects()

## Extract transect and observational data from DwC archive
LT_data <- wrangleData_LineTrans(DwC_archive_list = Rype_arkiv, 
                                 duplTransects = duplTransects,
                                 #localities = localities,
                                 areas = areas,
                                 areaAggregation = TRUE,
                                 minYear = minYear, maxYear = maxYear)


# WRANGLE KNOWN FATE CMR DATA #
#-----------------------------#

## Read in and reformat CMR data
d_cmr <- wrangleData_CMR(minYear = minYear)


# PREPARE INPUT DATA FOR INTEGRATED MODEL #
#-----------------------------------------#

## Reformat data into vector/array list for analysis with Nimble
input_data <- prepareInputData(d_trans = LT_data$d_trans, 
                               d_obs = LT_data$d_obs,
                               d_cmr = d_cmr,
                               #localities = localities, 
                               areas = areas,
                               areaAggregation = TRUE,
                               excl_neverObs = TRUE,
                               R_perF = R_perF,
                               R_parent_drop0 = R_parent_drop0,
                               dataVSconstants = TRUE,
                               save = TRUE)


# MODEL SETUP #
#-------------#

# Original version (zeroes-trick)
# model_setup <- setupModel(modelCode.path = "NIMBLE Code/RypeIDSM.R",
#                           customDist = FALSE,
#                           nim.data = input_data$nim.data,
#                           nim.constants = input_data$nim.constants,
#                           testRun = FALSE, initVals.seed = 0)
  
# Updated version (nimbleDistance::dHN)
model_setup <- setupModel(modelCode.path = "NIMBLE Code/RypeIDSM_multiArea_dHN.R",
                          customDist = TRUE,
                          nim.data = input_data$nim.data,
                          nim.constants = input_data$nim.constants,
                          testRun = TRUE, nchains = 3,
                          initVals.seed = 0)

# MODEL (TEST) RUN #
#------------------#
#t.start <- Sys.time()
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
#Sys.time() - t.start

saveRDS(IDSM.out, file = 'rypeIDSM_dHN_multiArea_realData_Lierne.rds')


# OPTIONAL: MCMC TRACE PLOTS #
#----------------------------#

MCMCvis::MCMCtrace(IDSM.out,
                   params = c("esw", "p", "D",
                              "R_year", "mu.R", "h.mu.R", "h.sigma.R", "sigmaT.R",
                              "sigma", "mu.dd", "sigma.dd",
                              "Mu.D1", "sigma.D",
                              "Mu.S1", "Mu.S2", "h.Mu.S1", "h.Mu.S2", "h.sigma.S1", "h.sigma.S2",
                              "ratio.JA1"))


# OPTIONAL: MODEL COMPARISON (PLOTS) #
#------------------------------------#

# modelComp <- plotModelComparison(modelPaths = c("rypeIDSM_realData_Lierne.rds", 
#                                                 "rypeIDSM_dHN_realData_Lierne.rds"), 
#                                  modelChars = c("Zeroes trick", "dHN"), 
#                                  N_sites = 58, N_years = 6,
#                                  plotPath = "Plots/ModelCompTest",
#                                  returnData = FALSE)
