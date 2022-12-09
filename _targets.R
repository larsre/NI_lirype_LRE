## Load packages required to define the pipeline
library(targets)
library(nimble)
library(nimbleDistance)

## Source all functions in "R" folder
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir('R')

## Set localities and time period of interest
#localities <- unname(read.csv("localities.csv", encoding = "utf-8"))[,1]
minYear <- 2007
maxYear <- 2020


## Set target-specific options such as packages.
tar_option_set(packages = c("LivingNorwayR", "tidyverse"))

## Define Targets List
list(
  tar_target(
    Rype_arkiv,
    downloadLN(datasets = c("Fjellstyrene", "Statskog", "FeFo"), versions = c(1.6, 1.7, 1.11), save = TRUE)
  ),
  
  tar_target(
    localities,
    listLocations()
  ),
  
  tar_target(
    areas,
    listAreas()
  ),
  
  tar_target(
    LT_data,
    wrangleData_LineTrans(DwC_archive_list = Rype_arkiv, 
                          #localities = localities,
                          areas = areas,
                          areaAggregation = TRUE,
                          minYear = minYear, maxYear = maxYear)
  ),
  
  tar_target(
    d_cmr,
    wrangleData_CMR(minYear = minYear)
  ),
  
  tar_target(
    input_data,
    prepareInputData(d_trans = LT_data$d_trans, 
                     d_obs = LT_data$d_obs,
                     d_cmr = d_cmr,
                     #localities = localities,
                     areas = areas,
                     areaAggregation = TRUE,
                     dataVSconstants = TRUE,
                     save = TRUE)
  ),
  
  tar_target(
    model_setup,
    setupModel(modelCode.path = "NIMBLE Code/RypeIDSM_multiArea_dHN.R",
               customDist = TRUE,
               nim.data = input_data$nim.data,
               nim.constants = input_data$nim.constants,
               testRun = TRUE, nchains = 3,
               initVals.seed = 0)
  ),
  
  tar_target(
    IDSM.out,
    nimbleMCMC(code = model_setup$modelCode,
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
  )
)

# Test by running tar_manifest(fields = all_of("command")) and tar_visnetwork() in the console

# Run workflow using tar_make() in the console
# Check with tar_network() in the console

# We can then use tar_read() and tar_load() to inspect and work with results
