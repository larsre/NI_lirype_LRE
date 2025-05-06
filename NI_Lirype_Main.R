#install.packages("remotes","nimble","EnvStats","extraDistr","reshape2")
#remotes::install_github("LivingNorway/LivingNorwayR")
#remotes::install_github("scrogster/nimbleDistance")
#devtools::install_github("NINAnor/NIcalc", build_vignettes = T)
library(tidyverse)
library(nimble)
library(nimbleDistance)
library(LivingNorwayR)
library(reshape2)


# SETUP #
#-------#

## Set seed
mySeed <- 32
set.seed(mySeed)

## Set number of chains
nchains <- 3

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

# (Re-)download data (TRUE) or use existing data in 'data' folder (FALSE)
downloadData <- TRUE

# Aggregation to area level
areaAggregation <- TRUE

# NI: estimate parameters by pre-2020 counties (FALSE; default NI) or by post-2020 counties (TRUE)
counties2020 <- FALSE

# Recruitment per adult pair (FALSE; default) or per adult female (TRUE)
R_perF <- FALSE

# Drop observations of juveniles with no adults present
R_parent_drop0 <- TRUE

# Aggregation level for reproduction data
# NOTE: if this is not defined, will default to group level
sumR.Level <- "line" # Summing at the line level

# Estimate time variation in survival (default FALSE)
survVarT <- FALSE

# Include rodent covariate on reproduction (default FALSE)
fitRodentCov <- FALSE

# Use of telemetry data from Lierne (default FALSE)
telemetryData <- FALSE

# Test run or not
testRun <- TRUE


# DOWNLOAD/FETCH DATA #
#---------------------#

if(downloadData) {
  # (Re-)download data from GBIF
  Rype_arkiv <- downloadLN(datasets = c("Fjellstyrene", "Statskog", "FeFo"), versions = c(1.9, 1.9, 1.14), save = TRUE)
} else {
  # Use existing data sets if available
  # NOTE: Importing existing data sets doesn't work properly with the newest
  #       versions of the data due to JSON parsing issues with 'purrr'.
  #       Use switch 'downloadData <- TRUE' for now.
  Rype_arkiv <- list()
  data.exists <- 0
  
  if(file.exists("data/FeFo_arkiv.zip")) {
    Rype_arkiv[[1]] <- suppressWarnings(LivingNorwayR::initializeDwCArchive("data/FeFo_arkiv.zip"))
    data.exists <- data.exists + 1
  }
  if(file.exists("data/Statskog_arkiv.zip")) {
    Rype_arkiv[[2]] <- suppressWarnings(LivingNorwayR::initializeDwCArchive("data/Statskog_arkiv.zip"))
    data.exists <- data.exists + 1
  }
  if(file.exists("data/Fjellstyrene_arkiv.zip")) {
    Rype_arkiv[[3]] <- suppressWarnings(LivingNorwayR::initializeDwCArchive("data/Fjellstyrene_arkiv.zip"))
    data.exists <- data.exists + 1
  }
  
  if(identical(data.exists, as.numeric(0))) {
    stop("No data sets found in /data/ folder. Set 'downloadData' to TRUE to (re-)download the data sets from GBIF.")
  } else {
    cat(paste("  ", data.exists, " data archives found and imported.", "\n", sep = ""))
  }
}


# WRANGLE LINE TRANSECT DATA #
#----------------------------#

## Set areas/counties of interest
# NOTE: listAreas is required to extract the specific ptarmigan sampling areas.
#       The list is conformed to areas that are present in the required small
#       rodent table.
areas <- listAreas()

if(counties2020) {
  counties <- listCounties2020()
} else {
  counties <- listCounties()
}

## Set time period of interest
minYear <- 2007
maxYear <- 2024

## List duplicate transects to remove
duplTransects <- listDuplTransects()

## Extract transect and observational data from DwC archive
# NOTE: for county aggregation, the data must first be wrangled using area
# aggregation/filtering
LT_data_orig <- wrangleData_LineTrans(DwC_archive_list = Rype_arkiv, 
                                 duplTransects = duplTransects,
                                 #localities = localities,
                                 areas = areas,
                                 areaAggregation = areaAggregation,
                                 minYear = minYear,
                                 maxYear = maxYear)

## Add last year of data (not available from GBIF due to quarantine period)
# NOTE: this data is not (and should not be) publicly available, and is a direct
#       export from the Hønsefuglportalen (HFP) database. Ask the administrator
#       of HFP for an export of data if there are any missing in the download
#       from GBIF. The data export is in a different format, but conformed here.
#       Note that this requires 'data/countData2024.rds' to be present.
LT_data_add <- addLastCount(LT_data_orig)


## Assign transect lines and observations to counties by geographical coordinates
# NOTE: Set 'counties = NULL' to avoid filtering on listCounties()/listCounties2020()
LT_data <- assignCounty(LT_data_add,
                        counties = counties,
                        counties2020 = counties2020)


# WRANGLE KNOWN FATE CMR DATA #
#-----------------------------#

## Read in and reformat CMR data
# NOTE: Even with 'survVarT' and 'telemetryData' set to FALSE, CMR-data must be
#       included for the model to work properly. Note that this requires
#       'data/CMR_Data.csv' to be present.
d_cmr <- wrangleData_CMR(minYear = minYear)
d_cmr$county_names <- c("Nord-Trøndelag")


# WRANGLE RODENT DATA #
#---------------------#

## Load and reformat rodent data
# NOTE: Even with 'fitRodentCov' set to FALSE, the rodent data must be included
#       for all the other functions to work properly. Note that this requires
#       'data/Rodent_data.rds' to be present, and that areas included in the
#       main data set conforms to the areas present in the rodent data.
d_rodent <- wrangleData_Rodent(duplTransects = duplTransects,
                               areas = areas,
                               areaAggregation = areaAggregation,
                               minYear = minYear,
                               maxYear = maxYear)


# PREPARE INPUT DATA FOR INTEGRATED MODEL #
#-----------------------------------------#

## Reformat data into vector/array list for analysis with Nimble
input_data <- prepareInputData(d_trans = LT_data$d_trans, 
                               d_obs = LT_data$d_obs,
                               d_cmr = d_cmr,
                               d_rodent = d_rodent,
                               counties = counties,
                               areaAggregation = areaAggregation,
                               countyAggregation = TRUE,
                               excl_neverObs = TRUE,
                               R_perF = R_perF,
                               R_parent_drop0 = R_parent_drop0,
                               sumR.Level = sumR.Level,
                               dataVSconstants = TRUE,
                               save = FALSE)


# MODEL SETUP #
#-------------#

## Write model code
modelCode <- writeModelCode(survVarT = survVarT,
                            telemetryData = telemetryData)

## Expand seeds for simulating initial values
MCMC.seeds <- expandSeed_MCMC(seed = mySeed, 
                              nchains = nchains)

## Setup for model using nimbleDistance::dHN
model_setup <- setupModel(modelCode = modelCode,
                          R_perF = R_perF,
                          survVarT = survVarT,
                          fitRodentCov = fitRodentCov,
                          nim.data = input_data$nim.data,
                          nim.constants = input_data$nim.constants,
                          testRun = testRun,
                          nchains = nchains,
                          initVals.seed = MCMC.seeds)

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
                       setSeed = MCMC.seeds)
Sys.time() - t.start

saveRDS(IDSM.out, file = 'rypeIDSM_dHN_multiArea_sepRE_test_2025-04-02.rds')


# TIDY UP POSTERIOR SAMPLES #
#---------------------------#
IDSM.out.tidy <- tidySamples(IDSM.out = IDSM.out,
                             save = TRUE,
                             fileName = "rypeIDSM_dHN_multiArea_sepRE_test_2025-04-02_tidy.rds")

# load from previous run
#IDSM.out.tidy <- readRDS("rypeIDSM_dHN_multiArea_sepRE_test_2025-04-02_tidy.rds")

# PREPARE OUTPUT DATA FOR NI #
#----------------------------#
# Format population densities from posterior samples
prepared.output <- prepareOutputNI(
  mcmc.out = IDSM.out.tidy,
  N_areas = input_data$nim.constant$N_areas,
  area_names = input_data$nim.constant$area_names,
  N_sites = input_data$nim.constant$N_sites,
  min_years = input_data$nim.constant$min_years,
  max_years = input_data$nim.constant$max_years,
  minYear = minYear,
  maxYear = maxYear,
  fitRodentCov = fitRodentCov,
  save = FALSE
)

# Calculate abundance (in 1000 ind.) based on estimated densities and habitat information
# NOTE: This part depends on 'data/ptarmigan_habitat.rds', which is the output from the
#       ptarmigan habitat model by Kvasnes et al. 2018 BMC Ecology (with updated data).
#       If the rds file is absent, the function will utilize the provided data table
#       ('data/county_habitat_data.csv') with values from the last run of the model.
output.data <- estimatePtarmiganAbundance(output = prepared.output)


# WRANGLE OUTPUT DATA TO FIT NI DATABASE #
#----------------------------------------#
# retrieve the currently stored values for Ptarmigan from the Nature Index (NI) database
species <- c("Lagopus lagopus")
indicators <- c("Lirype")
currentPtarmTable <- downloadData_NIdb(species = species, indicators = indicators, save = T, save_path = "data")

# create a backup copy of the current indicator value table
#saveRDS(currentPtarmTable, "data/oldIndicatorValues_2025.rds")

# check county nomenclature consistency
all(unique(output.data$Area) %in% unique(currentPtarmTable$Lirype$indicatorValues$areaName)) #should be TRUE

# update the values with the results from the population density model for Ptarmigan.
# time interval to update can be specified using the 'min_year' and 'max_year' arguments
newPtarmTable <- updateNItable(model.est = output.data,
                               cur.table = currentPtarmTable,
                               save = T,
                               save_path = "data",
                               min_year = NULL,
                               max_year = NULL)


# UPLOAD NEW DATA NI DATABASE #
#-----------------------------#

# upload the updated table and overwrite existing data
uploadData_NIdb(indicatorData = newPtarmTable)


# OPTIONAL: MCMC TRACE PLOTS #
#----------------------------#

library(MCMCvis)
plotMCMCTraces(mcmc.out = IDSM.out.tidy,
               fitRodentCov = fitRodentCov,
               survVarT = survVarT)


# OPTIONAL: TIME SERIES PLOTS #
#-----------------------------#

# NOTE: VitalRates = TRUE will only plot reproduction if survVarT is set to FALSE
plotTimeSeries(mcmc.out = IDSM.out.tidy, 
               N_areas = input_data$nim.constant$N_areas, 
               area_names = input_data$nim.constant$area_names, 
               N_sites = input_data$nim.constant$N_sites, 
               min_years = input_data$nim.constant$min_years, 
               max_years = input_data$nim.constant$max_years, 
               minYear = minYear,
               maxYear = maxYear,
               VitalRates = TRUE,
               DetectParams = TRUE,
               Densities = TRUE,
               survVarT = survVarT)
