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

# (Re-)download data (TRUE) or use existing data in 'data' folder (FALSE)
downloadData <- FALSE

# NI: estimate parameters by pre-2020 counties (FALSE; default NI) or by post-2020 counties (TRUE)
counties2020 <- FALSE

# Recruitment per adult pair (FALSE; default) or per adult female (TRUE)
R_perF <- FALSE

# Drop observations of juveniles with no adults present
R_parent_drop0 <- TRUE

# Aggregation level for reproduction data
# NOTE: if this is not defined, will default to group level
sumR.Level <- "line" # Summing at the line level

# Random effects shared across areas (default TRUE)
shareRE <- TRUE

# Estimate time variation in survival (default FALSE)
survVarT <- FALSE

# Include rodent covariate on reproduction (default FALSE)
fitRodentCov <- FALSE


# DOWNLOAD/FETCH DATA #
#---------------------#

if(downloadData) {
  # (Re-)download data from GBIF
  Rype_arkiv <- downloadLN(datasets = c("Fjellstyrene", "Statskog", "FeFo"), versions = c(1.7, 1.8, 1.12), save = FALSE)
} else {
  # Use existing data sets if available
  # NOTE: Warnings about existing temporary downloads in AppData have been suppressed
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

## Set localities/areas/counties and time period of interest
#localities <- listLocations() #not relevant for NI analysis
areas <- listAreas()

if(counties2020) {
  counties <- listCounties2020()
} else {
  counties <- listCounties()
}

minYear <- 2007
maxYear <- 2021

## List duplicate transects to remove
duplTransects <- listDuplTransects()

## Extract transect and observational data from DwC archive
# NOTE: for county aggregation, the data must first be wrangled using area
# aggregation/filtering
LT_data <- wrangleData_LineTrans(DwC_archive_list = Rype_arkiv, 
                                 duplTransects = duplTransects,
                                 #localities = localities,
                                 areas = areas,
                                 areaAggregation = TRUE,
                                 minYear = minYear,
                                 maxYear = maxYear)

## Assign transect lines and observations to counties by geographical coordinates
# NOTE: Set 'counties = NULL' to avoid filtering on listCounties()/listCounties2020()
LT_data <- assignCounty(LT_data,
                        counties = counties,
                        counties2020 = counties2020)


# WRANGLE KNOWN FATE CMR DATA #
#-----------------------------#

if(survVarT) {
  ## Read in and reformat CMR data
  d_cmr <- wrangleData_CMR(minYear = minYear)
}


# WRANGLE RODENT DATA #
#---------------------#

if(fitRodentCov) {
  ## Load and reformat rodent data
  d_rodent <- wrangleData_Rodent(duplTransects = duplTransects,
                                 #localities = localities,
                                 areas = areas,
                                 areaAggregation = TRUE,
                                 minYear = minYear, maxYear = maxYear)
}

# PREPARE INPUT DATA FOR INTEGRATED MODEL #
#-----------------------------------------#

## Reformat data into vector/array list for analysis with Nimble
# NOTE: set 'd_cmr' and/or 'd_rodent' to NULL if excluding survival and covariates
input_data <- prepareInputData(d_trans = LT_data$d_trans, 
                               d_obs = LT_data$d_obs,
                               d_cmr = NULL, #d_cmr
                               d_rodent = NULL, #d_rodent
                               #localities = localities, 
                               #areas = areas,
                               counties = counties,
                               areaAggregation = FALSE,
                               countyAggregation = TRUE,
                               excl_neverObs = TRUE,
                               R_perF = R_perF,
                               R_parent_drop0 = R_parent_drop0,
                               sumR.Level = sumR.Level,
                               dataVSconstants = TRUE,
                               save = TRUE)


# MODEL SETUP #
#-------------#

## Determine correct code path
code.path <- selectCodePath(shareRE = shareRE,
                            survVarT = survVarT)

## Setup for model using nimbleDistance::dHN
model_setup <- setupModel(modelCode.path = code.path,
                          customDist = TRUE,
                          nim.data = input_data$nim.data,
                          nim.constants = input_data$nim.constants,
                          R_perF = R_perF,
                          shareRE = shareRE,
                          survVarT = survVarT,
                          fitRodentCov = fitRodentCov,
                          testRun = TRUE,
                          nchains = 3,
                          initVals.seed = 0)

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

saveRDS(IDSM.out, file = 'rypeIDSM_dHN_multiArea_realData_testRun.rds')


# TIDY UP POSTERIOR SAMPLES #
#---------------------------#

IDSM.out.tidy <- tidySamples(IDSM.out = IDSM.out, save = TRUE)


# PREPARE OUTPUT DATA FOR NI #
#----------------------------#
output.data <- prepareOutputNI(mcmc.out = IDSM.out.tidy,
                               N_areas = input_data$nim.constant$N_areas, 
                               area_names = input_data$nim.constant$area_names, 
                               N_sites = input_data$nim.constant$N_sites, 
                               min_years = input_data$nim.constant$min_years, 
                               max_years = input_data$nim.constant$max_years, 
                               minYear = minYear,
                               maxYear = maxYear)
#list of data frames with (all values with lower/upper credible intervals)
#1) total density and recruitment
#2) total density and recruitment per year
#3) density and recruitment per area
#4) density and recruitment per area per year
#NB! values in database are in 1000 individuals per area (counties < 2020) - need to multiply density/km2 with available ptarmigan habitat
#
# Values extracted from the ptarmigan habitat model by Kvasnes et al. 2018:
# - Østfold (2010-2019 = 0); No data GBIF
# - Akershus (2010-2019 = 0); No data GBIF
# - Oslo (2010-2019 = 0); No data GBIF
# - Hedmark (2019 = 95)
# - Oppland (2019 = 67)
# - Buskerud (2019 = 57); Only 1 area (Øvre Numedal fjellstyre) in GBIF
# - Vestfold (2010-2019 = 0); No data GBIF
# - Telemark (2019 = 22); No data GBIF
# - Aust-Agder (2019 = 15); Only 1 area (Njardarheim) in GBIF
# - Vest-Agder (2019 = 18); No data GBIF
# - Rogaland (2019 = 4); No data GBIF
# - Hordaland (2019 = 19); Only 1 area (Eidfjord fjellstyre) in GBIF
# - Sogn og Fjordane (2019 = 10); No data GBIF
# - Møre og Romsdal (2019 = 14); No data GBIF
# - Sør-Trøndelag (2019 = 117)
# - Nord-Trøndelag (2019 = 129)
# - Nordland (2019 = 82)
# - Troms (2019 = 83)
# - Finnmark (2019 = 89)

# EXPORT OUTPUT DATA TO NI DATABASE? #
#------------------------------------#
#export.output(output.data)
# interact with the database via R?


# OPTIONAL: MCMC TRACE PLOTS #
#----------------------------#
# NOTE: Takes about 30+ minutes to generate PDF's

plotMCMCTraces(mcmc.out = IDSM.out.tidy,
               fitRodentCov = fitRodentCov)


# OPTIONAL: TIME SERIES PLOTS #
#-----------------------------#

# NOTE: VitalRates = TRUE will only plot reproduction if survVarT is set to FALSE
plotTimeSeries(mcmc.out = IDSM.out.tidy, 
               N_areas = input_data$nim.constant$N_areas, 
               area_names = input_data$nim.constant$area_names, 
               N_sites = input_data$nim.constant$N_sites, 
               min_years = input_data$nim.constant$min_years, 
               max_years = input_data$nim.constant$max_years, 
               minYear = minYear, maxYear = maxYear,
               VitalRates = TRUE, DetectParams = TRUE, Densities = TRUE, survVarT = survVarT)


# OPTIONAL: PLOT VITAL RATE POSTERIORS #
#--------------------------------------#

# plotPosteriorDens_VR(mcmc.out = IDSM.out.tidy,
#                      N_areas = input_data$nim.constant$N_areas, 
#                      area_names = input_data$nim.constant$area_names, 
#                      N_years = input_data$nim.constant$N_years,
#                      minYear = minYear,
#                      survAreaIdx = input_data$nim.constants$SurvAreaIdx,
#                      survVarT = survVarT,
#                      fitRodentCov = fitRodentCov) 


# OPTIONAL: PLOT COVARIATE PREDICTIONS #
#--------------------------------------#

# if(fitRodentCov){
#   plotCovPrediction(mcmc.out = IDSM.out.tidy,
#                     effectParam = "betaR.R",
#                     covName = "Rodent occupancy",
#                     minCov = 0, 
#                     maxCov = 1,
#                     N_areas = input_data$nim.constant$N_areas, 
#                     area_names = input_data$nim.constant$area_names,
#                     fitRodentCov = fitRodentCov)
# }


# OPTIONAL: MAP PLOTS #
#---------------------#

## Make map of Norwegian municipalities ("fylke")
# NorwayMunic.map <- setupMap_NorwayMunic(shp.path = "data/Kommuner_2018_WGS84/Kommuner_2018_WGS84.shp",
#                                         d_trans = LT_data$d_trans,
#                                         areas = areas, areaAggregation = TRUE)

## Plot population growth rate, density, and vital rates on map
# plotMaps(mcmc.out = IDSM.out.tidy, 
#          mapNM = NorwayMunic.map,
#          N_areas = input_data$nim.constant$N_areas, 
#          area_names = input_data$nim.constant$area_names, 
#          N_sites = input_data$nim.constant$N_sites, 
#          min_years = input_data$nim.constant$min_years, 
#          max_years = input_data$nim.constant$max_years, 
#          minYear = minYear, maxYear = maxYear,
#          fitRodentCov = fitRodentCov)


# OPTIONAL: MODEL COMPARISON (PLOTS) #
#------------------------------------#

# modelComp <- plotModelComparison(modelPaths = c("rypeIDSM_realData_Lierne.rds", 
#                                                 "rypeIDSM_dHN_realData_Lierne.rds"), 
#                                  modelChars = c("Zeroes trick", "dHN"), 
#                                  N_sites = 58, N_years = 6,
#                                  plotPath = "Plots/ModelCompTest",
#                                  returnData = FALSE)

