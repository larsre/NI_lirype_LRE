
# SETUP #
#-------#

## Set switches 

# (Re-)downloading data
# downloadData <- FALSE
downloadData <- TRUE


# DOWNLOAD/FETCH DATA #
#---------------------#

if(downloadData){
  Rype_arkiv <- downloadLN(version = 1.6, save = TRUE)
}else{
  stop("downloadData = FALSE not supported yet. There is an issue with encoding when using LivingNorwayR::initializeDwCArchive() that needs to be resolved first.")
  #Rype_arkiv <- initializeDwCArchive("data/Rype_arkiv.zip")
}


# WRANGLE LINE TRANSECT DATA #
#----------------------------#

## Set localities and time period of interest
localities <- "Lierne Fjellst. Vest"
minYear <- 2015
maxYear <- 2020

## Extract transect and observational data from DwC archive
LT_data <- wrangleData_LineTrans(DwC_archive = Rype_arkiv, 
                                 localities = localities,
                                 minYear = minYear, maxYear = maxYear)


# WRANGLE KNOWN FATE CMR DATA #
#-----------------------------#

## Read in and reformat CMR data
d_cmr <- wrangleData_CMR()


# PREPARE INPUT DATA FOR INTEGRATED MODEL #
#-----------------------------------------#

## Reformat data into vector/array list for analysis with Nimble
input_data <- prepareInputData(d_trans = LT_data$d_trans, 
                               d_obs = LT_data$d_obs,
                               d_cmr = d_cmr,
                               dataVSconstants = TRUE,
                               save = TRUE)

