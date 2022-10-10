
#' Extract and filter transect and observation data from DwC archive
#'
#' @param DwC_archive DwCArchive. A Darwin Core archive containing line transect 
#' data for willow ptarmigan (Lagopus lagopus) in Norway.
#' @param localities string or vector of strings. Names of localities to extract
#'data for. 
#' @param minYear integer. Earliest year of data to extract.
#' @param maxYear integer. Latest year of data to extract.  
#'
#' @return list of 2 tibbles. `d_trans` contains information on transects 
#' (events). `d_obs`contains information on observations made along transects 
#' (distance to transect line, numbers of birds in each age/sex class observed,
#' etc.)
#' 
#' @export
#'
#' @examples

wrangleData_LineTrans <- function(DwC_archive, localities, minYear, maxYear){
  
  ## Ensure necessary packages are available
  require('tidyverse')
  require('LivingNorwayR')
  
  
  ## Extract relevant parts from DwC_archive
  
  # Core table
  Core <- DwC_archive$getCoreTable()
  
  # Complete event table
  Eve_all <- tibble::as_tibble(Core$exportAsDataFrame()) %>%
    dplyr::mutate(sampleSizeValue = as.numeric(sampleSizeValue))
  
  # Complete occurrence table
  Occ <- tibble::as_tibble(DwC_archive$getExtensionTables()[[1]]$exportAsDataFrame())
  
  
  ## Filter event data by locality and year
  Eve <- Eve_all %>% 
    dplyr::mutate(eventDate = as.Date(eventDate)) %>%
    dplyr::mutate(Year = lubridate::year(eventDate)) %>%
    dplyr::filter(locality %in% localities) %>%
    dplyr::filter(between(Year, minYear, maxYear))
  
  
  ## Assemble transect level info
  d_trans <- Eve %>% 
    dplyr::select(locationID, eventDate, eventID, modified, 
                  samplingProtocol, eventRemarks, sampleSizeValue, 
                  stateProvince, municipality, locality, 
                  verbatimLocality, locationRemarks, Year) %>%
    dplyr::filter(eventRemarks == "Line transect") %>%
    dplyr::mutate(locationRemarks = gsub("In the original data this is known as lineID ", '', locationRemarks)) %>%
    tidyr::separate(., col = locationRemarks, sep = ",", into = c("LineName", "locationRemarks")) %>%
    dplyr::select(-locationRemarks) 
  
  
  ## Assemble observation level info
  
  # Observations: distance to transect lines
  d_obsTemp <- Eve %>% 
    dplyr::select(locationID, locality, parentEventID, eventID, eventRemarks, 
                  dynamicProperties, eventDate) %>%
    dplyr::filter(eventRemarks == "Human observation") %>%
    dplyr::mutate(dynamicProperties = purrr::map(dynamicProperties, ~ jsonlite::fromJSON(.) %>% as.data.frame())) %>%
    tidyr::unnest(dynamicProperties) %>% 
    dplyr::rename(DistanceToTransectLine = "perpendicular.distance.in.meters.from.transect.line.as.reported.by.the.field.worker") %>%
    dplyr::mutate(DistanceToTransectLine = as.numeric(DistanceToTransectLine), 
                  Year = lubridate::year(eventDate)) %>%
    dplyr::select(locationID, locality, parentEventID, eventID, DistanceToTransectLine, Year)
  
  # Observations: remaining information (willow ptarmigan only)
  d_obs <- Occ %>% 
    dplyr::select(eventID, scientificName, individualCount, sex, lifeStage) %>%
    #dplyr::mutate(lifeStage=if_else(is.na(lifeStage), "unknown", lifeStage)) %>%
    dplyr::mutate(SexStage = str_c(sex, lifeStage, sep = "")) %>%
    dplyr::select(-sex, -lifeStage) %>%
    tidyr::spread(key = "SexStage", value = "individualCount", fill = 0) %>%
    dplyr::right_join(., d_obsTemp, by = "eventID") %>%
    dplyr::filter(scientificName == "Lagopus lagopus")
  
  
  ## Collate and return data
  LT_data <- list(d_trans = d_trans, d_obs = d_obs)
}
