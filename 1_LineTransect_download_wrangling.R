
library(RJSONIO)
library(xml2)
library(tidyverse)
library(here)
library(lubridate)

# Set switch for (re-)downloading data
downloadData <- FALSE
#downloadData <- TRUE


##########################################################
#### Getting line transect data from 
#### GBIF, and preparing for distance sampling model


#########################################################

if(downloadData){
  
  # First - locate the data set ID from GBIF
  # Could be retreived also through data set search via
  # rgif 
  datasetID <- "b49a2978-0e30-4748-a99f-9301d17ae119"
  dataset <- RJSONIO::fromJSON(paste0("http://api.gbif.org/v1/dataset/",datasetID,"/endpoint"))
  # The data set is available as a Darwin COre Archive zip file
  # from the "endpoint"; 
  endpoint_url <- dataset[[1]]$url 
  download.file(endpoint_url, destfile="data/temp.zip", mode="wb")
  unzip ("data/temp.zip", exdir = "data")
  
}

############################################################
## Data is now located in ../data - folder. 
## The procedure could be started from here once data are in ../data forlder
############################################################
############################################################
############################################################

## Preparing the line transct data for further analyses; 

# Reading data; 
Eve2 <- as_tibble(read.csv("data/event.txt", sep="\t", stringsAsFactors = TRUE, encoding = "UTF-8")) %>% 
     select(-id)
Occ <- as_tibble(read.csv("data/occurrence.txt", sep="\t", stringsAsFactors = TRUE, encoding = "UTF-8")) %>% 
       select(-id) 

################################################################################
# Filtering the data - lines from the study area; 
# We use data from Lierne Fjellstyre - West. 

Eve <- Eve2 %>% 
  mutate(eventDate=as.Date(eventDate)) %>%
  mutate(Year=year(eventDate)) %>%
  filter(locality=="Lierne Fjellst. Vest") %>%
  filter(between(Year, 2015, 2020))


## Transect level info; 
d_trans <- Eve %>% select(locationID, eventDate, eventID, modified, samplingProtocol, 
                          eventRemarks, sampleSizeValue, stateProvince, municipality,
                          locality, verbatimLocality, locationRemarks, Year) %>%
  filter(eventRemarks=="Line transect") %>%
  mutate(locationRemarks=gsub("In the original data this is known as lineID ", '', locationRemarks)) %>%
  separate(., col=locationRemarks, sep=",", into=c("LineName", "locationRemarks")) %>%
  select(-locationRemarks) 

## Observation level info
## First - extracting info about distance to transect line from event table; 
## Observations; 
d_obsTemp <- Eve %>% select(locationID, parentEventID, eventID, eventRemarks, dynamicProperties, eventDate) %>%
  filter(eventRemarks=="Human observation") %>%
  mutate(dynamicProperties=gsub(" }", '', dynamicProperties)) %>%
  separate(., col=dynamicProperties, sep=":", into=c("temp1", "DistanceToTransectLine")) %>%
  mutate(DistanceToTransectLine=as.numeric(DistanceToTransectLine), Year=year(eventDate)) %>%
  select(locationID, parentEventID, eventID, DistanceToTransectLine, Year)

### Then - adding obervations from occurrence table - using only willow ptarmigan observations; 

d_obs <- Occ %>% select(eventID, scientificName, individualCount, sex,
                         lifeStage) %>%
  #mutate(lifeStage=if_else(is.na(lifeStage), "unknown", lifeStage)) %>%
  mutate(SexStage=str_c(sex, lifeStage, sep = "")) %>%
  select(-sex, -lifeStage) %>%
  spread(key="SexStage", value= "individualCount", fill=0) %>%
  right_join(., d_obsTemp) %>%
  filter(scientificName=="Lagopus lagopus")
  
################################################################################
################################################################################

###################################################################################
## Preparing data for Distance Sampling Analysis using the Open population demographic open distance sampling model; 

## N_sites; Number of unique survey lines
## N_year; Number of unique years

N_sites <- n_distinct(d_trans$locationID)
N_years <- n_distinct(d_trans$Year)

#############################################
## Observation-data; 
# W <- truncation dist

W <- 200

temp_dist <- d_obs %>% filter(between(DistanceToTransectLine, -0.1, W)) %>% 
  select(Year, DistanceToTransectLine) %>%
  mutate(Year2=Year-(min(Year))+1)

## Distance to transect line and year
y <- temp_dist$DistanceToTransectLine
Year_obs <- temp_dist$Year2

## Number of obs
N_obs <- length(y)
zeros_dist <- rep(0, length(y))

#############################################
#############################################
### Recruitment; 
### Might want to update this; 

temp_Rec <- d_obs %>% filter(between(DistanceToTransectLine, -0.1, W)) %>%
  
  mutate(R=unknownJuvenile+unknownunknown) %>%
  mutate(Maletemp=unknownJuvenile+unknownunknown+FemaleAdult, MaleIndeks=if_else(Maletemp==0, 1,0)) %>%
  select(Year, R, MaleIndeks) %>%
  mutate(Year2=Year-(min(Year))+1) %>%
  filter(MaleIndeks==0) # --> drop all observations of only males

R_obs <- temp_Rec$R
R_obs_year <- temp_Rec$Year2
N_R_obs <- length(R_obs)


##############################################
##############################################
### Matrix with Transect length - 
### Site on rows, years on columns; 

TransLen <- d_trans %>% select(locationID, Year, sampleSizeValue) %>%
  mutate(sampleSizeValue=sampleSizeValue/1000) %>%
  reshape2::dcast(locationID~Year, value.var="sampleSizeValue", sum) %>%
  arrange(locationID)

L <- TransLen %>% select(-locationID) %>% as.matrix() 
colnames(L) <- NULL

## Total covered area; 
A <- TransLen %>% select(-locationID) 
colnames(A) <- NULL
A <- colSums(A)*(W/1000)*2


##############################################
##############################################
### preparing matrix with number of birds/line

temp <- TransLen %>% select(locationID)

TaksObs <- d_obs %>% filter(between(DistanceToTransectLine, -0.1, W)) %>%
  mutate(cs=unknownJuvenile+unknownunknown+FemaleAdult+MaleAdult) %>%
  reshape2::dcast(locationID~Year, value.var="cs", sum) %>%
  right_join(., temp, by=c("locationID"="locationID")) %>%
  #as_tibble() %>%
  replace(., is.na(.), 0) %>%
  arrange(locationID)


N_line_year <- TaksObs %>% select(-locationID) %>% as.matrix() 
colnames(N_line_year) <- NULL


### The objects created here will form the different list-elements in the 
### jags-data under "6_RunningModels.R"


##############################################
##############################################
### Saving data for analysis

rype.data <- list(
  R_obs = R_obs, # Observed numbers of recruits
  R_obs_year = R_obs_year, # Year of observed numbers of recruits
  N_R_obs = N_R_obs, # Total number of observations of numbers of recruits
  
  y = y, # Distance to transect line for each individual observation
  zeros_dist = zeros_dist, # Vector of 0's of same length as y
  Year_obs = Year_obs, # Year of each observation
  N_obs = N_obs, # Total number of obsercations
  
  N_line_year = N_line_year, # Number of birds observed per site per year
  L = L, # Transect length per site and year
 
  N_years = N_years, # Number of years with data
  N_sites = N_sites, # Total number of monitored sites
  
  A = A, # Total covered area per year
  W = W # Truncation distance
)

saveRDS(rype.data, file = "RypeData_forIM.rds")
