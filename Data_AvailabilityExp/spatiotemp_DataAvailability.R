library(tidyr)

#------------------------------------------------#
# DOWNLOAD AND FORMAT DATASETS FROM LIVINGNORWAY #
#------------------------------------------------#

## Set dataset keys
key_Fjellstyrene <- "b49a2978-0e30-4748-a99f-9301d17ae119"
key_Statskog <- "6a948a1c-7e23-4d99-b1c1-ec578d0d3159"
key_FeFo <- "c47f13c1-7427-45a0-9f12-237aad351040"

## Download DwC archives
DwC_Fjellstyrene <- LivingNorwayR::getLNportalData(datasetKey = key_Fjellstyrene)
DwC_Statskog <- LivingNorwayR::getLNportalData(datasetKey = key_Statskog)
DwC_FeFo <- LivingNorwayR::getLNportalData(datasetKey = key_FeFo)

## Extract event and occurrence tables
event_Fjellstyrene <- tibble::as_tibble(DwC_Fjellstyrene$getCoreTable()$exportAsDataFrame()) %>%
  dplyr::mutate(sampleSizeValue = as.numeric(sampleSizeValue),
                eventDate = as.Date(eventDate))
event_Statskog <- tibble::as_tibble(DwC_Statskog$getCoreTable()$exportAsDataFrame()) %>%
  dplyr::mutate(sampleSizeValue = as.numeric(sampleSizeValue),
                eventDate = as.Date(eventDate))
event_FeFo <- tibble::as_tibble(DwC_FeFo$getCoreTable()$exportAsDataFrame()) %>%
  dplyr::mutate(sampleSizeValue = as.numeric(sampleSizeValue),
                eventDate = as.Date(eventDate))

occ_Fjellstyrene <- tibble::as_tibble(DwC_Fjellstyrene$getExtensionTables()[[2]]$exportAsDataFrame())
occ_Statskog <- tibble::as_tibble(DwC_Statskog$getExtensionTables()[[2]]$exportAsDataFrame())
occ_FeFo <- tibble::as_tibble(DwC_FeFo$getExtensionTables()[[2]]$exportAsDataFrame())

## Combine event and occurrence tables
event_Fjellstyrene$dataSet <- "Fjellstyrene"
event_Statskog$dataSet <- "Statskog"
event_FeFo$dataSet <- "FeFo"

occ_Fjellstyrene$dataSet <- "Fjellstyrene"
occ_Statskog$dataSet <- "Statskog"
occ_FeFo$dataSet <- "FeFo"

event <- dplyr::bind_rows(event_Fjellstyrene, event_Statskog, event_FeFo) 
occ <- dplyr::bind_rows(occ_Fjellstyrene, occ_Statskog, occ_FeFo)


#--------------------------------------------#
# CHECK FOR POTENTIAL ERRORS/INCONSISTENCIES #
#--------------------------------------------#

## Add year and extract line names
event <- event %>%
  dplyr::mutate(
    year = lubridate::year(eventDate),
    lineName = paste0(verbatimLocality, ' - ', stringi::stri_extract_first_regex(locationRemarks, "[0-9]+"))
  )

## Check for suspicious transect lengths
hist(event$sampleSizeValue, breaks = 200)
# --> Normal distribution with a heavy right tail
# --> But there seems to be some 0-inflation!

hist(subset(event, sampleSizeValue < 400)$sampleSizeValue, breaks = 100)
hist(subset(event, sampleSizeValue < 1000)$sampleSizeValue, breaks = 100)
# --> Transects with reported lengths < 200 may need a closer look

transect_issues <- subset(event, sampleSizeValue < 200)
nrow(transect_issues)
# --> 5 events are affected after data update

transect_issues <- transect_issues %>%
  dplyr::select(id, eventID, sampleSizeValue, sampleSizeUnit, eventDate, locationID, locality, verbatimLocality, lineName, dataSet)

#write.csv(transect_issues, file = "Hønsefugl_TransectLength_Issues.csv", row.names = FALSE)

## Check for duplicate transects within years
transect_duplicates <- event %>%
  dplyr::filter(eventRemarks == "Line transect") %>%
  dplyr::group_by(locationID, locality, verbatimLocality, dataSet, lineName, year) %>%
  dplyr::summarise(transectCount = dplyr::n(), .groups = 'keep')

table(transect_duplicates$transectCount)
# --> Vast majority only has one entry per year (as it should be)
# --> 5 transects have two entries per year, one has 3

transect_duplicates <- transect_duplicates %>%
  dplyr::filter(transectCount > 1)

#write.csv(transect_duplicates, file = "Hønsefugl_TransectDuplicates_Issues.csv", row.names = FALSE)

transect_duplicates$flag <- "x"
transect_duplicateSuspects <- dplyr::right_join(event, transect_duplicates[,c("locationID", "year", "flag")], by = c("locationID", "year")) %>%
  dplyr::filter(eventRemarks == "Line transect")

event <- dplyr::left_join(event, transect_duplicates[,c("locationID", "year", "flag")], by = c("locationID", "year"))


#-------------------------------------------#
# SUMMARIZE INFORMATION AT DIFFERENT LEVELS #
#-------------------------------------------#

## Count willow ptarmigan observations per event
occ_sum <- occ %>%
  dplyr::filter(!is.na(eventID) & scientificName == "Lagopus lagopus") %>%
  dplyr::group_by(eventID) %>%
  dplyr::summarise(n_occ = dplyr::n())

table(occ_sum$n_occ)
# --> The majority of events have between 1 and 3 observations (few have 4)

## Add counts to event table
event <- event %>%
  dplyr::inner_join(occ_sum, by = "eventID")

## Summarize events per transect
transect <- event %>%
  dplyr::group_by(dataSet, verbatimLocality, locality, lineName) %>%
  dplyr::summarise(minYear = min(lubridate::year(eventDate)),
                   maxYear = max(lubridate::year(eventDate)),
                   n_years = max(lubridate::year(eventDate)) - min(lubridate::year(eventDate)) + 1,
                   n_occ = sum(n_occ), .groups = "keep")
hist(transect$n_years)
table(transect$n_years)
mean(transect$n_years)

hist(transect$n_occ)
table(transect$n_occ)
mean(transect$n_occ)

# --> There are 2149 transects in total
# --> Study duration ranges from 1 to 23 years, average = 8.57 years
# --> Number of occurrences/observations ranges from 1 to 417, average = 38.78


## Summarize transects per locality
locality <- transect %>%
  dplyr::group_by(dataSet, verbatimLocality, locality) %>%
  dplyr::summarise(n_lines = dplyr::n(),
                   minYear_tot = min(minYear),
                   maxYear_tot = max(maxYear),
                   minYear_all = max(minYear),
                   maxYear_all = min(maxYear),
                   n_years_tot = max(maxYear) - min(minYear) + 1,
                   n_occ = sum(n_occ), .groups = "keep")
hist(locality$n_lines)
table(locality$n_lines)
mean(locality$n_lines)

hist(locality$n_occ)
table(locality$n_occ)
mean(locality$n_occ)

# --> There are 120 localities with 1-86 lines each (average = 17.91)
# --> Number of occurrences/observations ranges from 2 to 5671, average = 694.5

write.csv(locality, "Data_AvailabilityExp/localities.csv", row.names = FALSE)

## Summarize localities by (reporting) area 
areas <- locality %>%
  dplyr::group_by(dataSet, verbatimLocality) %>%
  dplyr::summarise(n_localities = dplyr::n(),
                   n_lines = sum(n_lines),
                   minYear_tot = min(minYear_tot),
                   maxYear_tot = max(maxYear_tot),
                   minYear_all = max(minYear_all),
                   maxYear_all = min(maxYear_all),
                   n_years_tot = max(maxYear_tot) - min(minYear_tot) + 1,
                   n_occ = sum(n_occ), .groups = "keep")

hist(areas$n_localities)
table(areas$n_localities)
mean(areas$n_localities)

hist(areas$n_lines)
table(areas$n_lines)
mean(areas$n_lines)

hist(areas$n_occ)
table(areas$n_occ)
mean(areas$n_occ)

# --> There are 44 areas with 1-11 localities (average 2.72) and 1-126 lines each (average = 48.84)
# --> Number of occurrences/observations ranges from 1 to 8560, average = 1894.091

write.csv(areas, "Data_AvailabilityExp/areas.csv")

#-------------------------------------------------#
# VISUALIZE DATA AVAILABILITY AT DIFFERENT LEVELS #
#-------------------------------------------------#

library(ggplot2)
library(viridis)

## Assemble transect data including years
transect_y <- event %>%
  dplyr::mutate(Year = lubridate::year(eventDate),
                lineNo = stringi::stri_extract_first_regex(lineName, "[0-9]+")) %>%
  dplyr::group_by(dataSet, verbatimLocality, locality, lineName, lineNo, Year) %>%
  dplyr::summarise(n_occ = sum(n_occ), .groups = "keep")

## Plot data availability over time by area, locality, and transect line
areaNames <- unique(transect_y$verbatimLocality)
localityNames <- unique(transect_y$locality)


for(a in 1:length(areaNames)){
  
  pdf(paste0('Data_AvailablityExp/DataAvailability_', areaNames[a], '.pdf'), height = 8, width = 5)
  
  for(l in 1:length(localityNames)){
    
    transect_y_sub <- subset(transect_y, verbatimLocality == areaNames[a] & locality == localityNames[l])
    
    if(nrow(transect_y_sub) > 0){
      print(
        ggplot(transect_y_sub, aes(x = Year, y = lineNo)) +
          #geom_line(aes(colour = as.factor(lineNo))) +
          geom_line(colour = "grey70") +
          #geom_point(aes(size = n_occ, colour = as.factor(lineNo))) + 
          geom_point(aes(colour = n_occ), size = 1.5) + 
          scale_x_continuous(limits = c(min(transect_y$Year), max(transect_y$Year)), breaks = min(transect_y$Year):max(transect_y$Year)) + 
          scale_color_viridis() + 
          ggtitle(paste0(areaNames[a], ': ', localityNames[l])) + 
          theme_bw() + 
          theme(panel.grid.minor.x = element_blank(), axis.text.x = element_text(angle = 90))
      )
    }
  }
  dev.off()
}
