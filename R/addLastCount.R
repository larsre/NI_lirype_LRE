#' Add data from years missing in GBIF due to quarantine, i.e. most recent year.
#'
#' @param df a list of 2 tibbles with the wrangled line transect and occurrence
#' data
#'
#' @return a list of 2 tibbles with the added missing data.
#' 
#' @export
#'
#' @examples
#'      addLastCount(LT_data)

addLastCount <- function(df = LT_data) {
  # import missing data
  dfLast <- readRDS("data/countData2024.rds")
  
  # create reference list for locationID
  refLoc <- df$d_trans %>% group_by(LineName, locationID) %>% ungroup() %>% distinct_at(vars(c(LineName, locationID)), .keep_all = F)
  refLoc <- refLoc[order(refLoc$LineName),]
  
  ## format line transect data
  dfl_trans <- subset(dfLast$d_trans24, select = c("TakseringID","LengdeTaksert","FK_Fylkekomnr","Kommunenavn","OmradeNavn","Rapporteringsniva","FK_LinjeID","Aar","STAsText"))
  dfl_trans$FK_LinjeID <- as.character(dfl_trans$FK_LinjeID)
  dfl_trans <- left_join(dfl_trans, refLoc, by = c("FK_LinjeID" = "LineName"))
  
  # remove line transects that are not available in the original data
  dfl_trans <- subset(dfl_trans, !is.na(dfl_trans$locationID))
  
  # rename columns
  colnames(dfl_trans) <- c("eventID","sampleSizeValue","stateProvince","municipality","locality","verbatimLocality","LineName","Year","footprintWKT","locationID")
  dfl_trans$eventDate <- NA
  dfl_trans$modified <- NA
  dfl_trans$samplingProtocol <- "Distance sampling based on line transects"
  dfl_trans$eventRemarks <- "Line transect"
  dfl_trans <- subset(dfl_trans, select = c("locationID","eventDate","eventID","modified","samplingProtocol","eventRemarks","sampleSizeValue","stateProvince","municipality","locality","verbatimLocality","LineName","Year","footprintWKT"))

  
  ## format occurrence data
  dfLast$d_obs24$scientificName <- "Lagopus lagopus"
  dfl_obs <- subset(dfLast$d_obs24, select = c("ObservasjonId","scientificName","AntallHunn","AntallHann","AntallKylling","AntallUkjent","FK_LinjeID","OmradeNavn","Rapporteringsniva","TakseringID","LinjeAvstand","Aar.x"))
  dfl_obs$FK_LinjeID <- as.character(dfl_obs$FK_LinjeID)
  dfl_obs <- left_join(dfl_obs, refLoc, by = c("FK_LinjeID" = "LineName"))
  
  # remove occurrences which are not attached to a line transect available in the original data
  dfl_obs <- subset(dfl_obs, !is.na(dfl_obs$locationID))
  
  # rename columns
  dfl_obs <- dfl_obs[,c(1:6,13,8:12)]
  colnames(dfl_obs) <- c("eventID","scientificName","FemaleAdult","MaleAdult","unknownJuvenile","unknownunknown","locationID","locality","verbatimLocality","parentEventID","DistanceToTransectLine","Year")
  
  
  # merge with original data
  df$d_trans <- rbind(df$d_trans, dfl_trans)
  df$d_obs <- rbind(df$d_obs, dfl_obs)
  
  return(df)
}