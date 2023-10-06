#' Assign each transect line to county by using the geographical data available in the data set
#'
#' @param df_trans a tibble with the wrangled line transect data
#'
#' @return a data frame with the assigned county name for each 'LineName' in the transect data set
#' @export
#'
#' @examples

assignCounty <- function(df_trans = LT_data$d_trans) {
  #load required libraries
  if(!require(sf)) { print("Package 'sf' not installed"); break; }
  if(!require(dplyr)) { print("Package 'dplyr' not installed"); break; }
  
  #convert WKT to sf objects
  wgs4326 <- st_as_sfc(df_trans$footprintWKT, crs = 4326)

  #import county map
  countyMap <- st_read("maps/counties2020.shp")
  Encoding(countyMap$name) <- "latin1"
  
  #convert map attributes to data frame
  county.names <- as.data.frame(countyMap)
  county.names <- county.names %>% dplyr::select(OBJECTID, name)

  #extract the county for each line transect
  countyIntersect <- st_intersects(wgs4326, countyMap, sparse = T)
  countyIntersect <- as.data.frame(countyIntersect)
  countyIntersect <- dplyr::left_join(countyIntersect, county.names, by = c("col.id" = "OBJECTID"))
  colnames(countyIntersect)[3] <- "county"
  
  #assign and return
  df_trans$county <- countyIntersect$county
  return(df_trans)
}