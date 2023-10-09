#' List counties to be included in the analysis
#' 
#' The current list contains counties with an adequate number of transect lines
#' and observations.
#'
#' @return a character vector with UTF-8 encoding. 
#' @export
#'
#' @examples

listCounties <- function(){
  
  counties <- c(
    #"Agder",
    #"Rogaland",
    #"Vestland",
    #"Viken",
    #"Møre og Romsdal",
    #"Vestfold og Telemark",
    "Innlandet",
    "Nordland",
    "Troms og Finnmark",
    "Trøndelag"
  )
  
  Encoding(counties) <- "UTF-8"
  return(counties)
}
