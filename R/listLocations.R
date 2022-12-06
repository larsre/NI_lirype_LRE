#' List the locations that should be used in the analysis
#'
#' @return a character vector with UTF-8 encoding. 
#' @export
#'
#' @examples

listLocations <- function(){
  
  localities <- c(
    "Lierne Fjellst. Vest",
    "Lierne Fjellst. øst",
    "Middagskneppen"
  )
  
  Encoding(localities) <- "UTF-8"
  return(localities)
}

