#' Download ptarmigan line transect data from Living Norway Data Portal
#'
#' @param version numeric. Version of the dataset to download. 
#' @param save logical. If TRUE (default), saves downloaded .zip archive into 
#' "data" folder. Folder is created if it does not exist yet.
#'
#' @return Darwin Core Archive (class DwCArchive, R6) containing line transect 
#' data for willow ptarmigan (Lagopus lagopus) in Norway. 
#' @export
#'
#' @examples

downloadLN <- function(version, save = TRUE){
  
  ## Define dataset key
  datasetKey <- "b49a2978-0e30-4748-a99f-9301d17ae119"
  
  ## Download specific data version from Living Norway
  Rype_arkiv <- LivingNorwayR::getLNportalData(datasetKey = datasetKey, version = version)
  
  ## Optional: Save the archive file to ...data/
  if(save){
    dir.create("data", showWarnings = FALSE)
    Rype_arkiv$exportAsDwCArchive(file.path("data/Rype_arkiv.zip"))
  }
  
  ## Return the archive
  return(Rype_arkiv)
}