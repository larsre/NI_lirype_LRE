#' List the areas that should be used in the analysis.
#' Areas are defined from 'Reporting levels', i.e. the spatial extent of the
#' area used for local density estimation
#' 
#' The current list contains valid areas (reporting levels) for the three data
#' sets openly available through GBIF
#'
#' @return a character vector with UTF-8 encoding. 
#' @export
#'
#' @examples
#'     listAreas()

listAreas <- function(){
  
  areas <- c(
    "Åfjord",
    "Ålen og Haltdalen Fjellstyre",
    "Budal Fjellstyre",
    "Dovre Fjellstyre",
    "Eidfjord Fjellstyre",
    "Engerdal Fjellstyre",
    "Fron Fjellstyre",
    "Gausdal Fjellstyre",
    "Helgeland 1",
    "Helgeland 2",
    "Helgeland 3",
    "Indre Finnmark",
    "Indre Salten",
    "Kongsvoll",
    #"Kvikne Fjellstyre", #Now part of region Nord-Østerdal
    "Lesja Fjellstyre",
    "Lierne Fjellstyre",
    "Midtre Salten",
    "Namdalseid Fjellstyre",
    "Namskogan Fjellstyre",
    "Njardarheim",
    "Nordre Salten",
    #"Os Fjellstyre", #Now part of region Nord-Østerdal
    "Osen Fjellstyre",
    #"Regionfelt Østlandet", #exclude from analysis - very few observations/sporadic counting last years
    "Ringebu Fjellstyre",
    "Røyrvik Fjellstyre",
    "Selbu Fjellstyre",
    "Snåsa Fjellstyre",
    "Soknedal Fjellstyre",
    "Sollia Fjellstyre",
    "Statskog Bangdalen/Klinga utmarkslag Rype", #Formerly 'Statskog og Klinga utm.'
    "Statskog Røros",
    "Steinkjer Lirype",
    "Stjørdal Fjellstyre",
    "Troms Ytre Nord",
    "Vest Finnmark kyst",
    "Vestre Slidre Fjellstyre",
    "Ytre Salten",
    #"Troms Midt",
    #"Troms Nord",
    #"Troms Sør",
    "Øst Finnmark",
    "Øvre Numedal Fjellstyre",
    "Øystre Slidre Fjellstyre"
    #added 2021:
    #"Høylandet Fjellstyre",
    #"Øyer Fjellstyre",
    #added 2022:
    #"Folldal Fjellstyre" #Now part of region Nord-Østerdal
    #added 2023:
    #"Sunndal Fjellstyre",
    #"Rendalen",
    #"Statskog Bangdalen/Klinga utmarkslag Rype" #new name for 'Statskog og Klinga utm.' from 2023
  )
  
  Encoding(areas) <- "UTF-8"
  return(areas)
}

