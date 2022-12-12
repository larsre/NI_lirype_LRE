#' List the areas that should be used in the analysis
#'
#' @return a character vector with UTF-8 encoding. 
#' @export
#'
#' @examples

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
    "Kvikne Fjellstyre",
    "Lesja Fjellstyre",
    "Lierne Fjellstyre",
    "Namdalseid Fjellstyre",
    "Namskogan Fjellstyre",
    "Os Fjellstyre",
    "Osen Fjellstyre",
    "Øvre Numedal Fjellstyre",
    "Øystre Slidre Fjellstyre",
    "Ringebu Fjellstyre",
    "Røyrvik Fjellstyre",
    "Selbu Fjellstyre",
    "Snåsa Fjellstyre",
    "Soknedal Fjellstyre",
    "Sollia Fjellstyre",
    "Steinkjer Lirype",
    "Stjørdal Fjellstyre",
    "Vestre Slidre Fjellstyre",
    #"Gjerstad",
    "Helgeland 1",
    "Helgeland 2",
    "Helgeland 3",
    "Indre Salten",
    "Indre Troms",
    "Kongsvoll",
    "Midtre Salten",
    "Njardarheim",
    "Nordre Salten",
    "Regionfelt Østlandet",
    "Statskog og Klinga utm.",
    "Statskog Røros",
    "Ytre Salten",
    "Ytre Troms",
    "Indre Finnmark",
    "Øst Finnmark",
    "Vest Finnmark kyst"
  )
  
  Encoding(areas) <- "UTF-8"
  return(areas)
}

