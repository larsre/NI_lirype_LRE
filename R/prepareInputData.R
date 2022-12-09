#' Prepare line transect and known fate CMR data for integrated analysis
#'
#' @param d_trans tibble containing information on transects (events). Output of
#' wrangleData_LineTrans(). 
#' @param d_obs tibble containing information on observations made along transects 
#' (distance to transect line, numbers of birds in each age/sex class observed,
#' etc.). Output of wrangleData_LineTrans(). 
#' @param d_cmr list with 2 elements. Surv1 and Surv2 are matrices of individuals 
#' released (column 1) and known to have survived (column 2) in each year (row)
#' for season 1 and season 2, respectively. Output of wrangleData_CMR().
#' @param localities vector of strings listing localities to consider. Either localities or areas must be provided. 
#' @param areas vector of strings listing areas to consider. Either localities or areas must be provided. 
#' @param areaAggregation logical. If TRUE, areas are used as smallest spatial unit. If FALSE, locations (within areas) are used as smallest spatial unit.
#' @param dataVSconstants logical. If TRUE (default) returns a list of 2 lists
#' containing data and constants for analysis with Nimble. If FALSE, returns a
#' list containing all data and constants. 
#' @param save logical. If TRUE (default) saves prepared data in working 
#' directory as .rds file.
#'
#' @return A list or list of lists, depending on argument `dataVSconstants` 
#' (see above).
#' @export
#'
#' @examples

prepareInputData <- function(d_trans, d_obs, d_cmr, localities = NULL, areas = NULL, areaAggregation, dataVSconstants = TRUE, save = TRUE){
  
  # Multi-area setup #
  #------------------#
  
  ## Assignment of spatial units
  if(areaAggregation){
    sUnits <- areas
  }else{
    sUnits <- localities
  }
  N_sUnits <- length(sUnits)
  
  ## Rename appropriate column in line transect data to reflect level of spatial aggregation
  if(areaAggregation){
    colnames(d_trans)[which(colnames(d_trans) == "verbatimLocality")] <- "spatialUnit"
    colnames(d_obs)[which(colnames(d_obs) == "verbatimLocality")] <- "spatialUnit"
  }else{
    colnames(d_trans)[which(colnames(d_trans) == "locality")] <- "spatialUnit"
    colnames(d_obs)[which(colnames(d_obs) == "locality")] <- "spatialUnit"
  }
  
  ## Variables shared across areas
  # Number of age classes
  N_ageC <- 2
  
  # (Number of) years
  range_yearsTot <- min(d_trans$Year):max(d_trans$Year)
  idx_yearsTot <- range_yearsTot - min(d_trans$Year) + 1
  N_yearsTot <- length(range_yearsTot)
  
  # Truncation distance
  W <- 200
  
  # Scaling parameter
  scale1 <- 1000
  
  ## Count sites per area (for defining correct array dimensions)
  site_count <- d_trans %>%                             
    dplyr::group_by(spatialUnit) %>%
    dplyr::summarise(count = n_distinct(locationID))
  
  ## Count observations per area
  obs_count <- d_obs %>% 
    dplyr::filter(between(DistanceToTransectLine, -0.1, W)) %>% 
    dplyr::group_by(spatialUnit) %>%
    dplyr::summarise(count = n())
    
  ## Set up arrays for area-specific data
  N_sites <- N_years <- min_years <- max_years <- N_obs <- N_R_obs <- rep(NA, N_sUnits)
  
  A <- matrix(NA, nrow = N_sUnits, ncol = N_yearsTot)
  y <- Year_obs <- zeros_dist <- matrix(NA, nrow = N_sUnits, ncol = max(obs_count$count))
  R_obs <- R_obs_year <- matrix(NA, nrow = N_sUnits, ncol = max(obs_count$count))
  
  L <- N_line_year <- N_J_line_year <- N_A_line_year <- array(0, dim = c(N_sUnits, max(site_count$count), N_yearsTot))
  
  N_a_line_year <- array(0, dim = c(N_sUnits, N_ageC, max(site_count$count), N_yearsTot))
  
  
  for(x in 1:N_sUnits){
    
    ## Subset data (specific area)
    # if(areaAggregation){
    #   d_trans_sub <- subset(d_trans, spatialUnit == areas[x])
    #   d_obs_sub <- subset(d_obs, spatialUnit == areas[x])
    # }else{
    #   d_trans_sub <- subset(d_trans, spatialUnit == localities[x])
    #   d_obs_sub <- subset(d_obs, spatialUnit == localities[x])
    # }
    d_trans_sub <- subset(d_trans, spatialUnit == sUnits[x])
    d_obs_sub <- subset(d_obs, spatialUnit == sUnits[x])
    
    # Constants #
    #-----------#
    
    ## Numbers of sites
    N_sites[x] <- n_distinct(d_trans_sub$locationID)
    
    ## (Number of) years with data
    years_mon <- sort(unique(d_trans_sub$Year)) - min(range_yearsTot) + 1
    min_years[x] <- min(years_mon)
    max_years[x] <- max(years_mon)

    # Transect characteristics #
    #--------------------------#
    
    ## Tansect lengths
    TransLen <- d_trans_sub %>% 
      dplyr::select(locationID, Year, sampleSizeValue) %>%
      dplyr::mutate(sampleSizeValue = sampleSizeValue/scale1) %>%
      reshape2::dcast(locationID~Year, value.var = "sampleSizeValue", sum) %>%
      arrange(locationID)
    
    TransLen_mat <- TransLen %>% select(-locationID) %>% as.matrix()
    colnames(TransLen_mat) <- as.numeric(colnames(TransLen_mat)) - min(range_yearsTot) + 1
    
    for(t in 1:N_yearsTot){
      
      if(t %in% years_mon){
        L[x, 1:N_sites[x], t] <- TransLen_mat[1:N_sites[x], which(years_mon == t)]
      }
    }
 
    
    ## Total covered area 
    A[x, 1:N_yearsTot] <- colSums(L[x,,])*(W/scale1)*2
    
    
    # Observation distance from transect #
    #------------------------------------#
    
    temp_dist <- d_obs_sub %>% 
      dplyr::filter(between(DistanceToTransectLine, -0.1, W)) %>% 
      dplyr::select(Year, DistanceToTransectLine) %>%
      dplyr::mutate(Year2 = Year - (minYear) + 1)
    
    ## Number of observations
    N_obs[x] <- length(temp_dist$DistanceToTransectLine)
    
    ## Distance to transect line
    y[x, 1:N_obs[x]] <- temp_dist$DistanceToTransectLine
    
    ## Observation year
    Year_obs[x, 1:N_obs[x]] <- temp_dist$Year2
    
    ## Vector of 0's of same length as y
    zeros_dist[x, 1:N_obs[x]] <- rep(0, N_obs[x])
    
    
    # Number of birds/line (pooled age classes) #
    #-------------------------------------------#
    
    temp <- TransLen %>% select(locationID)
    
    TaksObs <- d_obs_sub %>% filter(between(DistanceToTransectLine, -0.1, W)) %>%
      dplyr::mutate(cs = unknownJuvenile+unknownunknown+FemaleAdult+MaleAdult) %>%
      reshape2::dcast(locationID~Year, value.var = "cs", sum) %>%
      dplyr::right_join(., temp, by = c("locationID" = "locationID")) %>%
      replace(., is.na(.), 0) %>%
      dplyr::arrange(locationID)
    
    TaksObs_mat <- TaksObs %>% dplyr::select(-locationID) %>% as.matrix() 
    colnames(TaksObs_mat) <- as.numeric(colnames(TaksObs_mat)) - min(range_yearsTot) + 1
    
    for(t in 1:N_yearsTot){
      
      if(t %in% years_mon){
        N_line_year[x, 1:N_sites[x], t] <- TaksObs_mat[1:N_sites[x], which(years_mon == t)]
      }
    }
    
    # Number of birds/line (by age class) #
    #-------------------------------------#
    
    ## Juveniles (& unknowns)
    TaksObs_J <- d_obs_sub %>% filter(between(DistanceToTransectLine, -0.1, W)) %>%
      dplyr::mutate(cs = unknownJuvenile + unknownunknown) %>%
      reshape2::dcast(locationID~Year, value.var="cs", sum) %>%
      dplyr::right_join(., temp, by=c("locationID"="locationID")) %>%
      replace(., is.na(.), 0) %>%
      dplyr::arrange(locationID)
    
    TaksObs_J_mat <- TaksObs_J %>% dplyr::select(-locationID) %>% as.matrix() 
    colnames(TaksObs_J_mat) <- as.numeric(colnames(TaksObs_J_mat)) - min(range_yearsTot) + 1
    
    for(t in 1:N_yearsTot){
      
      if(t %in% years_mon){
        N_J_line_year[x, 1:N_sites[x], t] <- TaksObs_J_mat[1:N_sites[x], which(years_mon == t)]
      }
    }
    
    ## Adults 
    TaksObs_A <- d_obs_sub %>% filter(between(DistanceToTransectLine, -0.1, W)) %>%
      dplyr::mutate(cs = FemaleAdult + MaleAdult) %>%
      reshape2::dcast(locationID~Year, value.var="cs", sum) %>%
      dplyr::right_join(., temp, by=c("locationID"="locationID")) %>%
      replace(., is.na(.), 0) %>%
      dplyr::arrange(locationID)
    
    TaksObs_A_mat <- TaksObs_A %>% dplyr::select(-locationID) %>% as.matrix() 
    colnames(TaksObs_A_mat) <- as.numeric(colnames(TaksObs_A_mat)) - min(range_yearsTot) + 1
    
    for(t in 1:N_yearsTot){
      
      if(t %in% years_mon){
        N_A_line_year[x, 1:N_sites[x], t] <- TaksObs_A_mat[1:N_sites[x], which(years_mon == t)]
      }
    }
    
    ## Check and combine in array
    if(!(all(N_J_line_year + N_A_line_year == N_line_year))){
      warning("Number of observed adults and juveniles does not add up correctly. Double-check data.")
    }
    
    N_a_line_year[x, 1,,] <- N_J_line_year[x,,]
    N_a_line_year[x, 2,,] <- N_A_line_year[x,,]
    
    
    # Recruitment #
    #-------------#
    
    temp_Rec <- d_obs_sub %>% filter(between(DistanceToTransectLine, -0.1, W)) %>%
      dplyr::mutate(R = unknownJuvenile + unknownunknown) %>%
      dplyr::mutate(Maletemp = unknownJuvenile + unknownunknown + FemaleAdult, MaleIndeks = if_else(Maletemp == 0, 1, 0)) %>%
      dplyr::select(Year, R, MaleIndeks) %>%
      dplyr::mutate(Year2 = Year - (min(Year)) + 1) %>%
      dplyr::filter(MaleIndeks == 0) # --> drop all observations of only males
    
    N_R_obs[x] <- length(temp_Rec$R)
    R_obs[x, 1:N_R_obs[x]] <- temp_Rec$R
    R_obs_year[x, 1:N_R_obs[x]] <- temp_Rec$Year2

  } 
  
  ## Drop excess columns in recruitment data
  R_obs <- R_obs[1:N_sUnits, 1:max(N_R_obs)]
  R_obs_year <- R_obs_year[1:N_sUnits, 1:max(N_R_obs)]
  
  
  # Data assembly #
  #---------------#
  
  ## Set spatial index for CMR data
  if(areaAggregation){
    SurvAreaIdx <- which(sUnits == d_cmr$area_names)
  }else{
    SurvAreaIdx <- which(sUnits == d_cmr$locality_names)
  }
  
  ## Assembling all data in a list
  input.data <- list(
    R_obs = R_obs, # Observed numbers of recruits
    R_obs_year = R_obs_year, # Year of observed numbers of recruits
    N_R_obs = N_R_obs, # Total number of observations of numbers of recruits
    
    y = y, # Distance to transect line for each individual observation
    zeros_dist = zeros_dist, # Vector of 0's of same length as y
    Year_obs = Year_obs, # Year of each observation
    N_obs = N_obs, # Total number of observations
    
    N_line_year = N_line_year, # Number of birds observed per site per year
    N_a_line_year = N_a_line_year, # Number of birds observed per ageclass per site per year
    L = L, # Transect length per site and year
    
    N_years = N_yearsTot, # Max number of years with data
    min_years = min_years, # Earliest year (index) of monitoring per area
    max_years = max_years, # Last year (index) of monitoring per area
    N_sites = N_sites, # Total number of monitored sites per area
    
    A = A, # Total covered area per year
    W = W, # Truncation distance
    scale1 = scale1, # Scaling parameter
    N_ageC = N_ageC, # Number of age classes
    
    Survs1 = d_cmr$Survs1, # Season 1 releases & survivors (area 1)
    Survs2 = d_cmr$Survs2, # Season 2 releases & survivors (area 1)
    SurvAreaIdx = SurvAreaIdx,
    
    N_areas = N_sUnits,
    area_names = sUnits
  )
  
  ## Assembling Nimble data
  nim.data <- list(R_obs = input.data$R_obs, y = input.data$y, 
                   zeros.dist = input.data$zeros_dist, L = input.data$L, 
                   N_line_year = input.data$N_line_year, 
                   N_a_line_year = input.data$N_a_line_year, 
                   A = input.data$A,
                   Survs1 = input.data$Survs1, Survs2 = input.data$Survs2)
  
  ## Assembling Nimble constants
  nim.constants <- list(N_years = input.data$N_years, min_years = input.data$min_years, max_years = input.data$max_years,
                        W = input.data$W, scale1 = scale1,
                        N_obs = input.data$N_obs, Year_obs = input.data$Year_obs,
                        N_sites = input.data$N_sites, 
                        R_obs_year = input.data$R_obs_year, N_R_obs = input.data$N_R_obs,
                        N_ageC = N_ageC,
                        N_areas = input.data$N_areas, area_names = input.data$area_names,
                        SurvAreaIdx = input.data$SurvAreaIdx)
  
  ## Make final data list to return
  if(dataVSconstants){
    rype.data <- list(nim.data = nim.data,
                      nim.constants = nim.constants)
  }else{
    rype.data <- input.data
  }
  
  ## Optional: save data as .rds
  if(save){
    saveRDS(rype.data, file = "RypeData_forIM.rds")
  }
  
  ## Return data
  return(rype.data)
}



