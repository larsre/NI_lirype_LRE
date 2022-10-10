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
#' @param localities vector of strings listing areas to consider
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

prepareInputData <- function(d_trans, d_obs, d_cmr, localities, dataVSconstants = TRUE, save = TRUE){
  
  # Multi-area setup #
  #------------------#
  
  ## Number of areas
  N_areas <- length(localities)
  
  ## Variables shared across areas
  # Number of age classes
  N_ageC <- 2
  
  # Truncation distance
  W <- 200
  
  # Scaling parameter
  scale1 <- 1000
  
  ## Count sites per area (for defining correct array dimensions)
  site_count <- d_trans %>%                             
    dplyr::group_by(locality) %>%
    dplyr::summarise(count = n_distinct(locationID))
  
  ## Count observations per area
  obs_count <- d_obs %>% 
    dplyr::filter(between(DistanceToTransectLine, -0.1, W)) %>% 
    dplyr::group_by(locality) %>%
    dplyr::summarise(count = n())
    
  ## Set up arrays for are-specific data
  N_sites <- N_years <- N_obs <- N_R_obs <- rep(NA, N_areas)
  
  A <- matrix(NA, nrow = N_areas, ncol = n_distinct(d_trans$Year))
  y <- Year_obs <- zeros_dist <- matrix(NA, nrow = N_areas, ncol = max(obs_count$count))
  R_obs <- R_obs_year <- matrix(NA, nrow = N_areas, ncol = max(obs_count$count))
  
  L <- N_line_year <- array(NA, dim = c(N_areas, max(site_count$count), n_distinct(d_trans$Year)))
  
  N_a_line_year <- array(NA, dim = c(N_areas, N_ageC, max(site_count$count), n_distinct(d_trans$Year)))
  
  
  for(x in 1:N_areas){
    
    ## Subset data (specific area)
    d_trans_sub <- subset(d_trans, locality == localities[x])
    d_obs_sub <- subset(d_obs, locality == localities[x])
    
    # Constants #
    #-----------#
    
    ## Numbers of years and sites
    N_sites[x] <- n_distinct(d_trans_sub$locationID)
    N_years[x] <- n_distinct(d_trans_sub$Year)
    
    
    # Transect characteristics #
    #--------------------------#
    
    ## Tansect lengths
    TransLen <- d_trans_sub %>% 
      dplyr::select(locationID, Year, sampleSizeValue) %>%
      dplyr::mutate(sampleSizeValue = sampleSizeValue/scale1) %>%
      reshape2::dcast(locationID~Year, value.var = "sampleSizeValue", sum) %>%
      arrange(locationID)
    
    L[x, 1:N_sites[x], 1:N_years[x]] <- TransLen %>% select(-locationID) %>% as.matrix() 
    
    ## Total covered area 
    A_prel <- TransLen %>% 
      dplyr::select(-locationID) 

    A[x, 1:N_years[x]] <- colSums(A_prel)*(W/scale1)*2
    
    
    # Observation distance from transect #
    #------------------------------------#
    
    temp_dist <- d_obs_sub %>% 
      dplyr::filter(between(DistanceToTransectLine, -0.1, W)) %>% 
      dplyr::select(Year, DistanceToTransectLine) %>%
      dplyr::mutate(Year2 = Year - (min(Year)) + 1)
    
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
    
    N_line_year[x, 1:N_sites[x], 1:N_years[x]] <- TaksObs %>% 
      dplyr::select(-locationID) %>% 
      as.matrix() 
    
    
    # Number of birds/line (by age class) #
    #-------------------------------------#
    
    ## Juveniles (& unknowns)
    TaksObs_J <- d_obs_sub %>% filter(between(DistanceToTransectLine, -0.1, W)) %>%
      dplyr::mutate(cs = unknownJuvenile + unknownunknown) %>%
      reshape2::dcast(locationID~Year, value.var="cs", sum) %>%
      dplyr::right_join(., temp, by=c("locationID"="locationID")) %>%
      replace(., is.na(.), 0) %>%
      dplyr::arrange(locationID)
    
    N_J_line_year <- TaksObs_J %>% 
      dplyr::select(-locationID) %>% 
      as.matrix() 
    
    ## Adults 
    TaksObs_A <- d_obs_sub %>% filter(between(DistanceToTransectLine, -0.1, W)) %>%
      dplyr::mutate(cs = FemaleAdult + MaleAdult) %>%
      reshape2::dcast(locationID~Year, value.var="cs", sum) %>%
      dplyr::right_join(., temp, by=c("locationID"="locationID")) %>%
      replace(., is.na(.), 0) %>%
      dplyr::arrange(locationID)
    
    N_A_line_year <- TaksObs_A %>% 
      dplyr::select(-locationID) %>% 
      as.matrix() 
    
    ## Check and combine in array
    if(!(all(N_J_line_year + N_A_line_year == N_line_year[x, 1:N_sites[x], 1:N_years[x]]))){
      warning("Number of observed adults and juveniles does not add up correctly. Double-check data.")
    }
    
    N_a_line_year[x, 1, 1:N_sites[x], 1:N_years[x]] <- N_J_line_year
    N_a_line_year[x, 2, 1:N_sites[x], 1:N_years[x]] <- N_A_line_year
    
    
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
  R_obs <- R_obs[1:N_areas, 1:max(N_R_obs)]
  R_obs_year <- R_obs_year[1:N_areas, 1:max(N_R_obs)]
  
  
  # Data assembly #
  #---------------#
  
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
    
    N_years = N_years, # Number of years with data
    N_sites = N_sites, # Total number of monitored sites
    
    A = A, # Total covered area per year
    W = W, # Truncation distance
    scale1 = scale1, # Scaling parameter
    N_ageC = N_ageC, # Number of age classes
    
    Survs1 = d_cmr$Survs1, # Season 1 releases & survivors (area 1)
    Survs2 = d_cmr$Survs2, # Season 2 releases & survivors (area 1)
    
    area_names = localities
  )
  
  ## Assembling Nimble data
  nim.data <- list(R_obs = input.data$R_obs, y = input.data$y, 
                   zeros.dist = input.data$zeros_dist, L = input.data$L, 
                   N_line_year = input.data$N_line_year, 
                   N_a_line_year = input.data$N_a_line_year, 
                   A = input.data$A,
                   Survs1 = input.data$Survs1, Survs2 = input.data$Survs2)
  
  ## Assembling Nimble constants
  nim.constants <- list(N_years = input.data$N_years, W = input.data$W, scale1 = scale1,
                        N_obs = input.data$N_obs, Year_obs = input.data$Year_obs,
                        N_sites = input.data$N_sites, 
                        R_obs_year = input.data$R_obs_year, N_R_obs = input.data$N_R_obs,
                        N_ageC = N_ageC,
                        area_names = input.data$area_names)
  
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



