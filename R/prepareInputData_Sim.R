prepareInputData_Sim <- function(SimData, addDummyDim){
  
  ## Add dummy dimensions to observational data if requested
  if(addDummyDim){
    
    N_a_line_year <- addDummyDimension(SimData$DS.data$DS.count)
    L <- addDummyDimension(SimData$DS.data$L)
    y <- addDummyDimension(SimData$DS.data$d)
    Year_obs <- addDummyDimension(SimData$DS.data$d_year)
    sumR_obs <- addDummyDimension(SimData$Rep.data$sumR_obs)
    sumAd_obs <- addDummyDimension(SimData$Rep.data$sumAd_obs)
    sumR_obs_year <- addDummyDimension(SimData$Rep.data$sumR_obs_year)
    N_sites <- c(SimData$SimParams$Jmax, NA)
    N_obs <- c(length(SimData$DS.data$d), NA)
    N_sumR_obs <- c(SimData$Rep.data$N_sumR_obs, NA)
  
  }else{
    
    N_a_line_year <- SimData$DS.data$DS.count
    L <- SimData$DS.data$L
    y <- SimData$DS.data$d
    Year_obs <- SimData$DS.data$d_year
    sumR_obs <- SimData$Rep.data$sumR_obs
    sumAd_obs <- SimData$Rep.data$sumAd_obs
    sumR_obs_year <- SimData$Rep.data$sumR_obs_year
    N_sites <- SimData$SimParams$Jmax
    N_obs <- length(SimData$DS.data$d)
    N_sumR_obs <- SimData$Rep.data$N_sumR_obs
  }

  
  ## Reformat data into vector/array list for analysis with Nimble
  input_data <- list(
    nim.data = list(
      sumR_obs = sumR_obs,
      sumAd_obs = sumAd_obs,
      y = y,
      L = L,
      N_a_line_year = N_a_line_year,
      Survs1 = SimData$RT.data$Survs1,
      Survs2 = SimData$RT.data$Survs2
    ),
    
    nim.constants = list(
      N_areas = 1,
      SurvAreaIdx = 1,
      N_years = SimData$SimParams$Tmax,
      year_Survs = SimData$RT.data$year_Survs,
      N_years_RT = SimData$RT.data$N_years_RT,
      W = SimData$SimParams$W,
      N_obs = N_obs,
      Year_obs = Year_obs,
      N_sites = N_sites,
      sumR_obs_year = sumR_obs_year,
      N_sumR_obs = N_sumR_obs,
      N_ageC = SimData$SimParams$Amax
    )
  )
  
  ## Return data
  return(input_data)
}