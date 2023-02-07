#' Find and visualize area-years with (suspiciously) high detection probability estimates
#'
#' With the way the model is currently written, the distance sampling detection 
#' probability (p) is constrained to be 1 whenever the estimated effective strip
#' width (esw) is higher than the truncation distance (W). This can lead to a 
#' select few area-years with posterior distributions of p that have a 
#' substantial portion of 1's. This may or may not indicate an issue (violation
#' of model assumptions), and this is a helper function that identifies the
#' affected area-years and optionally visualized the data distribution for them. 
#' 
#' @param mcmc.out an mcmc list containing posterior samples from a model run.
#' @param N_areas integer. Number of areas included in analyses. 
#' @param area_names character vector containing area/location names.
#' @param N_years integer. Number of years included in analyses. 
#' @param minYear integer. The first year considered in analyses. 
#' @param threshold numeric. The threshold detection probability to consider as (too) high.
#' @param quantile numeric. The quantile of the posterior that needs to be above the threshold for getting flagged.
#' @param plot logical. If TRUE (default), observation distance distribution is plotted for flagged area-years.
#' @param obsDist matrix of area-specific detection distances.
#' @param obsDist_year matrix of year indices of area-specific detection distances.
#'
#' @return
#' @export
#'
#' @examples

findHighP <- function(mcmc.out, N_areas, area_names, N_years, minYear, threshold, quantile, plot = TRUE, obsDist, obsDist_year){
  
  ## List relevant parameters
  idxGrid <- expand.grid(c(1:N_areas), c(1:N_years))
  
  relParams <- paste0("p[", idxGrid[,1], ", ", idxGrid[,2], "]")
  
  ## Reduce posterior samples to relevant parameters only
  mcmc.out <- mcmc.out[, which(colnames(mcmc.out[[1]]) %in% relParams), drop = TRUE]
  
  ## Reformat samples
  out.data <- reshape2::melt(as.matrix(mcmc.out))
  colnames(out.data) <- c("Sample", "Parameter", "Value")
  
  ## Summarise posteriors (selected quantile) and filter using threshold
  sum.data <- out.data %>%
    dplyr::group_by(Parameter) %>%
    dplyr::summarise(quantile = quantile(Value, probs = quantile)) %>%
    dplyr::filter(quantile >= threshold)
  
  ## Extract area and year indices
  idx.data <- stringi::stri_extract_all_regex(sum.data$Parameter, "[0-9]+", simplify = TRUE)
  colnames(idx.data) <- c("AreaIdx", "YearIdx")
  idx.data <- data.frame(idx.data)
  idx.data$AreaIdx <- as.numeric(idx.data$AreaIdx)
  idx.data$YearIdx <- as.numeric(idx.data$YearIdx)
  idx.data$Area <- area_names[idx.data$AreaIdx]
  idx.data$Year <- idx.data$YearIdx + minYear - 1
  
  ## Assemble and count observations
  idx.data$ObsCount <- NA
  area_yr_obs <- list()
  for(i in 1:nrow(idx.data)){
    area_yr_obs[[i]] <- obsDist[idx.data$AreaIdx[i], which(obsDist_year[idx.data$AreaIdx[i],] == idx.data$YearIdx[i])]
    idx.data$ObsCount[i] <- length(area_yr_obs[[i]])
  }
  
  ## Optional: plot observation distances for relevant area-years
  if(plot){
    
    pdf("Plots/HighP_AreaYears.pdf", width = 6, height = 4.5)
    for(i in 1:nrow(idx.data)){
      
      if(idx.data$ObsCount[i] > 0){
        print(hist(area_yr_obs[[i]], prob = TRUE, breaks = length(area_yr_obs[[i]]),
                   xlab = "Observation distance",
                   main = paste0("Area = ", idx.data$Area[i], " (", idx.data$AreaIdx[i], "), ",
                                 "Year = ", idx.data$Year[i], " (", idx.data$YearIdx[i], ")")))
        print(lines(density(area_yr_obs[[i]]), col = 4, lwd = 2)) 
      }
    }
    dev.off()
  }
  
  ## Return list of affected area-years
  return(idx.data)
  
}