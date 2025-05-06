#' Prepare relevant output data from NIMBLE model for subsequent implementation
#' in the NI database
#'
#' @param mcmc.out an mcmc list containing posterior samples from a model run.
#' @param N_areas integer. Number of areas included in the analyses.
#' @param area_names character vector containing area/location names.
#' @param N_sites matrix containing the number of sites per location/area/county.
#' @param min_years integer vector. Indices of first year with available data
#' for each location/area/county.
#' @param max_years integer vector. Indices of last year with available data for
#' each location/area/county.
#' @param minYear integer. The first year considered in analyses.
#' @param maxYear integer. The last year considered in analyses.
#' @param fitRodentCov logical. If TRUE, makes posterior summaries for rodent 
#' effect in addition to other parameters.
#' @param save logical. If TRUE (default), saves the posterior summaries to an 
#' RDS file.
#'
#' @return a list of data frames with density estimates:
#'   total = average density estimates per area for the whole data period
#'   total.last10yr = average density estimates per area for the last 10 years of data
#'   annual.area = density estimates per area per year
#'   annual.area.adult = density estimates of adult ptarmigan per area per year
#' @export
#'
#' @examples

prepareOutputNI <- function(mcmc.out,
                            N_areas,
                            area_names,
                            N_sites,
                            min_years,
                            max_years,
                            minYear,
                            maxYear,
                            fitRodentCov,
                            save = TRUE) {
  
  ## Convert posterior samples to matrix format
  out.mat <- as.matrix(mcmc.out)
  
  ## Set year index ranges for plotting shared data collection period only
  minYearIdx_shared <- max(min_years)
  maxYearIdx_shared <- min(max_years)

  # Summarize posteriors for population density
  popDens.areayear <- data.frame() # density per area and year
  popDens.ad.areayear <- data.frame() #adult density per area and year
  popDens.tot <- data.frame() #total density for the whole time period
  popDens.10yr <- data.frame() #total density for the last 10 years

  ## Annual densities per area
  for (i in 1:N_areas) {
    # Prepare matrix for temporary storage of results
    popDens.sum <-  data.frame()
    popDens.ad.sum <- data.frame()

    # Determine area-specific year range
    area_yearIdxs <- (1:(maxYear - minYear + 1))
    area_years <- area_yearIdxs + (minYear - 1)

    for (t in 1:length(area_years)) {
      # Summarize annual average population densities
      popDens_juv_t0 <- out.mat[, paste0("meanDens[",  i, ", 1, ", area_yearIdxs[t], "]")]
      popDens_ad_t0 <- out.mat[, paste0("meanDens[",  i, ", 2, ", area_yearIdxs[t], "]")]
      popDens_mean_t0 <- popDens_juv_t0 + popDens_ad_t0

      # Merge average density back into posterior samples (for later calculation of total densities)
      out.mat <- cbind(out.mat, popDens_mean_t0)
      colnames(out.mat)[dim(out.mat)[2]] <- paste0("meanDens[", i, ", ", area_yearIdxs[t], "]")

      # Mean densities per year
      popDens_add <- data.frame(
        Area = area_names[i],
        Year = area_years[t],
        Median = median(popDens_mean_t0),
        lCI = unname(quantile(popDens_mean_t0, probs = 0.25)), #25% lower
        uCI = unname(quantile(popDens_mean_t0, probs = 0.75)), #75% upper
        Mean = mean(popDens_mean_t0),
        SD = sd(popDens_mean_t0)
      )
      popDens.sum <- rbind(popDens.sum, popDens_add)

      # Mean adult densities per year
      popDens_add.ad <- data.frame(
        Area = area_names[i],
        Year = area_years[t],
        Median = median(popDens_ad_t0),
        lCI = unname(quantile(popDens_ad_t0, probs = 0.25)), #25% lower
        uCI = unname(quantile(popDens_ad_t0, probs = 0.75)), #75% upper
        Mean = mean(popDens_ad_t0),
        SD = sd(popDens_ad_t0)
      )
      popDens.ad.sum <- rbind(popDens.ad.sum, popDens_add.ad)
    }

    popDens.areayear <- rbind(popDens.areayear, popDens.sum)
    popDens.ad.areayear <- rbind(popDens.ad.areayear, popDens.ad.sum)
  }

  ## Calculate average area-specific densities over different time periods
  for (i in 1:N_areas) {
    # Entire time period - average total density
    dens_mean_tot <- rowMeans(out.mat[, paste0("meanDens[", i, ", ", area_yearIdxs[1:(length(area_years) - 1)], "]")])
    out.mat <- cbind(out.mat, dens_mean_tot)
    colnames(out.mat)[dim(out.mat)[2]] <- paste0("densAvg_tot[", i, "]")

    # last 10 years - average total density
    dens_mean_10yr <- rowMeans(out.mat[, paste0("meanDens[", i, ", ", area_yearIdxs[(length(area_years) - 10):length(area_years) - 1], "]")])
    out.mat <- cbind(out.mat, dens_mean_10yr)
    colnames(out.mat)[dim(out.mat)[2]] <- paste0("densAvg_10yr[", i, "]")
  }

  ## Summarize average population densities
  for (i in 1:N_areas) {
    Dens_names <-
      c(paste0("densAvg_tot[",  i, "]"),
        paste0("densAvg_10yr[",  i, "]"))
    Dens_tot_add <- data.frame(
      Area = area_names[i],
      SummaryPeriod = paste0(minYear, "-", maxYear),
      Median = median(out.mat[, Dens_names[1]]),
      lCI = unname(quantile(out.mat[, Dens_names[1]], probs = 0.25)), #25% lower
      uCI = unname(quantile(out.mat[, Dens_names[1]], probs = 0.75)), #75% upper
      Mean = mean(out.mat[, Dens_names[1]]),
      SD = sd(out.mat[, Dens_names[1]])
    )
    Dens_10yr_add <- data.frame(
      Area = area_names[i],
      SummaryPeriod = paste0(maxYear - 10, "-", maxYear),
      Median = median(out.mat[, Dens_names[2]]),
      lCI = unname(quantile(out.mat[, Dens_names[2]], probs = 0.25)), #25% lower
      uCI = unname(quantile(out.mat[, Dens_names[2]], probs = 0.75)), #75% upper
      Mean = mean(out.mat[, Dens_names[2]]),
      SD = sd(out.mat[, Dens_names[2]])
    )
    popDens.tot <- rbind(popDens.tot, Dens_tot_add)
    popDens.10yr <- rbind(popDens.10yr, Dens_10yr_add)
  }

  return(
    list(
      total = popDens.tot,
      total.last10yr = popDens.10yr,
      annual.area = popDens.areayear,
      annual.area.adult = popDens.ad.areayear
    )
  )
}