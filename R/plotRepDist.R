#' Plot posterior distributions of annual reproductive rates
#'
#' @param mcmc.out an mcmc list containing posterior samples from a model run.
#' @param N_areas integer. Number of areas included in analyses. 
#' @param area_names character vector containing area/location names.
#' @param N_years integer. Number of years included in analyses. 
#' @param minYear integer. The first year considered in analyses. 
#' @param area character string. Optional name of the area for which to plot. If not provided, plots will be made for all areas. 
#'
#' @return
#' @export
#'
#' @examples

plotRepDist <- function(mcmc.out, N_areas, area_names, N_years, area = NULL){
  
  ## List relevant parameters
  if(is.null(area)){
    idxGrid <- expand.grid(c(1:N_areas), c(1:N_years))
  }else{
    idxGrid <- expand.grid(which(area_names == area), c(1:N_years))
  }
  
  relParams <- paste0("R_year[", idxGrid[,1], ", ", idxGrid[,2], "]")
  
  ## Reduce posterior samples to relevant parameters only
  mcmc.out <- mcmc.out[, which(colnames(mcmc.out[[1]]) %in% relParams), drop = TRUE]
  
  ## Reformat samples
  out.data <- reshape2::melt(as.matrix(mcmc.out))
  colnames(out.data) <- c("Sample", "Parameter", "Value")
  
  ## Add area and year indices
  idx.data <- stringi::stri_extract_all_regex(out.data$Parameter, "[0-9]+", simplify = TRUE)
  colnames(idx.data) <- c("AreaIdx", "YearIdx")
  idx.data <- data.frame(idx.data)
  idx.data$AreaIdx <- as.numeric(idx.data$AreaIdx)
  idx.data$YearIdx <- as.numeric(idx.data$YearIdx)
  idx.data$Year <- idx.data$YearIdx + minYear - 1
  
  out.data <- cbind(out.data, idx.data)
  
  ## Plot overlapping posterior densities
  plot.title <- ifelse(is.null(area), "Annual reproduction, all areas pooled", paste0("Annual reproduction, ", area))
  pdf.name <- ifelse(is.null(area), "AnnualReproduction_Distribution_allAreas.pdf", paste0("AnnualReproduction_Distribution_Area", which(area_names == area), ".pdf"))
  
  pdf(paste0("Plots/", pdf.name), height = 5, width = 7)
  ggplot(out.data, aes(x = Value)) + 
    geom_density(aes(group = Parameter, colour = Year, fill = Year), alpha = 0.1) + 
    geom_density(colour = "deeppink", fill = "deeppink", alpha = 0.5) + 
    scale_color_viridis() + 
    scale_fill_viridis() + 
    ggtitle(plot.title) + 
    theme_bw() + theme(panel.grid = element_blank())
  
  if(is.null(area)){
    for(i in 1:N_areas){
      out.data.sub <- out.data[which(out.data$AreaIdx == i), ]
      print(
        ggplot(out.data.sub, aes(x = Value)) + 
          geom_density(aes(group = Parameter, colour = Year, fill = Year), alpha = 0.1) + 
          geom_density(colour = "deeppink", fill = "deeppink", alpha = 0.5) + 
          scale_color_viridis() + 
          scale_fill_viridis() + 
          ggtitle(paste0("Annual reproduction, ", area_names[i])) + 
          theme_bw() + theme(panel.grid = element_blank())
        
      )
    }
  }
  
  dev.off()
  
  ## Set and return plot path
  plot.path <- paste0("Plots/", pdf.name)
  return(plot.path)
}
