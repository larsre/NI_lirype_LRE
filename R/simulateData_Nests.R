## Function for simulating nest survey data
simulateData_Nests <- function(nind.avg.NS, Tmin.NS, Tmax.NS, Tmax, RepRates){

  # Determine number of nests surveyed in each year
  n.nests <- rep(NA, Tmax)
  for(t in Tmin.NS:Tmax.NS){
    n.nests[t] <- rpois(1, nind.avg.NS)
  }

  # Set up empty vectors for storing data
  nChicks.obs <- year.obs <- vector()

  # Simulate nest survey data collection
  for(t in Tmin.NS:Tmax.NS){

    # - Simulate number of chicks in each nest surveyed in year t
    nChicks.obs.t <- rpois(n.nests[t], RepRates[t])

    # - Make a vector with year indices
    year.obs.t <- rep(t, n.nests[t])

    # - Merge data of year t with data from previous years
    nChicks.obs <- c(nChicks.obs, nChicks.obs.t)
    year.obs <- c(year.obs, year.obs.t)
  }

  # Assemble data in a list and return
  NS.data <- list(nChicks.obs = nChicks.obs, year.obs = year.obs)
  return(NS.data)
}