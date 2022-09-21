
setupModel <- function(modelCode.path, 
                       nim.data, nim.constants,
                       niter = 25000, nthin = 5, nburn = 5000, nchains = 3,
                       testRun = FALSE, initVals.seed){

  ## Load model code
  require('nimble')
  source(modelCode.path)
  
  ## Set parameters to monitor
  params <- c("esw", "R_year", "p", "D",
              "Density", "N_exp", 
              "mu.D1", "sigma.D", "mu.R", "sigma.R",
              "Mu.S1", "Mu.S2", "ratio.JA1")
  
  ## Simulate initial values
  set.seed(initVals.seed)
  initVals <- list()
  for(c in 1:nchains){
    initVals[[c]] <- simulateInits(nim.data = nim.data, nim.constants = nim.constants)
  }
  
  ## Adjust MCMC parameters if doing a test run
  if(testRun){
    niter <- 2
    nthin <- 1
    nburn <- 0
  }
  
  ## Collate model setup variables in a list
  setup <- list(
    modelCode = rypeIDSM,
    modelParams = params,
    initVals = initVals,
    mcmcParams = list(niter = niter, nthin = nthin, 
                      nburn = nburn, nchains = nchains)
  )
   
  return(setup) 
}