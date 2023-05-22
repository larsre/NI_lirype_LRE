
addDummyDimension <- function(x){
  
  if(length(dim(x)) > 3){
    stop("This function does not currently support adding dummy dimensions to input objects with more than 3 dimensions.")
  }
  
  if(is.vector(x)){
    x.out <- matrix(NA, nrow = 1, ncol = length(x))
    x.out[1,] <- x
  }else{
    x.out <- array(NA, dim = c(1, dim(x)))
    if(length(dim(x.out)) == 3){
      x.out[1,,] <- x
    }else{
      x.out[1,,,] <- x
    }
  }
  
  return(x.out)
}
