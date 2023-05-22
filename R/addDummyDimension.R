#' Add dummy dimensions to conform with dimensionality requirements of multi-area model
#'
#' @param x a vector, matrix, or 3-dimensional array that is lacking a dummy 
#' area dimension.
#'
#' @return
#' @export a matrix, 3- or 4-dimensional array constituting `x` with the added
#' dummy dimension. 
#'
#' @examples

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
