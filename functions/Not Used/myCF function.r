myCF <- function(x){
  x <- drop(as.matrix(x))
  if (length(dim(x)) > 1){ # if it's a matrix
    res <- (rowSums(x > 0)-1)/(rowSums(x)-1)
  } else {                 # if it's a vector
    res <- (sum(x > 0)-1)/(sum(x)-1)
  }
  return(res)
}