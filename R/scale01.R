#' Scale a numeric vector to between 0 & 1
#'
#' This can handle vectors full of zeros
#'
#' @param x Numeric vector.
#'
#' @return Numeric vector, scaled to 0,1
#'
#' @examples
#' x <- c(0,1,1,2,5,6,7,3)
#' scale01(x)
#'
#' @export
#'

scale01 <- function(x){
  if(max(x,na.rm=TRUE) == 0){return(x)} # if all zeros, do nothing
  if(max(x,na.rm=TRUE) > 0){return((x - min(x,na.rm=TRUE)) / (max(x,na.rm=TRUE) - min(x,na.rm=TRUE)))} # if some positive, scale as normal
  if(sum(x,na.rm=TRUE) < 0){x <- abs(x)
  scaled <- (x - min(x,na.rm=TRUE)) / (max(x,na.rm=TRUE) - min(x,na.rm=TRUE))
  return(-scaled)} # if all negative, use absolute values, then replace negative sign
}
