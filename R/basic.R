#' filterPer
#'
#' @param x  a data.frame object, species abundance profile
#' @param row a numberic, 1/2
#' @param percent percent to retain or remove
#' @param include a logistic , T/F
#'
#' @return
#' data.frame
#' @export
#'
#' @examples
filterPer <- function (x, row, percent, include = T)
{
  if (include) {
    index <- apply(x, row, function(x) {
      (sum(x != 0)/length(x)) > percent
    })
  }
  else {
    index <- apply(x, row, function(x) {
      (sum(x != 0)/length(x)) < percent
    })
  }
  if (row == 1) {
    out <- x[index, ]
  }
  else {
    out <- x[, index]
  }
  return(out)
}
