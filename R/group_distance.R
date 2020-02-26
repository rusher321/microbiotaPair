#' plot_pair_distance
#'
#' @param microbiota a dataframe object,species abundance, row is sample Id
#' @param metadata  a dataframe object, metadata, row is sample Id
#' @param method  a  character, distance method eg "bray"
#' @param percent.outlier a number, percent of the ourlier
#'
#' @return
#' ggplot object
#' @export
#'
#' @examples
plot_pair_distance <- function(microbiota, metadata, method = "bray", percent.outlier = 0.05){

  # match the ID , get the baseline & treatment data
  matchname <- names(table(metadata[, pairID_varname]))[table(metadata[, pairID_varname]) == 2]
  outconfig <- metadata[!is.na(match(metadata[, pairID_varname], matchname)), ]
  matchdat <-  metadata[order(metadata[, pairID_varname]), ]
  matchdat <- matchdat[order(matchdat[, time_varname]), ]
  number <- length(matchname)
  # to make sure the microbiota's sample ID is row
  matchmicrobiota <- microbiota[rownames(matchdat), ]
  if(nrow(matchmicrobiota)==0){
    matchmicrobiota <- t(microbiota)[rownames(matchdat), ]
  }
  # to get the distance
  distanceM <- as.matrix(vegan::vegdist(matchmicrobiota, method = method))

  pairdistance <- diag(distanceM[(number+1):(2*number), 1:number])

  intradistanceg1 <- distanceM[1:number, 1:number]
  intradistanceg1 <- intradistanceg1[lower.tri(intradistanceg1)]
  intradistanceg2 <- distanceM[(number+1):(2*number), (number+1):(2*number)]
  intradistanceg2 <- intradistanceg2[lower.tri(intradistanceg2)]

  exdistance <- distanceM[(number+1):(2*number), 1:number]
  exdistance <- c(exdistance[lower.tri(exdistance)], exdistance[upper.tri(exdistance)])

  # fiter the outlier
  filterE <- function(x){
    num <- length(x)
    y <- sort(x)
    y2 <- y[round(num*percent.outlier):round(num*(1-percent.outlier))]
    return(y2)
  }

  pairdistancef <- filterE(pairdistance)
  intradistanceg1f <- filterE(intradistanceg1)
  intradistanceg2f <- filterE(intradistanceg2)
  exdistancef <- filterE(exdistance)

  # plot
  gp <- levels(metadata[, time_varname])
  dist <- c(pairdistancef, intradistanceg1f, intradistanceg2f, exdistancef)
  group <- c(rep("pair", length(pairdistancef)), rep(gp[1], length(intradistanceg1f)),
             rep(gp[2], length(intradistanceg2f)), rep("nopair", length(exdistancef)))
  qdat <- data.frame(value = dist, type = group)
  p  <- ggplot(qdat, aes(value, colour=type))
  p  <- p + geom_density(size=1.5)+scale_color_grey()+mytheme

  return(p)
}
