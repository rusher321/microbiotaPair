#' pair_foldchange
#' The fold change value between the pair samples
#' @param microbiota a dataframe object,species abundance, row is sample Id
#' @param metadata  a dataframe object, metadata, row is sample Id
#' @param percent  numberic, to cutoff value of species occranece
#' @param order    logistic , to detemain the time order
#'
#' @return
#' list include ggplot object and the data.frame of foldchange
#' @export
#'
#' @examples
pair_foldchange <- function(microbiota, metadata, percent, order){

  # match the ID , get the baseline & treatment data
  matchname <- names(table(metadata[, pairID_varname]))[table(metadata[, pairID_varname]) == 2]
  outconfig <- metadata[!is.na(match(metadata[, pairID_varname], matchname)), ]
  matchdat <-  metadata[order(metadata[, pairID_varname]), ]
  if(order==F){
    matchdat <- matchdat[order(matchdat[, time_varname]), ]
  }else{
    matchdat <- matchdat[order(matchdat[, time_varname], decreasing = T), ]
  }
  number <- length(matchname)
  # to make sure the microbiota's sample ID is row
  matchmicrobiota <- microbiota[rownames(matchdat), ]

  # to remove low occrance feature

  matchmicrobiota <- filterPer(matchmicrobiota, row = 2, percent = percent)
  g1microbiota <- matchmicrobiota[1:number, ]
  g2microbiota <- matchmicrobiota[(number+1):(2*number), ]

  # pair fold change
  pairfd <- function(x, y){
     fd <- rep(NA, length(x))
     index <- which(x>0 & y>0)
     value <- y[index]/x[index]
     fd[index] <- value
     # filter the extreme vlaue
     value <- filterEx(value, percent.outlier = 0.05)
     if(x==0 & y>0){
        fd[x==0 & y>0] <- max(value)
     }else if(x>0 & y==0){
        fd[x>0 & y==0] <- min(value)
     }else if(x==0 & y==0){
        fd[x==0 & y==0] <- 1
     }
     return(fd)
  }

  # to get the pair fold change
  tmp <- c()
  x <- c()
  y <- c()
  for(i in 1:ncol(g1microbiota)){
    tmp[i] <- pairfd(g1microbiota[,i], g2microbiota[,i])
    x[i] <- mean(filterEx(g1microbiota[,i], percent.outlier = 0.05))
    y[i] <- mean(filterEx(g2microbiota[,i], percent.outlier = 0.05))
  }


  qdat <- data.frame(tax = colnames(g1microbiota), FoldChange = tmp,
                     before = x, treat = y)
  qdat$FoldChange2 <- abs(log(qdat$FoldChange))
  qdat$enrich <- ifelse(qdat$FoldChange>2, time_name[2],
                        ifelse(qdat$FoldChange<0.5, time_name[1], "None"))
  qdat$enrich <- factor(qdat$enrich, levels = c(time_name, "None"))
  p <- ggplot(qdat, aes(x=log(before), y=log(treat), size=FoldChange2, fill=enrich)) +
    geom_point(alpha=0.5, shape=21, color="black") +
    theme(legend.position="bottom") +
    scale_fill_manual(values = c(time_colour, "black"))+
    xlab("Baseline (Log)") +
    ylab("Treatment (Log)") +
      mytheme
  out <- list(qdat, p)
  return(out)

}
