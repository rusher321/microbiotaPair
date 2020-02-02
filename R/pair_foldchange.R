pair_foldchange <- function(microbiota, metadata, percent){

  # match the ID , get the baseline & treatment data
  matchname <- names(table(metadata[, pairID_varname]))[table(metadata[, pairID_varname]) == 2]
  outconfig <- metadata[!is.na(match(metadata[, pairID_varname], matchname)), ]
  matchdat <-  metadata[order(metadata[, pairID_varname]), ]
  matchdat <- matchdat[order(matchdat[, time_varname]), ]
  number <- length(matchname)
  # to make sure the microbiota's sample ID is row
  matchmicrobiota <- microbiota[rownames(matchdat), ]

  # to remove low occrance feature

  matchmicrobiota <- filterPer(matchmicrobiota, row = 2, percent = 0.2)
  g1microbiota <- matchmicrobiota[1:number, ]
  g2microbiota <- matchmicrobiota[(number+1):(2*number), ]

  # pair fold change
  pairfd <- function(x, y){
     index <- which(x>0 & y>0)
  }


}
