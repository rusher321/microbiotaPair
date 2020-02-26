#' two part compare
#' two part (paired test & McNemar’s Chi-Square Test) to maximum the data infomation
#' @param microbiota  a data.frame, microbiome data row is sample id, col is variable
#' @param metadata  a data.frame, metadata row is sample id ,col is variable
#' @param percent  numberic, 0-1, the cutoff of species occrance
#' @return
#' data.frame
#' @export
#'
#' @examples
two_part_compare <- function(microbiota , metadata, percent){
  # match the ID , get the baseline & treatment data
  matchname <- names(table(metadata[, pairID_varname]))[table(metadata[, pairID_varname]) == 2]
  outconfig <- metadata[!is.na(match(metadata[, pairID_varname], matchname)), ]
  matchdat <-  metadata[order(metadata[, pairID_varname]), ]
  # to order the varname ,keep same as the time_name
  matchdat[, time_varname] <- factor(matchdat[, time_varname], levels = time_name)
  matchdat <- matchdat[order(matchdat[, time_varname]), ]

  number <- length(matchname)
  # to make sure the microbiota's sample ID is row
  matchmicrobiota <- microbiota[rownames(matchdat), ]

  # to remove low occrance feature
  matchmicrobiota <- filterPer(matchmicrobiota, row = 2, percent = percent)
  g1microbiota <- matchmicrobiota[1:number, ]
  g2microbiota <- matchmicrobiota[(number+1):(2*number), ]

   # part1 transform the matrix to 0-1
  g1microbiota2 <- g1microbiota
  g2microbiota2 <- g2microbiota
  g1microbiota2[g1microbiota2 >0] <- 1
  g1microbiota2[g1microbiota2 == 0] <- 0
  g2microbiota2[g2microbiota2 >0] <- 1
  g2microbiota2[g2microbiota2 == 0] <- 0

  # filter all zero variable

  out <- matrix(NA, nrow = ncol(matchmicrobiota), ncol = 5+2+3)
  for(i in 1:ncol(matchmicrobiota)){
    x <- g1microbiota2[,i]
    y <- g2microbiota2[,i]
    out[i, 1:5] <- c(length(x), sum(x==0), sum(y==0), median(g1microbiota[,i]),
                     median(g2microbiota[,i]))
    x <- factor(x, levels = c(0, 1))
    y <- factor(y, levels = c(0, 1))
    tabletmp <- table(x, y)
    mcnemarres <- mcnemar.test(tabletmp)
    out[i, 6:7] <- as.numeric(c(mcnemarres$statistic, mcnemarres$p.value))
  }

  # part2 transform the matrix log
  minvalue <- min(matchmicrobiota[matchmicrobiota!=0])
  g1microbiota3 <- log(g1microbiota+minvalue)
  g2microbiota3 <- log(g2microbiota+minvalue)
  for(i in 1:ncol(matchmicrobiota)){
    x <- g1microbiota3[,i]
    y <- g2microbiota3[,i]
    if(x!=0)
    ttestres <- t.test(x, y ,paired = T)
    out[i, 8:10] <- as.numeric(c(ttestres$statistic, ttestres$estimate, ttestres$p.value))

  }

  # meta-analysis
  rownames(out) <- colnames(matchmicrobiota)
  colnames(out) <- c("No.sample",
                     paste0(rep(c("No.absent.", "medianAbundance."),each=2), time_name),
                     "Chisq_Mcnemar", "Pvalue_Mcnemar",
                     "T_Ttest", "Estimate_Ttest", "Pvalue_Ttest"
                     )

  return(out)

}

