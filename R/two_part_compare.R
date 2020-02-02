#' two part compare
#' From The Gut Microbiome Contributes to a Substantial Proportion of the Variation in Blood Lipids
#' @param dat , microbiome data row is sample id, col is variable
#' @param phe , metadata row is sample id ,col is variable
#' @param response , the variable to explain
#' @param cutoff , ensure to detect or undetect
#' @param number , when sample number is limited , default 10
#'
#' @return dataframe
#' @export
#'
#' @examples
two_part_compare <- function(dat , phe, response, cutoff, number=10){
  # match the sample ID
  id <- intersect(rownames(dat), rownames(phe))
  if(length(id)==0){
    stop("can't match the sample id")
  }
  dat <- dat[id, ]
  y <- phe[id, response]
  # part1 transform the matrix to 0-1
  dat2 <- dat
  dat2[dat2 >= cutoff] <- 1
  dat2[dat2 < cutoff ] <- 0

  # filter all zero variable

  out <- matrix(NA, nrow = ncol(dat), ncol = 3+4+4+2+2)
  for(i in 1:ncol(dat)){
    #print(i)
    x <- dat2[,i]
    out[i,1:3] <- c(sum(x==0), sum(x==1), median(dat[, i]))
    if(sum(x==0) < number |sum(x==1) < number){   # here set the cutoff 20, maybe need the sample size to adjust
      out[i,4:7] <- c(0, 0 , 0, 1)
    }else{
      res <- glm(y~x)
      out[i,4:7] <- as.numeric(summary(res)$coefficients[2,])
    }

  }

  # part2 transform the matrix log
  dat3 <- dat
  dat3[dat3 < cutoff] <- 0
  dat3 <- apply(dat3, 2, log10)
  for(i in 1:ncol(dat)){
    # print(i)
    x <- dat3[,i]
    id2 <- !is.infinite(x)
    x2 <- x[id2]
    y2 <- y[id2]
    if(length(x2) < number){   # here set the cutoff 20, maybe need the sample size to adjust
      out[i,8:11] <- c(0, 0 , 0, 1)
    }else{
      res <- glm(y2~x2)
      out[i,8:11] <- as.numeric(summary(res)$coefficients[2,])
    }

    # meta-analysis
    tmp_meta <- data.frame(beta = out[i, c(4,8)], se = out[i, c(5,9)])
    tmp_meta$se <- ifelse(tmp_meta$se == 0, 0.001, tmp_meta$se) # need to discuss for if se==0, can't get meta result
    res_meta <- rma(yi = beta, data = tmp_meta, sei = se, method = "DL")
    out[i, 12:13] <- c(res_meta$zval, res_meta$pval)

    # association value
    pvalue <- c(out[i, c(7,11,13)])
    estimate <- c(out[i, c(4,8,12)])

    out[i ,15] <- min(pvalue)
    zscore <- qnorm(1-(min(pvalue)/2))
    out[i, 14] <- ifelse(estimate[which.min(pvalue)] >0, zscore, -zscore)

  }

  # meta-analysis


  rownames(out) <- colnames(dat)
  colnames(out) <- c("No.absent", "No.Present", "medianAbundance",
                     paste0("binary", c("_estimate", "_se", "_tvalue", "_p")),
                     paste0("quantitative", c("_estimate", "_se", "_tvalue", "_p")),
                     paste0("Meta", c("_zvalue", "_p")),
                     paste0("Asso", c("_zvalue", "_p")))

  return(out)

}

