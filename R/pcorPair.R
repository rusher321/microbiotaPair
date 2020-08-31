#' @import ppcor
#' @importFrom phyloseq otu_table sample_data
#'
#'
#' @title pairwise partial correlations for Pair Study
#'
#' @description To estimate the relationship bwtween features/taxa and phenotypes
#' given others.
#'
#' @details 15/01/2020  ShenZhen China
#' @author  Huahui Ren
#'
#' @param microbiota (Required). data.frame, microbiota ,row is sample
#' @param metadata (Required). data.frame. phenotype, row is sample
#' @param metadataVar (Required). a set of characters
#' @param confounder (Required). A set of Characters string that are covariates.
#'        information.
#' @param method (Required). A chareacter string indicating which partial correlation
#'         coefficient is to be computed. "peason","kendall" or "spearman" can be abbreviated.
#' @return Returns a list with the results of pcc.
#'
#' @usage PcorPair(microbiota, metadata, metadataVar, confounder,method = "s")
#'
#' @example
#'
#' library(dplyr)
#' library(geepack)
#' data("physeq_data")
#' physeq <- physeq_data
#' microbitota <- otu_table(phyloseq)
#' metadata <- sample_table(phyloseq)
#' metadataVar <- c("BMI", "Glycine")
#' confounder <- c("Age", "Sex")
#'
#' PcorPair(microbiota, metadata, metadataVar, confounder, method = "s")
#'
#' @export
#'

PcorPair <- function (microbiota, metadata, metadataVar, confounder,
                      method = "s", time_varname){

  # match the sample ID
  id <- intersect(rownames(microbiota), rownames(metadata))
  dataset <- microbiota[id, ]
  metadata <- metadata[id, c(metadataVar, confounder, time_varname), drop=F]
  confounderindex <- which(colnames(metadata) %in% confounder)
  time_varnameindex <- which(colnames(metadata) %in% time_varname)

  if (length(confounderindex) != length(confounder)) {
    stop("please check the confounder variable")
  }
  # ready the data
  datacon <- metadata[, confounderindex]
  metadatafilter <- metadata[, -c(confounderindex, time_varnameindex)]

  # split the baseline & treatment based on your timevar
  id <- list()
  dat <- list()
  meta <- list()
  micro <- list()
  time_name <- metadata[,time_varnameindex]
  for(i in 1:length(time_name)){
    id[[i]] <- rownames(metadata[metadata[,time_varname]==time_name[i], ])
    dat[[i]] <- datacon[id[[i]], ]
    meta[[i]] <- metadatafilter[id[[i]], ]
    micro[[i]] <- dataset[id[[i]], ]
  }
  # output
  allresult <- list()
  result <- matrix(NA, nrow = ncol(metadatafilter), ncol = ncol(dataset) * 2)
  result <- as.data.frame(result)
  rownames(result) <- colnames(metadatafilter)

  # analysis
  for(m in 1:length(time_name)){
    for (i in 1:c(ncol(metadatafilter))) {
      for (j in 1:ncol(dataset)) {
        dat_com <- data.frame(x = meta[[m]][, i], y = micro[[m]][, j], dat[[m]])
        dat_com <- dat_com[!apply(dat_com, 1, function(x) {any(is.na(x))}), ]

      # if the variable is constant
        op1 <- sum(dat_com$x == dat_com$x[1]) == length(dat_com$x)
        op2 <- sum(dat_com$y == dat_com$y[1]) == length(dat_com$y)
        if (op1 | op2) {
          result[i, c((2 * j - 1):(2 * j))] <- c(0, 1)
        }else {
          tmp <- pcor.test(dat_com$x, dat_com$y, dat_com[,confounder],
                         method = method)
          result[i, c((2 * j - 1):(2 * j))] <- c(tmp$estimate,
                                               tmp$p.value)
      }
      colnames(result)[c((2 * j - 1):(2 * j))] <- paste0(colnames(dataset)[j],
                                                         c("_estimate", "_p.value"))
    }
    }
    allresult[[m]] <- result
  }
  names(allresult) <- time_name
  return(allresult)
}
