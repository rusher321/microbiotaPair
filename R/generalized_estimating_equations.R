#' @import geepack
#' @import dplyr
#' @importFrom phyloseq otu_table sample_data
#'
#'
#' @title Generalized Estimating Equations
#'
#' @description To estimate the relationship bwtween features/taxa and phenotypes.
#'
#' @details 15/01/2020  ShenZhen China
#' @author  Hua Zou
#'
#' @param physeq (Required).  A \code{phyloseq} object containing merged information of abundance,
#'        sample data including the measured variables and categorical information of the samples.
#' @param variable (Required). A character string as the measurement (repeated responses).
#' @param adjust_covariates (optional). A set of Characters string that are covariates.
#'        information.
#' @param PID (Required) the repeated subjects
#'
#' @return Returns a data.frame with the results of gee.
#'
#' @usage gee_test(physeq, variable, adjust_covariates, PID, pcutoff)
#'
#' @examples
#'
#' library(dplyr)
#' library(geepack)
#' data("physeq_data")
#' physeq <- physeq_data
#' variable <- c("BMI", "Glycine")
#' adjust_covariates <- c("Age", "Sex")
#' PID <- "ID"
#' pcutoff <- 0.05
#'
#' gee_test(physeq, variable, adjust_covariates, PID, pcutoff)
#'
#' @export plot_ordination
#'



gee_test <- function(physeq, variable="BMI", adjust_covariates=c("Age", "Sex"), PID="ID", pcutoff=0.05){

  phen <- data.frame(phyloseq::sample_data(physeq))
  prof <- data.frame(phyloseq::otu_table(physeq))

  sid <- intersect(rownames(phen), colnames(prof))
  phe <- phen[order(rownames(phen)%in%sid), ]
  prf <- prof[, order(colnames(prof)%in%sid)]


  gee_res <- function(measure, pvalue=0.05){
    phs <- phe[, c(measure, adjust_covariates, PID)]
    colnames(phs)[which(colnames(phs) == PID)] <- "id"

    res <- apply(prf, 1, function(x, datphe){

      df <- cbind(species=as.numeric(x), datphe)
      df$species <- scale(df$species, scale=T, center=T)
      fm <- formula(paste(measure,
                          paste("species",
                                paste(adjust_covariates, collapse = " + "),
                                sep = " + "),
                          sep = " ~ "))

      fit <- geeglm(fm, data=df, family=gaussian, id=id, corstr="exchangeable")
      res.fit <- data.frame(coef(summary(fit)))[2, ]

      res <- cbind(res.fit[,1], res.fit[,2],res.fit[,3], res.fit[,4], measure)
      return(res)
    }, phs) %>% t(.) %>% data.frame()

    colnames(res) <- c("Estimate", "Std.err", "Wald", "Pr(>|W|)", "Variable")
    res$FDR <- p.adjust(as.numeric(as.character(res[, 4])), method = "BH")

    if(pvalue > 0.05){
      res <- res
    }else{
      res <- res[which(res$FDR < pvalue), ]
    }
    return(res)
  }

  res <- NULL
  for(i in 1:length(variable)){
    tmp <- gee_res(variable[i], pcutoff)
    if(is.null(res)){
      res <- tmp
    }else{
      res <- rbind(res, tmp)
    }

    return(res)
  }
}

