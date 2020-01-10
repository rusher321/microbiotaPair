#' Permutational Multivariate Analysis of Variance Using Distance Matrices
#'
#' @description PERMANOVA1
#' The aim of PERMANOVA1 function is to asess the association between the overall profile and phenotype
#' @details 09/01/2020  ShenZhen China
#' @author  Hua Zou
#' @param x x is phenotype cantains the subjects' characteristic;
#' @param y y is the overall profile, colnames are samplenames, rownames are mearsures;
#' @param tag which connect x and y.
#'
#' @usage PERMANOVA1(x, y, tag)
#' @examples  result <- PERMANOVA1(phen, spf, tag="SampleID")
#'
#' @return  Return the F test result
#'
#' @export
#'
PERMANOVA1 <- function(x, y, tag) {

  # install_suggest_packages <- function(packages_name_list = c("randomForest")){
  #   usePackage <- function(p) {
  #     if (!is.element(p, installed.packages()[, 1]))
  #       install.packages(p, dep=TRUE, repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  #     suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
  #   }
  #   invisible(lapply(packages_name_list, usePackage))
  # }
  # packages_list <- c("dplyr", "vegan", "varhandle")
  # install_suggest_packages(packages_list)
  #
  # library(devtools)
  # library(roxygen2)
  # library(dplyr)
  # library(vegan)
  #
  # phen <- read.csv("inst/phen.csv")
  # spf <- read.table("inst/Species.profile")
  # amf <- read.table("inst/Amino.profile")
  # x <- phen
  # y <- spf
  # tag <- "SampleID"
  # use_data(phen, spf, amf)

  colnames(x)[which(colnames(x) == tag) ] <- "SampleID"
  sid <- intersect(as.character(x$SampleID), colnames(y))
  phe <- x %>% filter(SampleID %in% sid)
  prf <-  y %>% select(as.character(phe$SampleID)) %>%
    t() %>% data.frame()
  per <- apply(phe %>% select(-one_of("SampleID")), 2, function(a, pf){
    dat <- data.frame(value = a, pf)
    datphe <- dat$value %>% varhandle::unfactor()
    if (length(datphe) == 0 | unique(datphe) == 1) {
      res <- data.frame(length(datphe), rep(NA, 6))
      next
    }
    if (length(unique(datphe)) < 6) {
      datphe <- as.factor(datphe)
    }
    datprf <- dat[, -1, F]
    dis <- vegan::vegdist(datprf, method = "bray")
    set.seed(123)
    ad <- vegan::adonis(dis ~ datphe, permutations = 1000)
    tmp <- as.data.frame(ad$aov.tab) %>% slice(1)
    res <- c(length(datphe), as.numeric(tmp[, c(1:6)]))
    return(res)
  }, prf) %>% t() %>% data.frame()

  colnames(per) <- c("SumsOfSample", "Df", "SumsOfSqs",
                     "MeanSqs", "F.Model", "R2", "Pr(>F)")
  per$FDR <- p.adjust(per$`Pr(>F)`, method = "BH")
  return(per)
}
