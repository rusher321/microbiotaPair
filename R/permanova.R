#' A brief introduction to this function.
#'
#' @description PERMANOVA1
#' The aim of PERMANOVA1 function is to asess the association between the overall profile and phenotype
#' @details 09/01/2020  ShenZhen China
#' @author  Hua Zou
#' @param x x is phenotype cantains the subjects' characteristic;
#' @param y y is the overall profile
#'
#' @usage PERMANOVA1(x, y)
#' @examples  result <- PERMANOVA1(phen, spf)
#'
#' @return  Return the F test result
#' @export
#'

PERMANOVA1 <- function(x, y) {

  # install_suggest_packages <- function(packages_name_list = c("randomForest")){
  #   usePackage <- function(p) {
  #     if (!is.element(p, installed.packages()[, 1]))
  #       install.packages(p, dep=TRUE, repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  #     suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
  #   }
  #   invisible(lapply(packages_name_list, usePackage))
  # }
  #
  # install_suggest_packages(c("dplyr", "vegan"))
  #
  # library(devtools)
  # library(roxygen2)
  # library(dplyr)
  #
  # phen <- read.csv("data/phen.csv")
  # spf <- read.table("data/Species.profile")
  # amf <- read.table("data/Amino.profile")
  # use_data(phen, spf, amf)

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
    dis <- vegdist(datprf, method = "bray")
    set.seed(123)
    ad <- adonis(dis ~ datphe, permutations = 1000)
    tmp <- as.data.frame(ad$aov.tab) %>% slice(1)
    res <- c(length(datphe), as.numeric(tmp[, c(1:6)]))
    return(res)
  }, prf) %>% t() %>% data.frame()

  colnames(per) <- c("SumsOfSample", "Df", "SumsOfSqs",
                     "MeanSqs", "F.Model", "R2", "Pr(>F)")
  per$FDR <- p.adjust(per$`Pr(>F)`, method = "BH")
  return(per)
}
