#' @import dplyr
#' @importFrom vegan vegdist adonis
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom varhandle unfactor
#'
#'
#' @title Permutational Multivariate Analysis of Variance Using Distance Matrices
#'
#' @description The PERMANOVA1 aimed to assess the association between the overall profile and phenotype
#' @details 09/01/2020  ShenZhen China
#' @author  Hua Zou
#' @param physeq (Required). A \code{phyloseq} object containing merged information of abundance,
#'        sample data including the measured variables and categorical information of the samples.
#' @param sampleid (Required) A character of the sampleid to connect phenotype and  profiles.
#'
#' @usage PERMANOVA1(physeq, sampleid)
#' @examples
#'
#' data(physeq_data)
#' sampleid <- "SampleID"
#' PERMANOVA1(physeq_data, sampleid)
#'
#' @return  permanova result
#'
#' @export
#'
PERMANOVA1 <- function(physeq, sampleid) {

  # load("data/physeq_data.rda")
  # sampleid <- "sampleID"

  phen <- microbiome::meta(physeq) %>% tibble::rownames_to_column("SampleID")
  prof <- microbiome::abundances(physeq) %>% data.frame()

  colnames(phen)[which(colnames(phen) == sampleid) ] <- "SampleID"
  sid <- intersect(phen$SampleID, colnames(prof))
  phe <- phen[phen$SampleID%in%sid, ]
  prf <-  prof %>% select(phen$SampleID) %>%
    t() %>% data.frame()
  per <- apply(phe %>% select(-one_of("SampleID")), 2, function(x, pf){
    dat <- data.frame(value = x, pf)
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
