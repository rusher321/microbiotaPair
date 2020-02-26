#' mantel_test
#'
#' @param microbiota  a data.frame object, species profile, row is sample
#' @param metadata  a data.frame object, metadata inf
#' @param methodf  feature space distance method
#' @param methods sample space distance method, vegdist
#'
#' @return
#' figure
#' @export
#'
#' @examples
mantel_test <- function(microbiota,  metadata, methodf, methods){


  # match the ID , get the baseline & treatment data
  matchname <- names(table(metadata[, pairID_varname]))[table(metadata[, pairID_varname]) == 2]
  outconfig <- metadata[!is.na(match(metadata[, pairID_varname], matchname)), ]
  matchdat <-  metadata[order(metadata[, pairID_varname]), ]
  matchdat <- matchdat[order(matchdat[, time_varname]), ]
  number <- length(matchname)
  # to make sure the microbiota's sample ID is row
  matchmicrobiota <- microbiota[rownames(matchdat), ]
  matchmicrobiota <- filterPer(matchmicrobiota, row = 2, percent = 0.2)
  g1microbiota <- matchmicrobiota[1:number, ]
  g2microbiota <- matchmicrobiota[(number+1):(2*number), ]

  # to compute the feature space
  g1distf <- vegan::vegdist(g1microbiota, method = methodf)
  g2distf <- vegan::vegdist(g2microbiota, method = methodf)


  set.seed(000)

  par(mfrow = c(1,2))

  mantelf <- mantel(g1distf, g2distf, method = "spearman")
  sim <- mantelf$perm
  obs <- mantelf$statistic
  mr1 <- ade4::as.randtest(sim = sim, obs = obs)
  plot(mr1, main = paste0("Feature_Space p.value = ", round(mr1$pvalue, dig = 5)))

  # to compute the sample space
  g1dists <- vegan::vegdist(t(g1microbiota), method = methods)
  g2dists <- vegan::vegdist(t(g2microbiota), method = methods)
  mantels <- mantel(g1dists, g2dists, method = "spearman")
  sim <- mantels$perm
  obs <- mantels$statistic
  mr2 <- ade4::as.randtest(sim = sim, obs = obs)
  plot(mr2, main = paste0("Sample_Space p.value = ", round(mr2$pvalue, dig = 5)))

  # libray(customLayout)

}
