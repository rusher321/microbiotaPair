#' @import vegan
#' @import phyloseq
#' @importFrom tsnemicrobiota tsne_phyloseq
#'
#'
#' @title Ordination and beta dispersion
#'
#' @description This function performs ordination using a selected method. It also performs PERMANOVA using a distance matrix corresponding
#' to the taxa abundance under the different conditions/groups. It also computes  pairwise beta dispersion for
#' all conditions of a selected grouping variable, beta dispersion is measured as the average distance of group members to the group centroid.
#' The function a solution of ordination, PERMANOVA results and beta dispersion results. See value for details.
#'
#' @details 14/01/2020  ShenZhen China
#' @author  Hua Zou
#'
#' @param physeq (Required). A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param method (Optional). A character string specifying ordination method. All methods available to the \code{ordinate} function
#'               of \code{phyloseq} are acceptable here as well ( Default is "PCOA").
#'               "DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA", "PCA", "Tsne".
#' @param which_distance (Optional). A string character specifying dissimilarity index to be used in calculating pairwise distances (Default index is "bray".).
#'                       "unifrac","wunifrac","manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower",
#'                       "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao" or "mahalanobis".
#' @param grouping_column (Required). Character string specifying name of a categorical variable that is preffered for grouping the information.
#'        information.
#'
#' @return Returns a list of three items: \itemize{
#'         \item A solution of the ordination
#'         \item A \code{data.frame} of betadispersion results; compared groups, corresponding pairwise observed p-values
#'         and significant labels.
#'         \item an object of class "adonis" with all components ; see \code{adonis} for details.
#'         }
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#'
#' @usage ordination(physeq, method, grouping_column)
#' @examples
#' data(physeq_data)
#' physeq <- physeq_data
#' ord.res <- ordination(physeq, method="NMDS", grouping_column="Stage")
#'
#'
#' @export ordination
#'
ordination <- function(physeq, which_distance="bray", method="NMDS", grouping_column, pvalue.cutoff=0.05){

  meta_table <- data.frame(phyloseq::sample_data(physeq))
  otu_table <- data.frame(phyloseq::otu_table(physeq))
  meta_table$Groups <- meta_table[, grouping_column]

  beta_disper <- function(physeq,grouping_column,pvalue.cutoff=0.05,which_distance){

    meta_table <- data.frame(sample_data(physeq))
    meta_table$Groups <- meta_table[,grouping_column]
    # compute beta dispersion
    mod<- vegan::betadisper(phyloseq::distance(physeq,method=which_distance),meta_table$Groups,type="centroid")
    # compute pairwise beta dispersion for all levels in the grouping variable
    pmod <- vegan::permutest(mod, permutations = 99, pairwise = TRUE)
    p.values <- pmod$pairwise$observed
    # extract significantly dispersed pairs and assign significant labels
    p.values <- p.values[!is.na(p.values <= pvalue.cutoff)]
    signi_label <- paste(cut(p.values,breaks=c(-Inf,0.001,0.01,0.05, Inf), label=c("***", "**", "*", ".")))
    groups_compared <- names(p.values)

    betadisper_res <- data.frame(groups_compared, p.values, signi_label)
    out <- list("betadisper_res"=betadisper_res, "pmod"=pmod)
    return(out)
  }


  sol <- NULL
  if(method == "PCoA"){
    sol <- cmdscale(phyloseq::distance(physeq, which_distance),eig=T)
  }else if(method == "PCA"){
    sol <- stats::prcomp(otu_table, scale. = TRUE, rank. = 2)
  }else if(method == "Tsne"){
    sol <- tsnemicrobiota::tsne_phyloseq(physeq, distance=which_distance,
                                          perplexity = 8, verbose=0, rng_seed = 3901)
  }else{
    sol <- phyloseq::ordinate(physeq, method, distance=which_distance)
  }
  dist <- phyloseq::distance(physeq,which_distance)
  adonis_res <- vegan::adonis(dist ~ Groups, data=meta_table)

  betadisper_pw <- beta_disper(physeq, grouping_column, pvalue.cutoff, which_distance)
  betadisper_pw <- betadisper_pw$betadisper_res

  out <- list("solution"=sol,"betadispersion"=betadisper_pw,"adonis_res"=adonis_res, "groups"=meta_table$Groups)
  return(out)
}
