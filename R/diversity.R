#' @import dplyr
#' @import phyloseq
#' @import microbiome
#' @import ggplot2
#' @importFrom tibble column_to_rownames rownames_to_column
#'
#'
#' @title microbial diversity show the ecosystem state
#'
#' @description diversity
#' The aim of diversity function is to calculate the global indicators of the ecosystem state
#' @details 10/01/2020  ShenZhen China
#' @author  Hua Zou
#' @param x x is phenotype cantains the subjects' characteristic;
#' @param y y is the overall profile, colnames are samplenames, rownames are mearsures;
#' @param tag which connect x and y.
#' @param group beta diversity in group
#'
#' @usage diversity(x, y, tag, group)
#' @examples  result <- diversity(phen, spf, "SampleID", "Group")
#'
#' @return  Return diversity result
#'
#' @export
#'
diversity <- function(x, y, tag, group) {

  # install_suggest_packages <- function(packages_name_list = c("randomForest")){
  #   usePackage <- function(p) {
  #     if (!is.element(p, installed.packages()[, 1]))
  #       BiocManager::install(p, site_repository="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
  #     suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
  #   }
  #   invisible(lapply(packages_name_list, usePackage))
  # }
  # packages_list <- c("microbiome")
  # install_suggest_packages(packages_list)
  #
  # library(dplyr)
  #
  # phen <- read.csv("inst/phen.csv")
  # spf <- read.table("inst/Species.profile")
  # x <- phen
  # y <- spf
  # tag <- "SampleID"
  # group <- "Group"
  # usethis::use_data(phen, spf, amf)

  colnames(x)[which(colnames(x) == tag) ] <- "SampleID"
  colnames(x)[which(colnames(x) == group) ] <- "group"
  sid <- intersect(as.character(x$SampleID), colnames(y))
  phe <- x %>% filter(SampleID %in% sid) %>% tibble::column_to_rownames("SampleID")
  prf <-  y %>% select(as.character(rownames(phe))) %>% as.matrix()

  OTU <- phyloseq::otu_table(prf, taxa_are_rows = TRUE)
  PHE <- phyloseq::sample_data(phe)
  physeq <- phyloseq::phyloseq(OTU, PHE) # the rownames of two data must be same
  # metagMisc::phyloseq_to_df: Convert phyloseq object to data frame (for exporting)

  # alpha
  alpha_result <- microbiome::alpha(physeq, index = c("observed", "chao1", "diversity_inverse_simpson",
                  "diversity_shannon","diversity_coverage","evenness_camargo", "evenness_pielou",
                  "evenness_simpson","evenness_evar","evenness_bulla","dominance_dbp","dominance_dmn",
                  "dominance_absolute","dominance_relative", "dominance_simpson",
                  "dominance_core_abundance", "rarity_low_abundance","rarity_rare_abundance"), zeroes = TRUE)

  # diversity analysis from https://microbiome.github.io/tutorials/Betadiversity.html

  # Intra-individual divergence
  dat <- meta(physeq) %>% mutate(sample=rownames(meta(physeq)))
  betas_Intra <- list()
  groups <- as.character(unique(dat$group))
  for (g in groups) {

    df <- subset(dat, group == g)
    beta <- c()

    for (subj in df$ID) {
      # Pick the samples for this subject
      dfs <- subset(df, ID == subj)
      # Check that the subject has two time points
      if (nrow(dfs) == 2) {
        s <- as.character(dfs$sample)
        # Here with just two samples we can calculate the
        # beta diversity directly
        beta[[subj]] <- microbiome::divergence(
                        microbiome::abundances(physeq)[, s[[1]]],
                        microbiome::abundances(physeq)[, s[[2]]],
                                   method = "bray")
      }
    }
    betas_Intra[[g]] <- beta
  }

  return(list(alpha=alpha_result, beta=betas_Intra))

}


#' @import dplyr
#' @import ggplot2
#'
#'
#' @title The visualization of alpha diversity
#'
#' @description alpha_visualization
#' The aim of alpha_visualization function is to display the alpha diversity by boxplot
#' @details 10/01/2020  ShenZhen China
#' @author  Hua Zou
#' @param x x is phenotype cantains the subjects' characteristic;
#' @param y y is the result of diversity by diversity function;
#' @param index the type of diversity measures
#'
#' @usage alpha_visualization(x, y, index
#' @examples  result <- alpha_visualization(phen, diversity_result$alpha, "diversity_shannon")
#'
#' @return  Return boxplot
#'
#' @export
#'
alpha_visualization <- function(x,
                                y=diversity_result$alpha,
                                index="diversity_shannon"){

  y$SampleID <- rownames(y)
  mdat <- inner_join(x, y, by = "SampleID")
  colnames(mdat)[which(colnames(mdat) == index)] <- "Index"

  return(
    ggplot(mdat %>% mutate(Stage=factor(Stage, levels = c("Before", "After"))),
           aes(x = Stage, y = Index))+
      stat_boxplot(aes(color = Stage),
                   geom = "errorbar",
                   width = 0.15,
                   size = .3)+
      geom_boxplot(aes(fill = Stage),
                   width = .4,
                   outlier.shape = 1,
                   size = .3)+
      guides(fill = F, color = F)+
      facet_wrap(facets = "Group")+
      labs(x = "", y = "Shannon Index")+
      theme_classic()+
      theme(axis.title = element_text(size=10, color="black", face="bold"),
            axis.text = element_text(size=9, color="black"),
            text = element_text(size=8, color="black"),
            strip.text = element_text(size=9, color="black", face="bold"),
            panel.grid = element_blank(),
            legend.text=element_text(size=10, color = "black"),
            legend.position = c(1, 0),
            legend.justification = c(1, 0),
            legend.background = element_rect(color = "black", fill="white"))
  )
}
