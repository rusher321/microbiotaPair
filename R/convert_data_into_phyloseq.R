#' @importFrom data.table fread
#' @importFrom dplyr %>% inner_join
#' @importFrom phyloseq otu_table sample_data phyloseq
#' @importFrom tibble column_to_rownames rownames_to_column
#'
#' @title data2phyloseq
#' @description Convert all the variate such as phenotype metadata and the relative abundance profiles in the data.frame/matrix into the data of phyloseq type
#' @details 13/01/2020  ShenZhen China
#' @author  Hua Zou
#' @param x x  (Required) A table including multiple vectors of subjects' characteristic;
#' @param y y  (Required) A table including multiple metabolic measures of subjects;
#' @param z z  (Required) A table including multiple taxonomy, colnames are sampleid, rownames are taxonomy;
#' @param sampleid A  (Required) character of the subjects' ID to connect the three groups.
#'
#' @usage data2phyloseq(x, y, z, sampleid)
#' @examples
#' library(phyloseq)
#' library(dplyr)
#' library(data.table)
#' library(tibble)
#'
#'
#' x <- read.csv(system.file("extdata", "phenotype.csv", package="microbiotaPair"))
#' y <- fread(system.file("extdata", "metabolic.profile", package="microbiotaPair"))
#' z <- fread(system.file("extdata", "abundance.profile", package="microbiotaPair"))
#' sampleid <- "SampleID"
#'
#' physeq_data <- data2phyloseq(x, y, z, sampleid)
#' print(physeq_data)
#'
#' @return A phyloseq object
#'
#' @export
#'

data2phyloseq <- function(x, y, z, sampleid="SampleID"){

  # x=read.csv("temp/phen.csv")
  # y=fread("temp/Amino.profile")
  # z=fread("temp/Species.profile")
  # sampleid="SampleID"

  colnames(x)[which(colnames(x) == sampleid)] <- "Sample"
  meta <- y %>% column_to_rownames("V1") %>%
    t() %>% data.frame() %>%
    rownames_to_column("Sample")
  abun <- z %>% column_to_rownames("V1")

  # filter tax by occurrence
  remain_tax <- apply(abun, 1, function(x){ sum(x>0)/length(x)}) %>% data.frame() %>%
            setNames("occurrence") %>%
            rownames_to_column("tax") %>%
            filter(occurrence > 0.4)
  tax_table <- abun[rownames(abun)%in%remain_tax$tax, ] %>% as.matrix()
  meta_data <- inner_join(x, meta, by = "Sample") %>% column_to_rownames("Sample")
  otu <- otu_table(tax_table, taxa_are_rows = TRUE)
  sample_data <- sample_data(meta_data)
  physeq <- phyloseq(otu, sample_data)

  return(physeq)

}
