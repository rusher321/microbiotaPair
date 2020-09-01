#' @importFrom data.table fread
#' @importFrom dplyr %>% inner_join
#' @importFrom phyloseq otu_table sample_data phyloseq
#' @importFrom tibble column_to_rownames rownames_to_column
#'
#' @title data2phyloseq
#' @description Convert all the variate such as phenotype metadata and the relative abundance profiles in the data.frame/matrix into the data of phyloseq type
#' @details 13/01/2020  ShenZhen China
#' @author  Hua Zou
#' @param metadata x  (Required) A table including multiple vectors of subjects' characteristic;
#' @param micro y  (Required) A table including multiple taxonomy, colnames are sampleid, rownames are taxonomy;
#' @param  occ A numberic of the occurrence above to restain
#' @usage data2phyloseq(x, y, z, sampleid)
#' @examples
#' library(phyloseq)
#' library(dplyr)
#' library(data.table)
#' library(tibble)
#'
#' metadata <- read.csv(system.file("extdata", "phenotype.csv", package="microbiotaPair"))
#' micro <- read.table(system.file("extdata", "abundance.profile", package="microbiotaPair"),sep="\t")
#' sampleid <- "SampleID"
#'
#' physeq_data <- data2phyloseq(metadata, micro)
#' #print(physeq_data)
#'
#' @return A phyloseq object
#'
#' @export data2phyloseq
#'

data2phyloseq <- function(metadata, micro, occ = 0.2){

  # filter tax by occurrence
  remain_tax <- apply(micro, 1, function(x){ sum(x>0)/length(x)}) %>% data.frame() %>%
            setNames("occurrence") %>%
            rownames_to_column("tax") %>%
            dplyr::filter(occurrence >= occ)
  tax_table <- micro[rownames(micro)%in%remain_tax$tax, ] %>% as.matrix()
  otu <- otu_table(tax_table, taxa_are_rows = TRUE)
  sample_data <- sample_data(metadata)
  physeq <- phyloseq(otu, sample_data)

  return(physeq)

}


