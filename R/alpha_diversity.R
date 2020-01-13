#' @importFrom dplyr inner_join
#' @importFrom vegan diversity
#' @importFrom microbiome alpha
#' @importFrom phyloseq otu_table
#'
#'
#' @title Alpha diversity
#'
#' @description The alpha_diversity aims to calculate the global indicators of the ecosystem state
#' @details 10/01/2020  ShenZhen China
#' @author  Hua Zou
#' @param physeq (Required). A \code{phyloseq} object containing merged information of abundance,
#'        sample data including the measured variables and categorical information of the samples.
#' @param method (Required). A list of character strings specifying \code{method} to be used to calculate for alpha diversity
#'        in the data. Available methods are: "observed", "chao1", "simpson", "shannon" and "evenness" .
#'
#' @usage alpha_diversity(physeq, method)
#' @examples
#' data(physeq_data)
#' method <- "all"
#' alpha_diversity(physeq_data, method)
#'
#' @return  It returns a data frame of diversity measure and corresponding indices/methods
#'
#' @export alpha_diversity
#'
alpha_diversity <- function(physeq, method){

  #==check for validity of selected methods
  method <- match.arg(method, c("observed", "chao1", "simpson", "shannon", "evenness", "all"), several.ok = TRUE)

  microbiome_alpha <- function(method){
    abund_table <- phyloseq::otu_table(physeq_data)
    message(method)
    malpha <- microbiome::alpha(abund_table, index = method)
    df_malpha <- data.frame(sample=rownames(malpha), value=malpha)
    return(df_malpha)
  }

  microbiome_alpha2 <- function(method){
    abund_table <- t(phyloseq::otu_table(physeq_data))
    message(method)
    malpha <- vegan::diversity(abund_table, method)
    df_malpha <- data.frame(sample=names(malpha), value=malpha)
    colnames(df_malpha)[2] <- method
    return(df_malpha)
  }

  df <- NULL
  if(any(c("all", "observed") %in% method)){
    observed <- microbiome_alpha("observed")
    if(is.null(df)){df <- observed}
    else{df <- dplyr::inner_join(df, observed, by = "sample")}
   }

  if(any(c("all", "chao1") %in% method)){
    chao1 <- microbiome_alpha("chao1")
    if(is.null(df)){df <- chao1}
    else{df <- dplyr::inner_join(df, chao1, by = "sample")}
  }

  if(any(c("all", "simpson") %in% method)){
    simpson <- microbiome_alpha2("simpson")
    if(is.null(df)){df <- simpson}
    else{df <- dplyr::inner_join(df, simpson, by = "sample")}
    }

  if(any(c("all", "shannon") %in% method)){
    shannon <- microbiome_alpha2("shannon")
    if(is.null(df)){df <- shannon}
    else{df <- dplyr::inner_join(df, shannon, by = "sample")}
    }

  if(any(c("all", "evenness") %in% method)){
      message(method)
      H <- vegan::diversity(t(phyloseq::otu_table(physeq_data)))
      S <- vegan::specnumber(t(phyloseq::otu_table(physeq_data)))
      J <- H/log(S)
      evenness <- data.frame(sample=names(J), evenness=J)
      if(is.null(df)){
        df <- evenness}
      else {
        df <- dplyr::inner_join(df, evenness, by = "sample")}
    }

  return(df)
}
