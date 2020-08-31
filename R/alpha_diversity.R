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
#' @author(s)  Hua Zou, Huahui Ren huahui.ren@bio.ku.dk
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
alpha_diversity <- function(physeq_data, method, paired = T){

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

  # rbind the sample ID
  if(paired){
    sampledata <- sample_data(physeq_data) %>% select(sample_varname, time_varname, pairID_varname) %>%
    rename(sample = sample_varname) %>% inner_join(df, by = "sample")
  }else{
    sampledata <- sample_data(physeq_data)
    sampledata$sample <- rownames(sampledata)
    sampledata <- inner_join(df, sampledata, by = "sample")
}
  return(sampledata)
}

#############################################

#' @title hill_number
#'
#' @description The hill_number aims to calculate the hill number
#' @details 10/01/2020  ShenZhen China
#' @author(s)  Hua Zou, Huahui Ren
#'
#' @param phyloseq ,  A \code{phyloseq} object containing merged information of abundance,
#'        sample data including the measured variables and categorical information of the samples.
#' @param q , A vevtor eg c(0:1), c(0:3)
#' @param samID , A character indicate the sample ID is in row or col; include "row" or "col"
#'
#' @usage hill_number(physeq, c(1:3))
#' @examples
#' data(physeq_data)
#' hill_number(physeq_data, c(1:3))
#' @return
#' data.frame
#' @export
#'

hill_number <- function(phyloseq, q, samID = "row"){

  # to get the inf from phyloseq
  otu <- as.data.frame(otu_table(phyloseq))
  sam <- as.data.frame(sample_data(phyloseq))
  # to generate the hill result
  if(samID=="row"){
    res <- sapply(q, function(x){res <- hill_taxa(t(otu), q=x)})
  }else{
    res <- sapply(q, function(x){res <- hill_taxa(otu, q=x)})
  }
  # to rbind the sample inf
  colnames(res) <- paste0("q=", q)
  hillres <- cbind(res, sam[rownames(res), c(time_varname, pairID_varname)])

  return(hillres)

}

###############################################
#' ggplot2 extension for an iNEXT object
#' this function focked from the iNEXT package
#' \code{ggiNEXT}: the \code{\link[ggplot2]{ggplot}} extension for \code{\link{iNEXT}} Object to plot sample-size- and coverage-based rarefaction/extrapolation curves along with a bridging sample completeness curve
#' @param x an \code{iNEXT} object computed by \code{\link{iNEXT}}.
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1});
#' sample completeness curve (\code{type = 2}); coverage-based rarefaction/extrapolation curve (\code{type = 3}).
#' @param se a logical variable to display confidence interval around the estimated sampling curve.
#' @param facet.var create a separate plot for each value of a specified variable:
#'  no separation \cr (\code{facet.var="none"});
#'  a separate plot for each diversity order (\code{facet.var="order"});
#'  a separate plot for each site (\code{facet.var="site"});
#'  a separate plot for each combination of order x site (\code{facet.var="both"}).
#' @param color.var create curves in different colors for values of a specified variable:
#'  all curves are in the same color (\code{color.var="none"});
#'  use different colors for diversity orders (\code{color.var="order"});
#'  use different colors for sites (\code{color.var="site"});
#'  use different colors for combinations of order x site (\code{color.var="both"}).
#' @param grey a logical variable to display grey and white ggplot2 theme.
#' @param usercolor a vector including the color that user define
#' @param ... other arguments passed on to methods. Not currently used.
#' @return a ggplot2 object
#' @examples
#' library(iNEXT)
#' data(spider)
#' # single-assemblage abundance data
#' out1 <- iNEXT(spider$Girdled, q=0, datatype="abundance")
#' ggiNEXT(x=out1, type=1)
#' ggiNEXT(x=out1, type=2)
#' ggiNEXT(x=out1, type=3)
#'
#'\dontrun{
#' # single-assemblage incidence data with three orders q
#' data(ant)
#' size <- round(seq(10, 500, length.out=20))
#' y <- iNEXT(ant$h500m, q=c(0,1,2), datatype="incidence_freq", size=size, se=FALSE)
#' ggiNEXT(y, se=FALSE, color.var="order")
#'
#' # multiple-assemblage abundance data with three orders q
#' z <- iNEXT(spider, q=c(0,1,2), datatype="abundance")
#' myggiNEXT.iNEXT(z, facet.var="site", color.var="order")
#'}
#' @export
#'



myggiNEXT.iNEXT <- function(x, type=1, se=TRUE, facet.var="none", color.var="site", grey=FALSE,
                          usercolor=NULL){
  TYPE <-  c(1, 2, 3)
  SPLIT <- c("none", "order", "site", "both")
  if(is.na(pmatch(type, TYPE)) | pmatch(type, TYPE) == -1)
    stop("invalid plot type")
  if(is.na(pmatch(facet.var, SPLIT)) | pmatch(facet.var, SPLIT) == -1)
    stop("invalid facet variable")
  if(is.na(pmatch(color.var, SPLIT)) | pmatch(color.var, SPLIT) == -1)
    stop("invalid color variable")

  type <- pmatch(type, 1:3)
  facet.var <- match.arg(facet.var, SPLIT)
  color.var <- match.arg(color.var, SPLIT)
  if(facet.var=="order") color.var <- "site"
  if(facet.var=="site") color.var <- "order"

  options(warn = -1)
  z <- fortify(x, type=type)
  options(warn = 0)
  if(ncol(z) ==7) {se <- FALSE}
  datatype <- unique(z$datatype)
  if(color.var=="none"){
    if(levels(factor(z$order))>1 & "site"%in%names(z)){
      warning("invalid color.var setting, the iNEXT object consists multiple sites and orders, change setting as both")
      color.var <- "both"
      z$col <- z$shape <- paste(z$site, z$order, sep="-")

    }else if("site"%in%names(z)){
      warning("invalid color.var setting, the iNEXT object consists multiple orders, change setting as order")
      color.var <- "site"
      z$col <- z$shape <- z$site
    }else if(levels(factor(z$order))>1){
      warning("invalid color.var setting, the iNEXT object consists multiple sites, change setting as site")
      color.var <- "order"
      z$col <- z$shape <- factor(z$order)
    }else{
      z$col <- z$shape <- rep(1, nrow(z))
    }
  }else if(color.var=="order"){
    z$col <- z$shape <- factor(z$order)
  }else if(color.var=="site"){
    if(!"site"%in%names(z)){
      warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
      z$col <- z$shape <- factor(z$order)
    }
    z$col <- z$shape <- z$site
  }else if(color.var=="both"){
    if(!"site"%in%names(z)){
      warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
      z$col <- z$shape <- factor(z$order)
    }
    z$col <- z$shape <- paste(z$site, z$order, sep="-")
  }
  zz=z
  z$method[z$method=="observed"]="interpolated"
  z$lty <- z$lty <- factor(z$method, levels=unique(c("interpolated", "extrapolated"),
                                                   c("interpolation", "interpolation", "extrapolation")))
  z$col <- factor(z$col)
  data.sub <- zz[which(zz$method=="observed"),]

  g <- ggplot(z, aes_string(x="x", y="y", colour="col")) +
    geom_point(aes_string(shape="shape"), size=5, data=data.sub)


  g <- g + geom_line(aes_string(linetype="lty"), lwd=1.5) +
    guides(linetype=guide_legend(title="Method"),
           colour=guide_legend(title="Guides"),
           fill=guide_legend(title="Guides"),
           shape=guide_legend(title="Guides")) +
    theme(legend.position = "bottom",
          legend.title=element_blank(),
          text=element_text(size=18),
          legend.key.width = unit(1.2,"cm"))

  if(type==2L) {
    g <- g + labs(x="Number of sampling units", y="Sample coverage")
    if(datatype=="abundance") g <- g + labs(x="Number of individuals", y="Sample coverage")
  }
  else if(type==3L) {
    g <- g + labs(x="Sample coverage", y="Species diversity")
  }
  else {
    g <- g + labs(x="Number of sampling units", y="Species diversity")
    if(datatype=="abundance") g <- g + labs(x="Number of individuals", y="Species diversity")
  }

  if(se)
    g <- g + geom_ribbon(aes_string(ymin="y.lwr", ymax="y.upr", fill="factor(col)", colour="NULL"), alpha=0.2)


  if(facet.var=="order"){
    if(length(levels(factor(z$order))) == 1 & type!=2){
      warning("invalid facet.var setting, the iNEXT object do not consist multiple orders.")
    }else{
      g <- g + facet_wrap(~order, nrow=1, scales = "free_y")
      if(color.var=="both"){
        g <- g + guides(colour=guide_legend(title="Guides", ncol=length(levels(factor(z$order))), byrow=TRUE),
                        fill=guide_legend(title="Guides"))
      }
    }
  }

  if(facet.var=="site"){
    if(!"site"%in%names(z)) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple sites.")
    }else{
      g <- g + facet_wrap(~site, nrow=1)
      if(color.var=="both"){
        g <- g + guides(colour=guide_legend(title="Guides", nrow=length(levels(factor(z$order)))),
                        fill=guide_legend(title="Guides"))
      }
    }
  }

  if(facet.var=="both"){
    if(length(levels(factor(z$order))) == 1 | !"site"%in%names(z)){
      warning("invalid facet.var setting, the iNEXT object do not consist multiple sites or orders.")
    }else{
      g <- g + facet_wrap(site~order)
      if(color.var=="both"){
        g <- g +  guides(colour=guide_legend(title="Guides", nrow=length(levels(factor(z$site))), byrow=TRUE),
                         fill=guide_legend(title="Guides"))
      }
    }
  }

  if(grey){
    g <- g + theme_bw(base_size = 18) +
      scale_fill_grey(start = 0, end = .4) +
      scale_colour_grey(start = .2, end = .2) +
      guides(linetype=guide_legend(title="Method"),
             colour=guide_legend(title="Guides"),
             fill=guide_legend(title="Guides"),
             shape=guide_legend(title="Guides")) +
      theme(legend.position="bottom",
            legend.title=element_blank())
  }
  if(!is.null(usercolor)){
    g <- g + theme_bw(base_size = 18) +
      scale_fill_manual(values = usercolor) +
      scale_colour_manual(values = usercolor) +
      guides(linetype=guide_legend(title="Method"),
             colour=guide_legend(title="Guides"),
             fill=guide_legend(title="Guides"),
             shape=guide_legend(title="Guides")) +
      theme(legend.position="bottom",
            legend.title=element_blank(),
            #panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
           )
  }
  g <- g + theme(legend.box = "vertical")
  return(g)

}



