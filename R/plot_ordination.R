#' @import ggplot2
#' @import phyloseq
#' @import dplyr
#'
#'
#' @title Ordination plots
#'
#' @description This function produces visualisation of ordination and beta dispersion results.
#'
#' @details 14/01/2020  ShenZhen China
#' @author  Hua Zou
#'
#' @param ordination.res (Required).  A solution of ordination results returned from \link[microbiotaPair]{ordination}.
#' @param method (Optional). A character string specifying ordination method. All methods available to the \code{ordinate} function
#'               of \code{phyloseq} are acceptable here as well ( Default is "PCOA").
#'               "DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA", "PCA", "Tsne".
#' @param grouping_column (Required). Character string specifying name of a categorical variable that is preffered for grouping the information.
#'        information.
#' @param PID   id (Required) for paired test
#'
#' @return Returns a ggplot object.
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#'
#' @usage plot_ordination(physeq, method, grouping_column)
#' @examples
#'
#' data(physeq_data)
#' physeq <- physeq_data
#' ord.res <- ordination(physeq, method="NMDS", grouping_column="Stage")
#' plot_ordination(ord.res, method="PCoA", grouping_column="Stage", PID="ID")
#'
#'
#' @export plot_ordination
#'
plot_ordination <- function(ordination.res, method, grouping_column="Stage", PID="ID"){

  sol <- ordination.res$solution
  adn_res <- ordination.res$adonis_res
  betadisper_res <- ordination.res$betadispersion

  df_points <- data.frame(sol$points) %>% setNames(c("Axis1", "Axis2"))
  df_meta <- phyloseq::sample_data(physeq) %>% select(c(PID, grouping_column)) %>%
    data.frame() %>% setNames(c("PID", "grouping"))
  df_mdat <- inner_join(df_points %>% rownames_to_column("sample"),
                        df_meta %>% rownames_to_column("sample"),
                        by = "sample")

  #coloring function
  gg_color_hue <- function(n){
    hues <- seq(15,375,length=n+1)
    hcl(h=hues,l=65,c=100)[1:n]
  }

  cols <- gg_color_hue(length(unique(df_mdat$grouping)))

  p <- ggplot(df_mdat, aes(x=Axis1, y=Axis2))+
    geom_point(aes(color = grouping))+
    scale_color_manual(values = cols)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'),
          legend.key = element_rect(fill = 'transparent'))

  group_label <- cbind(Axis1=tapply(df_mdat$Axis1, df_mdat$grouping, mean),
                       Axis2=tapply(df_mdat$Axis2, df_mdat$grouping, mean)) %>%
    data.frame() %>%
    rownames_to_column("grouping")
  group_border <- plyr::ddply(df_mdat, 'grouping', function(x)x[chull(x[[2]], x[[3]]), ])

  p <- p + geom_line(aes(group=PID), linetype = "dashed", alpha = 0.3) +
    geom_text(data = group_label, aes(x=Axis1, y=Axis2, label=grouping, color=grouping)) +
    geom_polygon(data = group_border, aes(fill = grouping), color = "black", alpha = 0.1, show.legend = FALSE)+
    scale_color_manual(values = cols)+
    guides(group=F, fill=F, color=F)+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")

  #axis labels
  if(method=="NMDS"){
    #annotate plot with stress value from NMDS results
    stress.value <- sol$stress
    stress.label <- paste("STRESS=",round(stress.value,4))
    p <- p +     p2 <- p + annotate("text", x = max(df_mdat$Axis1),
                                    y = max(df_mdat$Axis2),
                                    label = stress.label)
    p <- p + xlab("NMDS1")+ ylab("NMDS2")
  }else if(method == "PCoA"){
    eig_var <- (ord.res$solution$eig[1:2]/sum(ord.res$solution$eig))*100
    p <- p + xlab(paste("PCoA1 (",sprintf("%.4g", eig_var[1]),"%)",sep="")) +
             ylab(paste("PCoA2 (",sprintf("%.4g",eig_var[2]),"%)",sep=""))
  }

  #add the adonis results on the plot using custom annoatation
  #this only happens if adonis results turn out significant
  if(!is.null(adn_res)){
    adn_pvalue<-adn_res[[1]][["Pr(>F)"]][1]
    adn_rsquared<-round(adn_res[[1]][["R2"]][1],3)
    #use the bquote function to format adonis results to be annotated on the ordination plot.
    signi_label <- paste(cut(adn_pvalue,breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ".")))
    adn_res_format <- bquote(atop(atop("PERMANOVA",R^2==~.(adn_rsquared)), atop("p-value="~.(adn_pvalue)~.(signi_label), phantom())))

    #annotate plot with adonis results
    p <- p + annotate("text", x = max(df_mdat$Axis1),
                       y = min(df_mdat$Axis2),
                       label = adn_res_format)
  }

  # add a table of beta dispersion results
  if(!is.null(betadisper_res)){
    anova_label <- NULL
    for(i in 1:nrow(betadisper_res)){
      tmp <- paste("groups=", betadisper_res[i, 1], " p_value=", signif(betadisper_res[i, 2], 3), " label=", betadisper_res[i, 3], sep = " ")
      if(is.null(anova_label)){
        anova_label <- tmp
      }else{
        anova_label <- paste(anova_label, tmp, sep = "\n")
      }
      p <- p + ggtitle(anova_label)
    }

}
  return(p)
}
