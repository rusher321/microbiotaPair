#' @importFrom dplyr %>% rename
#' @importFrom tidyr gather
#' @importFrom microbiomeSeq perform_anova
#' @importFrom phyloseq otu_table sample_data
#'
#'
#' @title Plot alpha diversity
#' @description This function calculates alpha diversity of provided community data using
#' selected method(s). It performs pair-wise ANOVA of diversity measures between groups
#' and outputs a plot for each of the selected methods(indices) annotated with significance levels.
#'
#' @details 14/01/2020  ShenZhen China
#' @author  Hua Zou
#' @param physeq (Required). A \code{phyloseq} object containing merged information of abundance,
#'        sample data including the measured variables and categorical information of the samples.
#' @param method (Required). A list of character strings specifying \code{method} to be used to calculate for alpha diversity
#'        in the data. Available methods are: "observed", "chao1", "simpson", "shannon", "evenness" and "all".
#' @param grouping_column (Required). A character string specifying the name of a categorical variable containing grouping information.
#' @param pid (Required) A character string for paired test
#'
#' @usage plot_alpha_diversity(physeq, method, grouping_column, pid)
#' @examples
#' library(dplyr)
#' library(ggplot2)

#' data(physeq_data)
#' method <- c("shannon", "simpson")
#' grouping_column <- "Stage"
#' pid <- "ID"
#' plot_alpha_diversity(physeq_data, method, grouping_column, pid)
#'
#' @return  Returns a ggplot object
#'
#' @export plot_alpha_diversity
#'
plot_alpha_diversity <- function(physeq, method, grouping_column, pid){

  abund_table <- phyloseq::otu_table(physeq)
  meta_table <- phyloseq::sample_data(physeq)
  div.df <- microbiotaPair::alpha_diversity(physeq, method) %>%
    tidyr::gather(key="measure", value=value, -sample)

  df <- data.frame(div.df, (meta_table[, c(grouping_column, pid)])[as.character(div.df$sample),])

  pr <- levels(df[, grouping_column])
  cmp <- NULL
  for(i in 1:(length(pr) -1 )){
    for(j in (i+1):length(pr)){
      tmp <- c(pr[i], pr[j])
      if(is.null(cmp)){
        cmp <- tmp
      }else{
        cmp <- list(cmp, tmp)
      }
    }
  }

  anova_res <- microbiomeSeq::perform_anova(df, meta_table, grouping_column, 1)
  df_pw <- anova_res$df_pw

  df <- df %>% dplyr::rename("Group"=grouping_column) %>% mutate(Group=factor(Group, levels = rev(pr)))
  p <- ggplot(df, aes(x = Group, y = value)) +
    stat_boxplot(geom = "errorbar",width = 0.15) +
    geom_boxplot(aes(fill = Group), width = 0.4,
                 outlier.colour = "black", outlier.shape = 21, outlier.size = 1) +
    stat_summary(fun.y = mean, geom = "point", shape = 16, size = 2, color = "black") +
    geom_jitter(position = position_jitter(height = 0, width=0), shape = 21) +
    geom_line(aes(group=ID)) +
    facet_wrap(~measure, scales="free_y", nrow=1) +
    labs(x = treat, y = "Observed Values") +
    guides(color=F, fil=F) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 10,color = "black", face = "bold", hjust = 0.5),
          axis.title = element_text(size = 10, color = "black",face = "bold"),
          axis.text = element_text(size = 9, color = "black"), axis.ticks.length = unit(-0.05, "in"),
          axis.text.y = element_text(margin = unit(c(0.3,0.3, 0.3, 0.3), "cm"), size = 9),
          axis.text.x = element_text(margin = unit(c(0.3,0.3, 0.3, 0.3), "cm")),
          text = element_text(size = 8,color = "black"),
          strip.text = element_text(size = 9, color = "black", face = "bold"),
          panel.grid = element_blank())

  if(!is.null(df_pw)){
    for(i in 1:dim(df_pw)[1]){
      p <- p + geom_path(inherit.aes=F,aes(x,y),
                       data=data.frame(x = c(which(levels(df[,"Group"])==as.character(df_pw[i,"from"])),
                                               which(levels(df[,"Group"])==as.character(df_pw[i,"to"]))),
                                         y = c(as.numeric(as.character(df_pw[i,"y"])),as.numeric(as.character(df_pw[i,"y"]))),
                                         measure=c(as.character(df_pw[i,"measure"]), as.character(df_pw[i,"measure"]))),
                       color="black",lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.05, "inches")))
      p <- p + geom_text(inherit.aes=F,aes(x=x,y=y,label=label),
                       data=data.frame(x=(which(levels(df[,"Group"])==as.character(df_pw[i,"from"]))+
                                            which(levels(df[,"Group"])==as.character(df_pw[i,"to"])))/2,
                                       y=as.numeric(as.character(df_pw[i,"y"])),
                                       measure=as.character(df_pw[i,"measure"]),
                                       label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),
                                              breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),
                                              label=c("***", "**", "*", "ns")))))
    }
  }

  return(p)
}
