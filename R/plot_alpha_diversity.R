#' @importFrom tidyr gather
#' @importFrom phyloseq otu_table sample_data
#'
#'
#' @title Plot alpha diversity
#' @description This function calculates alpha diversity of provided community data using
#' selected method(s). It performs paired wilcox or t-test of diversity measures between groups
#' and outputs a plot for each of the selected methods(indices) annotated with significance levels.
#'
#' @details 22/01/2020  ShenZhen China
#' @author  Huahui Ren, zouhua
#' @param alphares (Required). A result from \code{alpha_diversity} function
#'        sample data including the measured variables and categorical information of the samples.
#' @param method (Required). A list of character strings specifying \code{method} to be used to calculate for alpha diversity
#'        in the data. Available methods are: "t.test", "wilcox.test".
#'
#' @usage plot_alpha_diversity(alphares)
#' @examples
#' library(dplyr)
#' library(ggplot2)

#' source("R/config.R")
#' data(physeq_data)
#' alphares <- alpha_diversity(physeq_data)
#' plot_alpha_diversity(physeq_data, method, grouping_column, pid)
#'
#' @return  Returns a ggplot object
#'
#' @export plot_alpha_diversity
#'

plot_alpha_diversity <- function(alphares, method = "wilcox.test"){


    qdat <- gather(alphares, key = "alpha_index", value = "index_value",
                   -c(sample, pairID_varname, time_varname))
    # order the time_varname
    qdat[, time_varname] <- factor(qdat[, time_varname], levels = time_name)

    pr <- levels(qdat[, time_varname])
    mycmp <- list(c(pr))

    p <- ggboxplot(qdat, x = time_varname, y = "index_value", fill = time_varname) +
    geom_line(aes(group = ID), colour = "grey") +
    stat_compare_means(comparisons = mycmp, method = method, paired = TRUE) +
   # stat_boxplot(geom = "errorbar",width = 0.15) +
    #stat_summary(fun.y = mean, geom = "point", shape = 16, size = 2, color = "black") +
    geom_jitter(position = position_jitter(height = 0, width=0), shape = 21) +
    facet_wrap(~alpha_index, scales="free_y", nrow=1) +
    scale_y_continuous(expand = expand_scale(mult=c(0.1,0.1)))+ # add y axis
    guides(color=F, fil=F) + scale_fill_manual(values = time_colour)+mytheme

    return(p)

  }

