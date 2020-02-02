#' plot_centrel
#'
#' @param pro a data.frame or matrix, species abundance
#' @param method a character, distance method,eg "bray", "euclidean"
#' @param config a data.frame include sample metadata
#' @param color vector, group color
#'
#' @return
#' ggplot object
#' @export
#'
#' @examples
plot_centrel <- function(pro, method , config, color = time_colour){

  id <- intersect(rownames(pro), rownames(config))
  pro <- pro[id,]
  config <- config[id, time_varname]
  # compute the distance

  prodis <- vegan::vegdist(pro, method = method)
  mod <- betadisper(prodis, config)

  qdata <- data.frame(dis = mod$distance, label = config)

  # plot
  my_comparisons = list()
  num <- combn(length(unique(config)),2)
  for(i in 1:ncol(num)){my_comparisons[[i]] <- num[,i]}

  p <- ggboxplot(qdata, x="label", y = "dis", color = "label" ,add = "jitter",alpha=0.6,size = 0.5)+
    stat_compare_means(comparisons = my_comparisons)+
    scale_color_manual(values=color)+
    theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 45),legend.position = "none")+xlab("")+
    ylab("Distance to centroid")

  return(p)

}
