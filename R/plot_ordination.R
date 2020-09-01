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
#' @author  Hua Zou; hua huiren
#'
#' @param ordination.res (Required).  A solution of ordination results returned from \link[microbiotaPair]{ordination}.
#' @param method (Optional). A character string specifying ordination method. All methods available to the \code{ordinate} function
#'               of \code{phyloseq} are acceptable here as well ( Default is "PCOA").
#'               "DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA", "PCA", "Tsne".
#' @param grouping_column (Required). Character string specifying name of a categorical variable that is preffered for grouping the information.
#'        information.
#' @param PID   id (Required) for paired test
#'
#' @param time_colour  a vector include the color
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
plot_ordination <- function(ordination.res, phylores, method, grouping_column="Stage",
                            PID="ID" , paired = F, time_colour = NULL){

  if(method == "Tsne"){
    temp <- data.frame(ordination.res$solution$tsne$par)
    rownames(temp) <- ordination.res$solution$names
    sol <- list(solution = temp)
  }else{
    sol <- ordination.res$solution
  }
  adn_res <- ordination.res$adonis_res
  betadisper_res <- ordination.res$betadispersion
  df_points <- data.frame(sol$points) %>% setNames(c("Axis1", "Axis2"))
  if(paired){
    df_meta <- phyloseq::sample_data(phylores) %>% select(c(PID, grouping_column)) %>%
      data.frame() %>% setNames(c("PID", "grouping"))
  }else{
    df_meta <- phyloseq::sample_data(phylores)[, c(grouping_column), drop=F] %>%
      data.frame() %>% setNames(c("grouping"))
  }

  df_mdat <- inner_join(df_points %>% rownames_to_column("sample"),
                        df_meta %>% rownames_to_column("sample"),
                        by = "sample")

  #coloring function
  if(is.null(time_colour)){
    gg_color_hue <- function(n){
      hues <- seq(15,375,length=n+1)
      hcl(h=hues,l=65,c=100)[1:n]
    }

    cols <- gg_color_hue(length(unique(df_mdat$grouping)))
  }else{
    cols <- time_colour
  }

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

  if(paired){
    p <- p + geom_line(aes(group=PID), linetype = "dashed", alpha = 0.3)
  }
    p <- p + geom_text(data = group_label, aes(x=Axis1, y=Axis2, label=grouping, color=grouping)) +
    geom_polygon(data = group_border, aes(fill = grouping), color = "black", alpha = 0.1, show.legend = FALSE)+
    #scale_color_manual(values = cols)+
    guides(group=F, fill=F, color=F)+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")

  #axis labels
  if(method=="NMDS"){
    #annotate plot with stress value from NMDS results
    stress.value <- sol$stress
    stress.label <- paste("STRESS=",round(stress.value,4))
    p <- p + annotate("text", x = max(df_mdat$Axis1),
                                    y = max(df_mdat$Axis2),
                                    label = stress.label)
    p <- p + xlab("NMDS1")+ ylab("NMDS2")
  }else if(method == "PCoA"){
    eig_var <- (ord.res$solution$eig[1:2]/sum(ord.res$solution$eig))*100
    p <- p + xlab(paste("PCoA1 (",sprintf("%.4g", eig_var[1]),"%)",sep="")) +
             ylab(paste("PCoA2 (",sprintf("%.4g",eig_var[2]),"%)",sep=""))
  }else if(method == "Tsne"){
    p <- p + xlab("Tsne1")+ ylab("Tsne2")
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
  return(p)
}


############# plot the tax composition  ###########


#' comtaxTop
#'
#' @param dat
#' @param group
#' @param top
#' @param group_var
#'
#' @return
#' @export
#'
#' @examples
comtaxTop <- function(dat, group, top, group_var){

  # order features
  library(tidyverse)
  b <- data.frame(
    Max = apply(dat,1,mean)
  ) %>%
    rownames_to_column(var="Species") %>%
    arrange(desc(Max))
  order_feature <- b$Species[1:top]
  if(nrow(b)>top){
    dat <- dat %>%
      rownames_to_column(var="Species") %>%
      mutate(
        Species = replace(Species, !Species %in% order_feature, "Other")
      ) %>%
      group_by(Species)  %>%
      summarise_all(sum) %>%
      column_to_rownames(var="Species")
    order_feature <- c(order_feature, "Other")
  }

  # order samples
  a <- data.frame(Sum = as.numeric(dat[order_feature[1],]))
  rownames(a) <- colnames(dat)
  group <- group[rownames(a), group_varname, drop=F]
  colnames(group) <- "Group"
  a <- cbind(a,group) %>% rownames_to_column(var="SampleID")
  a <- a[order(a$Sum,decreasing = T),]
  a <- a[order(a$Group),]
  order_sampels <- as.character(a$SampleID)
  a$y = ""
  a$SampleID <- factor(a$SampleID, levels = order_sampels)

  # Dataframe transform
  dat2 <- dat %>%
    rownames_to_column(var="Species") %>%
    gather(Samples,Abundance,-Species) %>%
    mutate(
      Species = factor(Species,levels = order_feature),
      Samples = factor(Samples,levels = order_sampels)
    )

  # barplot
  # Colors <- c("#000080","#0029FF","#00D5FF","#7DFF7A","#FFE600","#FF4700","#800000","#808080")
  Colors <- brewer.pal(12, "Paired")  #
  if(top>11){
    stop("add the colors scale!")
  }
  p1 <- ggplot(dat2,aes(Samples,Abundance,fill=Species))+
    geom_bar(stat = "identity",color="black")+
    scale_y_continuous(expand = expand_scale(mult = c(0,0.1)))+
    scale_fill_manual(values = Colors)+
    xlab("") + ylab("Relative abundance")+
   # labs(title = FullName)+
    theme_bw()+
    theme(
      axis.title = element_text(size = 14,color = "black"),
      axis.text = element_text(size = 13,color = "black"),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text(size = 14,color = "black"),
      legend.text = element_text(size = 13,color = "black",face = "italic"),
      panel.grid = element_blank(),
      plot.title = element_text(size = 14,color = "black",hjust = 0.5)
    )

  p2 <- ggplot(a,aes(SampleID,y,fill=Group))+
    geom_tile()+
    xlab("Samples")+ylab("")+
    scale_y_discrete(expand = c(0,0))+
    scale_fill_manual(values = time_colour)+
    theme_bw()+
    theme(
      axis.title = element_text(size = 14,color = "black"),
      axis.text = element_text(size = 13,color = "black"),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text(size = 14,color = "black"),
      legend.text = element_text(size = 13,color = "black"),
      panel.grid = element_blank()
    )

  ggarrange(p1,p2,nrow = 2,ncol = 1,align = "v",heights = c(5,1),legend = "right")
}























