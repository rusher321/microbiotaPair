#' @title filterPer
#'
#' @param x , a data.frame object, species abundance profile
#' @param row, a numberic, 1/2
#' @param percent, percent to retain or remove
#' @param include, a logistic , T/F
#'
#' @return
#' data.frame
#'
#' @examples
#'
filterPer <- function (x, row, percent, include = T)
{
  if (include) {
    index <- apply(x, row, function(x) {
      (sum(x != 0)/length(x)) > percent
    })
  }
  else {
    index <- apply(x, row, function(x) {
      (sum(x != 0)/length(x)) < percent
    })
  }
  if (row == 1) {
    out <- x[index, ]
  }
  else {
    out <- x[, index]
  }
  return(out)
}




#' @title filterEx
#' @description filter the extreme value
#' @param x vector
#' @param percent.outlier numberic,defalut 0.05
#'
#' @return vector
#'
#' @examples
filterEx <- function(x, percent.outlier=0.05){
  num <- length(x)
  y <- sort(x)
  y2 <- y[round(num*percent.outlier):round(num*(1-percent.outlier))]
  return(y2)
}


#' @title multiplot
#' @description combine multi plot to one
#' @param plotlist
#' @param file
#' @param cols
#' @param layout
#'
#' @return
#' @export
#'
#' @examples
multiplot <- function(plotlist, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- plotlist

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}




#' @Title rocPlot
#' @description
#' @param response
#' @param predict
#'
#' @return
#' @export
#'
#' @examples
rocPlot <- function(response, predict, label = "AUC(low vs high) = "){

  roc1 <- pROC::roc(response, predict)
  roc1.ci <- round(as.numeric(ci.auc(response , predict)),2)*100
  title1 <- paste0(label, roc1.ci[2],"% (",roc1.ci[1],"-",roc1.ci[3],"%)")
  rocdat1 <- data.frame(roc1$sensitivities, 1-roc1$specificities, title1)
  names(rocdat1) <- c('sen', 'spe', 'group')

  label1 <- title1

  annosize <- 5.5

  p <- ggplot(rocdat1, aes(x=spe, y=sen, colour=group, group=group)) +
    xlab('1-Specificity') +
    ylab('Sensitivity') +
    geom_line(size=1.5) +
    mytheme +
    theme(legend.position = c(0.7,0.25)) +
    scale_colour_manual(values = c("#F3B670"),
                        labels=c(label1)) +
    guides(colour=guide_legend(title = NULL)) +
    geom_abline(linetype = 'dashed', colour = 'grey')

  return(p)

}



#' @title rocPlot2
#' @description  use the gg_roc to plot the roc curve
#' @param response
#' @param predict
#' @param title
#'
#' @return
#' @export
#'
#' @examples
rocPlot2 <- function(response, predict, title){

  roc1 <- pROC::roc(response, predict)
  roc1.ci <- round(as.numeric(ci.auc(response , predict)),2)*100
  title1 <- paste0("AUC(low vs high) = ",roc1.ci[2],"% (",roc1.ci[1],"-",roc1.ci[3],"%)")

  qdat <- data.frame(response = response, predict = predict, group = title1)

  label1 <- title1

  annosize <- 5.5

  p <- ggplot(qdat, aes(m = predict, d = as.numeric(qdat$response)-1, colour = group)) +
    geom_roc()+
    mytheme +
    theme(legend.position = c(0.7,0.25)) +
    scale_colour_manual(values = c("#000000"),
                        labels=c(label1)) +
    guides(colour=guide_legend(title = NULL)) +
    geom_abline(linetype = 'dashed', colour = 'grey')+ggtitle(title)

  return(p)

}


#' @Title myspearman
#' @description  to compute the correlation between the multi species and metadata
#' @param phe
#' @param species
#' @param method
#'
#' @return
#' @export
#'
#' @examples
myspearman <- function(phe, species, method = "s"){


  id <- intersect(rownames(phe), colnames(species))
  phe <- phe[id, ]
  species.g <- species[ ,id]

  pro.n <- nrow(species.g)
  phe.n <- ncol(phe)
  out.s <- matrix("NA", pro.n, phe.n*2)
  rownames(out.s) <- rownames(species.g)
  colnames(out.s) <- rep("pvalue", phe.n*2)

  for (i in 1:phe.n){

    x <- as.numeric(phe[,i] )
    print(length(x))
    colnames(out.s)[i*2-1] <- paste0(colnames(phe)[i], "_estimate")
    colnames(out.s)[i*2] <-   paste0(colnames(phe)[i], "_p.value")

    for (j in 1:pro.n){
      y <- as.numeric(species.g[j,])
      # rm NA
      id <- !is.na(x) & !is.na(y)
      cor.s <- cor.test(x[id],y[id],method = method)
      out.s[j,(i*2-1):(i*2)] <- c(cor.s$estimate, cor.s$p.value)

    }

  }

  out2 <- apply(out.s, 2, as.numeric)
  rownames(out2) <- rownames(out.s)
  return(out2)


}










