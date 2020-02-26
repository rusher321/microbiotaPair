#' @import dplyr
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom varhandle unfactor
#'
#'
#' @title Wilcoxon Sign-Rank Test
#'
#' @description
#' Tests for the difference between two related variables; takes into account the magnitude and direction of difference
#'
#' @details 01/10/2020 ShenZhen China
#' @author  Hua Zou/ Huahui Ren
#'
#' @param physeq (Required). A \code{phyloseq} object containing merged information of abundance,
#'        sample data including the measured variables and categorical information of the samples.
#' @param sampleid (Required) A character of the sampleid to connect phenotype and  profiles.
#' @param PID   id (Required) for paired test
#' @param GROUP names (Required) of group information, only contain two levels if grp1 or grp2 haven't been provided
#' @param grp1  one of groups to be converted into 0 (optional)
#' @param grp2  one of groups to be converted into 1 (optional)
#'
#' @usage wilcox_sign(physeq, sampleid, PID, GROUP, grp1, grp2)
#' @examples
#'
#' data(physeq_data)
#' sampleid <- "SampleID"
#' pid <- "ID"
#' group <- "Stage"
#' grp1 <- "Before"
#' grp2 <- "After"
#'
#' result <- wilcox_sign(physeq_data, sampleid, pid, group, grp1, grp2)
#'
#'
#' @return  Returns a result of Wilcoxon Sign-Rank Test
#' @return  type:       kind of data
#' @return  Block:      group information
#' @return  Num:        number of group
#' @return  P-value:    P by Wilcoxon Sign-Rank Test
#' @return  FDR:        adjusted by BH
#' @return  Enrichment: directory by median or directory by rank
#' @return  Occurence:  occurence of two groups
#' @return  median:     both or each group
#' @return  rank:       each group
#' @return  FDR:        adjusted P value by BH
#' @return  Odds Ratio:     95% Confidence interval
#'
#' @export
#'
wilcox_sign <- function(physeq, DNAID, PID, GROUP,
                        grp1=NULL, grp2=NULL){

  phen <- microbiome::meta(physeq) %>% tibble::rownames_to_column(DNAID)
  prof <- microbiome::abundances(physeq) %>% data.frame()

  # determine x with two cols and names are corret
  colnames(phen)[which(colnames(phen) == DNAID)] <- "SampleID"
  colnames(phen)[which(colnames(phen) == PID)] <- "ID"
  colnames(phen)[which(colnames(phen) == GROUP)] <- "Stage"
  if (length(which(colnames(phen)%in%c("SampleID","ID","Stage"))) != 3){
    warning("x without 2 cols: DNAID, ID, GROUP")
  }
  phe <- phen %>% select(c("SampleID", "ID", "Stage"))

  # select groups
  if(length(grp1)){
    phe.cln <- phe %>% filter(Stage%in%c(grp1, grp2)) %>%
      mutate(Stage=factor(Stage, levels = c(grp1, grp2))) %>%
      arrange(ID, Stage)
    pr <- c(grp1, grp2)
  } else {
    phe.cln <- phe %>% mutate(Stage=factor(Stage)) %>%
      arrange(ID, Stage)
    pr <- levels(phe.cln$Stage)
  }

  if (length(levels(phe.cln$Stage)) > 2) {
    stop("The levels of `group` are more than 2")
  }

  # profile
  sid <- intersect(phe.cln$SampleID, colnames(prof))
  prf <- prof %>% select(sid) %>%
    rownames_to_column("tmp") %>%
    # occurrence of rows more than 0.3
    filter(apply(select(., -one_of("tmp")), 1, function(x){sum(x[!is.na(x)] != 0)/length(x)}) > 0.3) %>%
    data.frame() %>% column_to_rownames("tmp") %>%
    t() %>% data.frame()

  # judge no row of profile filter
  if (ncol(prf) == 0) {
    stop("No row of profile to be choosed\n")
  }

  # determine the right order and group levels
  for(i in 1:nrow(prf)){
    if ((rownames(prf) != phe.cln$SampleID)[i]) {
      stop(paste0(i, " Wrong"))
    }
  }

  # merge phenotype and profile
  mdat <- inner_join(phe.cln %>% filter(SampleID%in%sid),
                     prf %>% rownames_to_column("SampleID"),
                     by = "SampleID")

  dat.phe <- mdat %>% select(c(1:3))
  dat.prf.tmp <- mdat %>% select(-c(1:3))
  dat.prf <- apply(dat.prf.tmp, 2, function(x){
    as.numeric(as.character(x))}) %>% data.frame()

  res <- apply(dat.prf, 2, function(x, grp){

    origin <- data.frame(value=as.numeric(x), grp)
    number <- tapply(origin$value, origin$Stage, function(x){sum(!is.na(x))})
    Num <- paste0(pr[1], number[1], "_vs_",
                  pr[2], number[2])
    # remove NA data
    dat <- origin %>% na.omit()
    intersectFun <- function(x){
      tmp <- x %>% mutate(Stage = factor(Stage))
      id <- unique(as.character(tmp$ID))
      for (i in 1:length(levels(tmp$Stage))) {
        id <- intersect(id,
                        unlist(tmp %>% filter(Stage == levels(Stage)[i]) %>% select(ID)))
      }
      return(id)
    }
    dat.cln <- dat %>% filter(ID%in%intersectFun(dat)) %>%
      arrange(ID, Stage)

    p <- signif(wilcox.test(value ~ Stage, data=dat.cln, paired=T)$p.value, 6)

    # median
    md <- signif(median(dat.cln$value), 4)
    mdn <- signif(tapply(dat.cln$value, dat.cln$Stage, median), 4)
    if ( mdn[1] > mdn[2] & p < 0.05) {
      enrich1 <- pr[1]
    } else if (mdn[1] < mdn[2] & p < 0.05) {
      enrich1 <- pr[2]
    } else if (p > 0.05 | mdn[1] == mdn[2]){
      enrich1 <- "No significance"
    }

    # rank
    rk <- rank(dat.cln$value)
    rnk <- signif(tapply(rk, dat.cln$Stage, mean), 4)
    if ( rnk[1] > rnk[2] & p < 0.05) {
      enrich2 <- pr[1]
    } else if (rnk[1] < rnk[2] & p < 0.05) {
      enrich2 <- pr[2]
    } else if (p > 0.05 | rnk[1] == rnk[2]){
      enrich2 <- "No significance"
    }
    occ <- signif(tapply(dat.cln$value, dat.cln$Stage, function(x){
      round(sum(x > 0)/length(x), 4)}), 4)
    Pair <- nrow(dat.cln)

    res <- c(Num,Pair,p,enrich1,enrich2,occ,md,mdn,rnk)

    return(res)
  }, dat.phe) %>%
    t(.) %>% data.frame(.) %>%
    rownames_to_column("type") %>%
    varhandle::unfactor(.)

  colnames(res)[2:ncol(res)] <- c("Number","Paired","Pvalue",
                                  "Enrich_median", "Enrich_rank",
                                  paste0(pr, "_occurence"), "median_all",
                                  paste0(pr, "_median"), paste0(pr, "_rank"))
  res$Block <- paste0(pr[1], "_vs_", pr[2])
  res.cln <- res %>% select(c(1,14,2:13)) %>%
    mutate(Pvalue=as.numeric(Pvalue)) %>%
    mutate(FDR=p.adjust(Pvalue, method = "BH")) %>%
    arrange(FDR, Pvalue)
  res2 <- res.cln[,c(1:5,15,6:14)]


  # scale profile
  dat.prf.cln <- dat.prf[, -1]
  dat.phe.cln <- dat.phe %>% mutate(Group=ifelse(Stage==pr[1], 0, 1))
  idx <- which(colnames(dat.phe.cln) == "Group")

  # glm result for odd ratios 95%CI
  glmFun <- function(m, n){
    # calculate the glm between profile and group information
    #
    # Args:
    #   m:  result of group information which must to be numeric
    #   n:  taxonomy to be glm
    #
    # Returns:
    #   the glm result of between taxonomy group
    n[n==0] <- min(n[n!=0])
    dat.glm <- data.frame(group=m, marker=log(n)) %>% na.omit()
    model <- summary(rlm(group ~ marker, data = dat.glm,
                         family = binomial(link = "logit")))
    res <- signif(exp(model$coefficients["marker",1]) +
                    qnorm(c(0.025,0.5,0.975)) * model$coefficients["marker",1], 2)

    return(res)
  }

  glm_res <- t(apply(dat.prf.cln, 2, function(x, group){
    res <- glmFun(group, as.numeric(x))
    return(res)
  }, group = dat.phe.cln[, idx]))
  Odd <- glm_res %>% data.frame() %>%
    setNames(c("upper", "expected","lower")) %>%
    mutate("Odds Ratio (95% CI)" = paste0(expected, " (", lower, ";", upper, ")"))
  Odd$type <- rownames(glm_res)

  res_merge <- inner_join(res2,
                          Odd[, c(4:5)], by = "type")

  return(res_merge)
}
