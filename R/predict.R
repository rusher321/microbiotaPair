#' LassoFeaturePlot
#' plot the no-zeor lasso score of feature
#' @param lassoresult
#' @param top
#' @param type
#'
#' @return
#' @export
#'
#' @examples
LassoFeatruePlot <- function(lassoresult ,top=NULL, type, title, optimal){

    # feature importance
    featurescore <- lassoresult$feature
    # get the max & min
    minmax <- c()
    for(i in 1:length(featurescore)){
      tmpminmax <- c(max(featurescore[[i]][,2]), min(featurescore[[i]][,2]))
      minmax <- c(minmax, tmpminmax)
    }
    xmin <- min(minmax)
    xmax <- max(minmax)

    # plot feature impotance

    qplot <- function(qdat, top= NULL){

      qdat <- qdat[order(qdat$Score), ]
      if(is.null(top)){
        order_Type <- as.character(qdat$Type)
        qdat$Type <- factor(qdat$Type, levels = order_Type)
      }else{
        qdat <- tail(qdat, top)
        qdat <- arrange(qdat, Score)
        order_Type <- as.character(qdat$Type)
        qdat$Type <- factor(qdat$Type, levels = order_Type)
      }

      qdat$dir <- factor(qdat$dir, levels = c("pos", "neg"))
      p <- ggplot(qdat, aes(y=Type,x=Score,color=dir)) +
        scale_x_continuous(limits =c(xmin,xmax)) +
        geom_segment(xend=0,aes(yend=Type),size=5)  +
        geom_vline(xintercept = 0) +
        scale_color_manual(values = c("#ef3b2c", "#2171b5"))+
        xlab("Coefficient estimate")+ylab("")+ggtitle(title)+mytheme


      return(p)
    }

    qlist <- list()
    for(i in 1:length(featurescore)){
      qlist[[i]] <- qplot(featurescore[[i]])
    }

    # predict performance
    model <- lassoresult$model
    predict <- predict(lassoresult$model, lassoresult$x, type = type, s=optimal)
    qdat <- data.frame(Ture = lassoresult$y, Predict = predict[,1])

    if(type == "auc" | type=="response"){
       p2 <- rocPlot(response = as.factor(qdat$Ture), predict = qdat$Predict)
    }else{
      text <- round(cor(qdat[,1], qdat[,2], method = "s"), 2)
      grob <- grobTree(textGrob(paste0(title,"Spearman = ", text), x=0.1,  y=0.95, hjust=0,
          gp=gpar(col="red", fontsize=13, fontface="italic")))

       p2 <- ggplot(qdat, aes(x = Ture, y = Predict))+geom_point()+xlab("Actual Value(scale)")+
         ylab("Predict Value(scale)")+
        annotation_custom(grob)
    }

  return(list(qlist, p2))

}


#' CvLasso
#' select feature & predict performance using lasso method
#'
#' @param responsetag
#' @param response
#' @param microbio
#' @param transformM
#' @param normalization
#' @param family
#' @param numfold
#' @param ...
#'
#' @return
#' @import glmnet
#' @import dplyr
#' @import tibble
#' @export
#'
#' @examples
CvLasso <- function(responsetag, response, microbio, transformM = "IQN",
                            normalization=NULL, family, numfold, optimal,...){
  # notice the type parameter depend on your family
  # match the sample ID
  id <- intersect(rownames(response), rownames(microbio))
  print(paste0("the sample size is ", length(id)))

  # to transform the data
  tmp <- microbio[id, ]
  invt <- function(x){qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
  AST <- function(x) {return(sign(x) * asin(sqrt(abs(x))))}

  if(transformM == "IQN"){
    mdat <- as.data.frame(apply(tmp, 2, invt))
  }else if(transformM == "AST"){
    mdat <- as.data.frame(apply(tmp, 2, AST))
  }else{
    mdat <- tmp
  }

  if(is.null(normalization)){
    mdat2 <- data.frame(cbind(response[id, responsetag, drop =F] , mdat))
  }else if(normalization=="IQN"){
    mdat2 <- data.frame(tmpvar <- invt(response[id, responsetag, drop =F]) , mdat)
  }else if(normalization=="scale"){
    mdat2 <- data.frame(tmpvar <- scale(response[id, responsetag, drop =F]) , mdat)
  }

  # rm the na row
  mdat2 <- mdat2[!is.na(mdat2[,1]), ]

  # lasso model
  if(nrow(mdat2) < 10){
    stop("the sample number is below 10, please check the number!")
  }else{
    # model
    set.seed(123)
    if(family == "binomial"){
      lasso <- cv.glmnet(x=as.matrix(mdat2[,-1]),
                         y=as.factor(mdat2[,1]),
                         #family='gaussian',
                         family = family,
                         nfolds = numfold,
                         alpha = 1,
                         nlambda = 100,
                         keep = T,...)
    }else{
      lasso <- cv.glmnet(x=as.matrix(mdat2[,-1]),
                       y=mdat2[,1],
                       #family='gaussian',
                       family = family,
                       nfolds = numfold,
                       alpha = 1,
                       nlambda = 100,
                       keep = T,...)
    }

    # feature importance
    featurescore <- list()
    if(optimal=="1se"){
      lasso.mk <- coef(lasso, lasso$lambda.1se)
    }else{
      lasso.mk <- coef(lasso, lasso$lambda.min)
    }
    if(family != "multinomial"){

      featurescore[[1]] <- data.frame(as.matrix(lasso.mk)) %>%
      setNames("Score") %>%
      rownames_to_column("Type") %>%
      slice(-c(1:2)) %>%
      filter(Score!=0) %>%
      mutate(dir = ifelse(Score >0 ,"pos", "neg"))%>%
      arrange(abs(Score))

    }else if(family == "multinomial"){


      for(i in 1:length(lasso.mk)){
        featurescore[[i]] <- data.frame(as.matrix(lasso.mk[[i]])) %>%
            setNames("Score") %>%
            rownames_to_column("Type") %>%
            slice(-c(1:2)) %>%
            filter(Score!=0) %>% mutate(dir = ifelse(Score >0 ,"pos", "neg"))%>%
            arrange(abs(Score))

        }
    }

    out <- list(model = lasso, x=as.matrix(mdat2[,-1]),
                y=mdat2[,1], feature=featurescore)
    return(out)

  }
}

##  random foreset

#' randomForestTwo
#' two class randomforest model
#' @param data  data.frame or matrix,species data or metabolism data, row is sample ID
#' @param metadata data.frame or matrix, 0/1 response dataset, row is sample ID
#' @param response character, response feature name
#' @param repeatNum numberic, repeat number of feeature selection produce
#' @param foldNum numberic, cross validation
#' @param factorLev vector, 0/1 factor name
#'
#' @return list
#' @import caret
#' @export
#'
#' @examples
randomForestTwo <- function(data, metadata, response, repeatNum, foldNum, factorLev){


#  feature selection function

features <- function(trainx, trainy, foldNum, repeatNum){
  # make sure data & config 's sample is same

  index <- nearZeroVar(trainx)  # find the zero-variance var
  if(length(index) == 0){data.cle <- trainx}else{data.cle <- trainx[, -index]}

  # if feature over 100, only add 10 repeats compution
  # if the sample number is small, set the feature number not higher the sample number

  numsample <- nrow(trainx)
  if(ncol(data.cle) >= 100 & numsample>=100){
    subsets <- c(1:100,seq(101, ncol(data.cle), by=ceiling(ncol(data.cle)-100)/10))
  }else if(numsample<100 & ncol(data.cle) > numsample){
    subsets <- c(1:numsample)
  }else{
    subsets <- c(1:ncol(data.cle))
  }

  if(foldNum == "leaveone"){
　　ctrl = rfeControl(functions = rfFuncs,
　　                      #  rerank = TRUE, # iteration produce to recaluate the feature importance
　　                      method = "cv", number = numsample, repeats = repeatNum,
　　                      verbose = FALSE, returnResamp = "all", saveDetails = T,
　　                      allowParallel=T)
  }else{
    ctrl = rfeControl(functions = rfFuncs,
                      #  rerank = TRUE, # iteration produce to recaluate the feature importance
                      method = "cv", number = foldNum, repeats = repeatNum,
                      verbose = FALSE, returnResamp = "all", saveDetails = T,
                      allowParallel=T)
  }

  Profile = rfe(data.cle, trainy, sizes = subsets, rfeControl = ctrl)
  return(Profile)
}

# model produce

model <- function(trainx, trainy, foldNum){
  # set the mtry
  if(ncol(trainx) <= 4 || is.null(ncol(trainx))){
    grid <- expand.grid(mtry = 2)
  }else{
    grid <- expand.grid(mtry = seq(2,floor(sqrt(ncol(trainx))), by=1))
  }

  if(foldNum == "leaveone"){
    ctrl <- trainControl(method="LOOCV",
                         summaryFunction=twoClassSummary,  # Use AUC to pick the best model
                         savePredictions=T,
                         classProbs=TRUE,
                         sampling = "up")

  }else{
    ctrl <- trainControl(method="repeatedcv", number = foldNum, repeats = 10,
                       summaryFunction=twoClassSummary,  # Use AUC to pick the best model
                       savePredictions=T,
                       classProbs=TRUE,
                       sampling = "up")
  }

  rf  <- train(trainx, trainy,
               method = "rf",   # Radial kernel
               tuneGrid = grid, #
               trControl=ctrl)

  return(rf)
}

# wilcox test to get direction

testdir <- function(data){
  dir <- c()
  for(i in 1:(ncol(data)-1)){
    x <- data[,i]
    y <- data$y
    meanr <- tapply(rank(x), y, mean)
    dir[i] <- ifelse(meanr[1] > meanr[2], "low", "high")
    pvalue <- wilcox.test(x~y)$p.value
    dir[i] <- ifelse(pvalue <0.05, dir[i], "none")
  }
  return(dir)
}

# ready data

 id <- intersect(rownames(data), rownames(metadata))
 trainy <- metadata[id , response]
 trainx <- data[id, ]
 naindex <- which(is.na(trainy)) # rm the NA response
 if(length(naindex)!=0){
    trainx <- trainx[-naindex, ]
    trainy <- trainy[-naindex]
  }
 trainy <- as.factor(trainy)
 levels(trainy) <- factorLev

# feature selection , can impove the repeatNum to makesure get more rubust the feature
featureRes <- features(trainx = trainx, trainy = trainy, foldNum = foldNum, repeatNum = repeatNum)
opt <- featureRes$optVariables

mod <- model(trainx = trainx[,opt,drop=F], trainy = trainy, foldNum = foldNum)


feature <- as.data.frame(mod$finalModel$importance)
data <- as.data.frame(trainx[,rownames(feature), drop=F])
data$y <- trainy
dir <- testdir(data)
feature$dir <- dir

out <- list(varSel = featureRes, feature = feature ,  model = mod)

return(out)

}


# model result plot

#' perfPlot
#' performance of randomForesttwo
#' @param randomRes obeject of randomForesttwo
#' @param title charcter, title of figure
#'
#' @return
#' @import pROC
#' @import plotROC
#' @import ggpubr
#' @import ggplot2
#' @export
#'
#' @examples
perfPlot <- function(randomRes, title){

  model <- randomRes$model
  besttune <- as.numeric(model$bestTune)
  resample <- model$pred[model$pred$mtry == besttune, ]
  # roc plot
  roc <- rocPlot2(resample$obs, resample[,4], title=title)
  # feature plot
  feature <- randomRes$feature
  feature$feature <- rownames(feature)
  featureim <- ggpubr::ggbarplot(feature, x = "feature", y = "MeanDecreaseGini",
            fill = "dir",
            color = "white",
            palette = "jco",
            sort.val = "desc",
            sort.by.groups = FALSE,
            x.text.angle = 90)

  out <- list(roc, featureim)
  return(out)

}



<<<<<<< Updated upstream
randomForestTworobust <- function(data, metadata, response, repeatNum, foldNum, factorLev){


  #  feature selection function

  features <- function(trainx, trainy, repeatNum){
    # make sure data & config 's sample is same

    index <- nearZeroVar(trainx)  # find the zero-variance var
    if(length(index) == 0){data.cle <- trainx}else{data.cle <- trainx[, -index]}
    # if feature over 100, only add 10 repeats compution
    if(ncol(data.cle) >= 100){
      subsets <- c(1:100,seq(101, ncol(data.cle), by=ceiling(ncol(data.cle)-100)/10))
    }else{
      subsets <- c(1:ncol(data.cle))
    }
    ctrl = rfeControl(functions = rfFuncs, method = "cv", repeats = repeatNum,
                      verbose = FALSE, returnResamp = "final", saveDetails = T,
                      allowParallel=T)
    Profile = rfe(data.cle, trainy, sizes = subsets, rfeControl = ctrl)
    return(Profile)
  }


  model <- function(trainx, trainy, foldNum){
    # set the mtry
    if(ncol(trainx) <= 4 || is.null(ncol(trainx))){
      grid <- expand.grid(mtry = 2)
    }else{
      grid <- expand.grid(mtry = seq(2,floor(sqrt(ncol(trainx))), by=1))
    }
    ctrl <- trainControl(method="repeatedcv", number = foldNum, repeats = 10,
                         summaryFunction=twoClassSummary,  # Use AUC to pick the best model
                         savePredictions=T,
                         classProbs=TRUE,
                         sampling = "up")

    rf  <- train(trainx, trainy,
                 method = "rf",   # Radial kernel
                 tuneGrid = grid, #
                 trControl=ctrl)

    return(rf)
  }

  # ready data
  id <- intersect(rownames(data), rownames(metadata))
  trainy <- metadata[id , response]
  trainx <- data[id, ]
  naindex <- which(is.na(trainy)) # rm the NA response
  if(length(naindex)!=0){
    trainx <- trainx[-naindex, ]
    trainy <- trainy[-naindex]
  }
  trainy <- as.factor(trainy)
  levels(trainy) <- factorLev

  # split the data & get the cv dataset
  if(foldNum == "leaveone"){
    foldNum <- length(trainy)
  }
  foldlist <- createFolds(trainy, foldNum)

  # output
  featurelist <- list()
  pred <- c()
  ture <- c()

  for(i in 1:foldNum){
    trainxsub <- trainx[-foldlist[[i]], ]
    trainysub <- trainy[-foldlist[[i]]]
    testxsub <- trainx[foldlist[[i]], ]
    testysub <- trainy[foldlist[[i]]]

    featureRes <- features(trainx = trainxsub, trainy = trainysub, repeatNum = repeatNum)
    opt <- featureRes$optVariables
    mod <- model(trainx = trainxsub[,opt, drop=F], trainy = trainysub, foldNum = 5)

    featurelist[[i]] <- mod$finalModel$importance
    pred <- c(pred, predict(mod, testxsub[, opt, drop=F],type="prob")[,2])
    ture <- c(ture, testysub)

  }

  preddat <- data.frame(predV = pred, obsV = ture)
  out <- list(feature = featurelist,  cvres = preddat)

  return(out)

}

  randomForestTworobust <- function(data, metadata, response, repeatNum, foldNum, factorLev){


  #  feature selection function

  features <- function(trainx, trainy, repeatNum){
    # make sure data & config 's sample is same

    index <- nearZeroVar(trainx)  # find the zero-variance var
    if(length(index) == 0){data.cle <- trainx}else{data.cle <- trainx[, -index]}
    # if feature over 100, only add 10 repeats compution
    if(ncol(data.cle) >= 100){
      subsets <- c(1:100,seq(101, ncol(data.cle), by=ceiling(ncol(data.cle)-100)/10))
    }else{
      subsets <- c(1:ncol(data.cle))
    }
    ctrl = rfeControl(functions = rfFuncs, method = "cv", repeats = repeatNum,
                      verbose = FALSE, returnResamp = "final", saveDetails = T,
                      allowParallel=T)
    Profile = rfe(data.cle, trainy, sizes = subsets, rfeControl = ctrl)
    return(Profile)
  }


  model <- function(trainx, trainy, foldNum){
    # set the mtry
    if(ncol(trainx) <= 4 || is.null(ncol(trainx))){
      grid <- expand.grid(mtry = 2)
    }else{
      grid <- expand.grid(mtry = seq(2,floor(sqrt(ncol(trainx))), by=1))
    }
    ctrl <- trainControl(method="repeatedcv", number = foldNum, repeats = 10,
                         summaryFunction=twoClassSummary,  # Use AUC to pick the best model
                         savePredictions=T,
                         classProbs=TRUE,
                         sampling = "up")

    rf  <- train(trainx, trainy,
                 method = "rf",   # Radial kernel
                 tuneGrid = grid, #
                 trControl=ctrl)

    return(rf)
  }

  # ready data
  id <- intersect(rownames(data), rownames(metadata))
  trainy <- metadata[id , response]
  trainx <- data[id, ]
  naindex <- which(is.na(trainy)) # rm the NA response
  if(length(naindex)!=0){
    trainx <- trainx[-naindex, ]
    trainy <- trainy[-naindex]
  }
  trainy <- as.factor(trainy)
  levels(trainy) <- factorLev

  # split the data & get the cv dataset
  if(foldNum == "leaveone"){
    foldNum <- length(trainy)
  }
  foldlist <- createFolds(trainy, foldNum)

  # output
  featurelist <- list()
  pred <- c()
  ture <- c()

  for(i in 1:foldNum){
    trainxsub <- trainx[-foldlist[[i]], ]
    trainysub <- trainy[-foldlist[[i]]]
    testxsub <- trainx[foldlist[[i]], ]
    testysub <- trainy[foldlist[[i]]]

    featureRes <- features(trainx = trainxsub, trainy = trainysub, repeatNum = repeatNum)
    opt <- featureRes$optVariables
    mod <- model(trainx = trainxsub[,opt, drop=F], trainy = trainysub, foldNum = 5)

    featurelist[[i]] <- mod$finalModel$importance
    pred <- c(pred, predict(mod, testxsub[, opt, drop=F],type="prob")[,2])
    ture <- c(ture, testysub)

  }

  preddat <- data.frame(predV = pred, obsV = ture)
  out <- list(feature = featurelist,  cvres = preddat)

  return(out)

}

=======
# randomForestTworobust <- function(data, metadata, response, repeatNum, foldNum, factorLev){
#
#
#   #  feature selection function
#
#   features <- function(trainx, trainy, repeatNum){
#     # make sure data & config 's sample is same
#
#     index <- nearZeroVar(trainx)  # find the zero-variance var
#     if(length(index) == 0){data.cle <- trainx}else{data.cle <- trainx[, -index]}
#     # if feature over 100, only add 10 repeats compution
#     if(ncol(data.cle) >= 100){
#       subsets <- c(1:100,seq(101, ncol(data.cle), by=ceiling(ncol(data.cle)-100)/10))
#     }else{
#       subsets <- c(1:ncol(data.cle))
#     }
#     ctrl = rfeControl(functions = rfFuncs, method = "cv", repeats = repeatNum,
#                       verbose = FALSE, returnResamp = "final", saveDetails = T,
#                       allowParallel=T)
#     Profile = rfe(data.cle, trainy, sizes = subsets, rfeControl = ctrl)
#     return(Profile)
#   }
#
#
#   model <- function(trainx, trainy, foldNum){
#     # set the mtry
#     if(ncol(trainx) <= 4 || is.null(ncol(trainx))){
#       grid <- expand.grid(mtry = 2)
#     }else{
#       grid <- expand.grid(mtry = seq(2,floor(sqrt(ncol(trainx))), by=1))
#     }
#     ctrl <- trainControl(method="repeatedcv", number = foldNum, repeats = 10,
#                          summaryFunction=twoClassSummary,  # Use AUC to pick the best model
#                          savePredictions=T,
#                          classProbs=TRUE,
#                          sampling = "up")
#
#     rf  <- train(trainx, trainy,
#                  method = "rf",   # Radial kernel
#                  tuneGrid = grid, #
#                  trControl=ctrl)
#
#     return(rf)
#   }
#
#   # ready data
#   id <- intersect(rownames(data), rownames(metadata))
#   trainy <- metadata[id , response]
#   trainx <- data[id, ]
#   naindex <- which(is.na(trainy)) # rm the NA response
#   if(length(naindex)!=0){
#     trainx <- trainx[-naindex, ]
#     trainy <- trainy[-naindex]
#   }
#   trainy <- as.factor(trainy)
#   levels(trainy) <- factorLev
#
#   # split the data & get the cv dataset
#   if(foldNum == "leaveone"){
#     foldNum <- length(trainy)
#   }
#   foldlist <- createFolds(trainy, foldNum)
#
#   # output
#   featurelist <- list()
#   pred <- c()
#   ture <- c()
#
#   for(i in 1:foldNum){
#     trainxsub <- trainx[-foldlist[[i]], ]
#     trainysub <- trainy[-foldlist[[i]]]
#     testxsub <- trainx[foldlist[[i]], ]
#     testysub <- trainy[foldlist[[i]]]
#
#     featureRes <- features(trainx = trainxsub, trainy = trainysub, repeatNum = repeatNum)
#     opt <- featureRes$optVariables
#     mod <- model(trainx = trainxsub[,opt, drop=F], trainy = trainysub, foldNum = 5)
#
#     featurelist[[i]] <- mod$finalModel$importance
#     pred <- c(pred, predict(mod, testxsub[, opt, drop=F],type="prob")[,2])
#     ture <- c(ture, testysub)
#
#   }
#
#   preddat <- data.frame(predV = pred, obsV = ture)
#   out <- list(feature = featurelist,  cvres = preddat)
#
#   return(out)
#
# }randomForestTworobust <- function(data, metadata, response, repeatNum, foldNum, factorLev){
#
#
#   #  feature selection function
#
#   features <- function(trainx, trainy, repeatNum){
#     # make sure data & config 's sample is same
#
#     index <- nearZeroVar(trainx)  # find the zero-variance var
#     if(length(index) == 0){data.cle <- trainx}else{data.cle <- trainx[, -index]}
#     # if feature over 100, only add 10 repeats compution
#     if(ncol(data.cle) >= 100){
#       subsets <- c(1:100,seq(101, ncol(data.cle), by=ceiling(ncol(data.cle)-100)/10))
#     }else{
#       subsets <- c(1:ncol(data.cle))
#     }
#     ctrl = rfeControl(functions = rfFuncs, method = "cv", repeats = repeatNum,
#                       verbose = FALSE, returnResamp = "final", saveDetails = T,
#                       allowParallel=T)
#     Profile = rfe(data.cle, trainy, sizes = subsets, rfeControl = ctrl)
#     return(Profile)
#   }
#
#
#   model <- function(trainx, trainy, foldNum){
#     # set the mtry
#     if(ncol(trainx) <= 4 || is.null(ncol(trainx))){
#       grid <- expand.grid(mtry = 2)
#     }else{
#       grid <- expand.grid(mtry = seq(2,floor(sqrt(ncol(trainx))), by=1))
#     }
#     ctrl <- trainControl(method="repeatedcv", number = foldNum, repeats = 10,
#                          summaryFunction=twoClassSummary,  # Use AUC to pick the best model
#                          savePredictions=T,
#                          classProbs=TRUE,
#                          sampling = "up")
#
#     rf  <- train(trainx, trainy,
#                  method = "rf",   # Radial kernel
#                  tuneGrid = grid, #
#                  trControl=ctrl)
#
#     return(rf)
#   }
#
#   # ready data
#   id <- intersect(rownames(data), rownames(metadata))
#   trainy <- metadata[id , response]
#   trainx <- data[id, ]
#   naindex <- which(is.na(trainy)) # rm the NA response
#   if(length(naindex)!=0){
#     trainx <- trainx[-naindex, ]
#     trainy <- trainy[-naindex]
#   }
#   trainy <- as.factor(trainy)
#   levels(trainy) <- factorLev
#
#   # split the data & get the cv dataset
#   if(foldNum == "leaveone"){
#     foldNum <- length(trainy)
#   }
#   foldlist <- createFolds(trainy, foldNum)
#
#   # output
#   featurelist <- list()
#   pred <- c()
#   ture <- c()
#
#   for(i in 1:foldNum){
#     trainxsub <- trainx[-foldlist[[i]], ]
#     trainysub <- trainy[-foldlist[[i]]]
#     testxsub <- trainx[foldlist[[i]], ]
#     testysub <- trainy[foldlist[[i]]]
#
#     featureRes <- features(trainx = trainxsub, trainy = trainysub, repeatNum = repeatNum)
#     opt <- featureRes$optVariables
#     mod <- model(trainx = trainxsub[,opt, drop=F], trainy = trainysub, foldNum = 5)
#
#     featurelist[[i]] <- mod$finalModel$importance
#     pred <- c(pred, predict(mod, testxsub[, opt, drop=F],type="prob")[,2])
#     ture <- c(ture, testysub)
#
#   }
#
#   preddat <- data.frame(predV = pred, obsV = ture)
#   out <- list(feature = featurelist,  cvres = preddat)
#
#   return(out)
#
# }
#
>>>>>>> Stashed changes

#### 
   randomForesttwo <- function(x){
   
   # feature selection 
   # rm the zero 
   # select the lowest wilcox pvalue,  & rmmove the species with the highest correlation value  >  cutoff  (o.6 )
   # then the second run , until get the top X highest feature 
     
       
   # model on leave one 
        
   
   
   
   
   }



