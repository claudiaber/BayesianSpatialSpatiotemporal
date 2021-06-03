ForecastandMeasure <- function(future_dist, x1_fore, x2_fore, Y_fore,
                               fitglm, fitcarstan, fitSTcar, fitSTtrend,
                               fitstan=NULL,
                               fitSTcarfix=NULL,
                               W,Y_matrix,
                               PredictedTime,
                               K=length(x1_fore),
                               Perc_REC,#= c(0.03, 0.05, 0.1, 0.5),
                               Thresh_REC ){#= 0.5){
  # future_dist <- 4# 2018 09 --> 2019 01
  newdata <- data.frame("x1_matrix"=x1_fore, 
                        "x2_matrix"=x2_fore)
  predictionGLM <- predict(fitglm,type = "response",newdata = newdata)
  
  Theta_CARStan <- Theta_STcar <- Theta_STtrend <- Theta_glm <- Theta_STfix <-  NULL
  # for(i in 1:(ncol(x1_matrix))){
  # print(i)
  # print("GLM")
  if(!is.null(fitstan)){
    print("GLMstan")
    extglmstan <- extract(fitstan)
    Logit_singletime <-  matrix(extglmstan$b_0,
                                nrow = length(extglmstan$b_0),
                                ncol = K, byrow = FALSE) +
      extglmstan$b_x1 %*% t(x1_fore) + 
      extglmstan$b_x2 %*% t(x2_fore)
    eLogit <- exp(Logit_singletime)
    rm(Logit_singletime)
    Theta_samples_Stan_singletime <- eLogit / (1 + eLogit)
    rm(eLogit)
    Theta_glm <- rbind(Theta_glm, cbind(
      Mean=apply(Theta_samples_Stan_singletime,2,mean),
      t(apply(Theta_samples_Stan_singletime,2,quantile, 
              probs=c(0.025,0.975)))
    ))
    rm(fitstan, extglmstan, Theta_samples_Stan_singletime)
    gc()
  }
  
  print("CAR")
  extcarstan <- extract(fitcarstan)
  Logit_singletime <-  matrix(extcarstan$b_0,
                              nrow = length(extcarstan$b_0),
                              ncol = K, byrow = FALSE) +
    extcarstan$b_x1 %*% t(x1_fore) + 
    extcarstan$b_x2 %*% t(x2_fore)+ 
    extcarstan$phi 
  eLogit <- exp(Logit_singletime)
  rm(Logit_singletime)
  Theta_samples_Stan_singletime <- eLogit / (1 + eLogit)
  rm(eLogit)
  Theta_CARStan <- rbind(Theta_CARStan, cbind(
    Mean=apply(Theta_samples_Stan_singletime,2,mean),
    t(apply(Theta_samples_Stan_singletime,2,quantile, 
            probs=c(0.025,0.975)))
  ))
  rm(fitcarstan, extcarstan, Theta_samples_Stan_singletime)
  gc()
  
  
  if(!is.null(fitSTcarfix)){
    print("STfix")
    extSTfix <- extract(fitSTcarfix)
    Logit_singletime <-  matrix(extSTfix$b_0,
                                nrow = length(extSTfix$b_0),
                                ncol = K, byrow = FALSE) +
      extSTfix$b_x1 %*% t(x1_fore) + 
      extSTfix$b_x2 %*% t(x2_fore)+ 
      extSTfix$phi
    eLogit <- exp(Logit_singletime)
    rm(Logit_singletime)
    Theta_samples_Stan_singletime <-  eLogit / (1 + eLogit)
    rm(eLogit)
    Theta_STfix <- rbind(Theta_STfix, cbind(
      Mean=apply(Theta_samples_Stan_singletime,2,mean),
      t(apply(Theta_samples_Stan_singletime,2,quantile, 
              probs=c(0.025,0.975)))
    ))
    rm(extSTfix,fitSTcarfix,  Theta_samples_Stan_singletime)
    gc()
    
  }
    
  print("AST")
  extSTcar <- extract(fitSTcar)
  Chi_future <- extSTcar$ar1_par^future_dist *extSTcar$chi[,ncol(extSTcar$chi)]
  Logit_singletime <-  matrix(extSTcar$b_0,
                              nrow = length(extSTcar$b_0),
                              ncol = K, byrow = FALSE) +
    extSTcar$b_x1 %*% t(x1_fore) + 
    extSTcar$b_x2 %*% t(x2_fore)+ 
    extSTcar$phi + 
    matrix(Chi_future, 
           ncol=K,
           nrow=length(extSTcar$chi[,ncol(extSTcar$chi)]))
  eLogit <- exp(Logit_singletime)
  rm(Logit_singletime)
  Theta_samples_Stan_singletime <-  eLogit / (1 + eLogit)
  rm(eLogit)
  Theta_STcar <- rbind(Theta_STcar, cbind(
    Mean=apply(Theta_samples_Stan_singletime,2,mean),
    t(apply(Theta_samples_Stan_singletime,2,quantile, 
            probs=c(0.025,0.975)))
  ))
  rm(extSTcar,fitSTcar,  Theta_samples_Stan_singletime)
  gc()
  
  print("ASTtrend")
  extSTtrend <- extract(fitSTtrend)
  Chi_future <- extSTtrend$ar1_par^future_dist *extSTtrend$chi[,ncol(extSTtrend$chi)]
  Logit_singletime <-  matrix(extSTtrend$b_0,
                              nrow = length(extSTtrend$b_0),
                              ncol = K, byrow = FALSE) +
    extSTtrend$b_x1 %*% t(x1_fore) + 
    extSTtrend$b_x2 %*% t(x2_fore)+ 
    matrix( extSTtrend$b_3 *PredictedTime,
            ncol=K,
            nrow=length(extSTtrend$b_3)) + 
    extSTtrend$phi + 
    matrix(Chi_future, 
           ncol=K,
           nrow=length(extSTtrend$chi[,ncol(extSTtrend$chi)]))
  eLogit <- exp(Logit_singletime)
  rm(Logit_singletime)
  Theta_samples_Stan_singletime <- eLogit / (1 + eLogit)
  rm(eLogit)
  Theta_STtrend <- rbind(Theta_STtrend, cbind(
    Mean=apply(Theta_samples_Stan_singletime,2,mean),
    t(apply(Theta_samples_Stan_singletime,2,quantile, 
            probs=c(0.025,0.975)))
  ))
  rm(extSTtrend, fitSTtrend, Theta_samples_Stan_singletime)
  gc()
  
  #   
  # }
  REC <- AUC <- NULL
  REC_bayes <- AUC_bayes <- NULL
  REC_CAR <- AUC_CAR <- NULL
  REC_fix <- AUC_fix <- NULL
  REC_ST <- AUC_ST <- NULL
  REC_STtrend <- AUC_STTrend <- NULL
  
  
  print("GLM")
  pred <- ROCR::prediction( predictions = predictionGLM, labels =Y_fore )
  AUC <- c(AUC , ROCR::performance( pred, "auc")@y.values[[1]])
  REC <-c(unlist(lapply(Perc_REC, FUN = function(j) {sum(Y_fore[order(predictionGLM,
                                                                      decreasing = TRUE)[1:ceiling(K*j)]])/
      sum(Y_fore)})),
      sum(Y_fore[predictionGLM>Thresh_REC])/sum(Y_fore),
      sum(Y_fore[predictionGLM>Thresh_REC])/sum(predictionGLM>Thresh_REC)
  )
  if(!is.null(Theta_glm)){
    print("GLMstan")
    predstan <- ROCR::prediction( predictions = Theta_glm[,"Mean"], labels =Y_fore )
    AUC_bayes <- c(AUC_bayes , ROCR::performance( predstan, "auc")@y.values[[1]])
    REC_bayes <-c(unlist(lapply(Perc_REC, FUN = function(j) {sum(Y_fore[order(Theta_glm[,"Mean"],
                                                                              decreasing = TRUE)[1:ceiling(K*j)]])/
        sum(Y_fore)})),
        sum(Y_fore[Theta_glm[,"Mean"]>Thresh_REC])/sum(Y_fore),
        sum(Y_fore[Theta_glm[,"Mean"]>Thresh_REC])/sum(Theta_glm[,"Mean"]>Thresh_REC) )
  }
  
  print("CAR")
  pred1 <- ROCR::prediction(predictions = Theta_CARStan[, "Mean"], labels = Y_fore)
  AUC_CAR <-c(AUC_CAR,  ROCR::performance( pred1, "auc")@y.values[[1]])
  REC_CAR <-  c(unlist(lapply(Perc_REC, FUN = function(j) {sum(Y_fore[order(Theta_CARStan[,"Mean"],
                                                                            decreasing = TRUE)[1:ceiling(K*j)]])/
      sum(Y_fore)})),
      sum(Y_fore[Theta_CARStan[,"Mean"]>Thresh_REC])/sum(Y_fore),
      sum(Y_fore[Theta_CARStan[,"Mean"]>Thresh_REC])/sum(Theta_CARStan[,"Mean"]>Thresh_REC)  )
  
  
  if(!is.null(Theta_STfix)){
    print("STfix")
    pred2 <- ROCR::prediction(predictions = Theta_STfix[, "Mean"], labels = Y_fore)
    AUC_fix <-c(AUC_ST,  ROCR::performance( pred2, "auc")@y.values[[1]])
    REC_fix <- c(unlist(lapply(Perc_REC, FUN = function(j) {sum(Y_fore[order(Theta_STfix[,"Mean"],
                                                                            decreasing = TRUE)[1:ceiling(K*j)]])/
        sum(Y_fore)})),
        sum(Y_fore[Theta_STfix[,"Mean"]>Thresh_REC])/sum(Y_fore),
        sum(Y_fore[Theta_STfix[,"Mean"]>Thresh_REC])/sum(Theta_STfix[,"Mean"]>Thresh_REC)  )  
  }

  
  print("AST")
  pred2 <- ROCR::prediction(predictions = Theta_STcar[, "Mean"], labels = Y_fore)
  AUC_ST <-c(AUC_ST,  ROCR::performance( pred2, "auc")@y.values[[1]])
  REC_ST <- c(unlist(lapply(Perc_REC, FUN = function(j) {sum(Y_fore[order(Theta_STcar[,"Mean"],
                                                                          decreasing = TRUE)[1:ceiling(K*j)]])/
      sum(Y_fore)})),
      sum(Y_fore[Theta_STcar[,"Mean"]>Thresh_REC])/sum(Y_fore),
      sum(Y_fore[Theta_STcar[,"Mean"]>Thresh_REC])/sum(Theta_STcar[,"Mean"]>Thresh_REC)  )
  
  print("ASTtrend")
  pred3 <- ROCR::prediction(predictions = Theta_STtrend[, "Mean"], labels = Y_fore)
  AUC_STTrend <-c(AUC_STTrend,  ROCR::performance( pred3, "auc")@y.values[[1]])
  REC_STtrend <-  c(unlist(lapply(Perc_REC, FUN = function(j) {sum(Y_fore[order(Theta_STtrend[,"Mean"],
                                                                                decreasing = TRUE)[1:ceiling(K*j)]])/
      sum(Y_fore)})),
      sum(Y_fore[Theta_STtrend[,"Mean"]>Thresh_REC])/sum(Y_fore),
      sum(Y_fore[Theta_STtrend[,"Mean"]>Thresh_REC])/sum(Theta_STtrend[,"Mean"]>Thresh_REC)  )
  
  print("NetNaive")
  predNN <- W %*% Y_matrix[,ncol(Y_matrix)]
  Futuretime <- future_dist+ncol(Y_matrix)
  
  
  AUC_NNaive <- NULL
  REC_NNaive <- rep(0,times=length(Perc_REC)+2)
  
  # for(i in 4:10){
  # print(i)
  
  pred3 <- ROCR::prediction(predictions =predNN, labels = Y_fore)
  # perf3 <- ROCR::performance( pred3, "tpr", "fpr" )
  AUC_NNaive <-c(AUC_NNaive,  ROCR::performance( pred3, "auc")@y.values[[1]])
  for(j in 1:length(Perc_REC)){
    ONES_NN <- order(predNN, decreasing = TRUE)[1:ceiling(K*Perc_REC[j])]
    singlePredNN <- rep(0, K)
    singlePredNN[ONES_NN] <- 1
    REC_NNaive[j] <-sum(Y_fore[order(singlePredNN, decreasing = TRUE)[1:ceiling(K*Perc_REC[j])]])/
      sum(Y_fore)
  }
  REC_NNaive[length(Perc_REC)+1] <- NA
  REC_NNaive[length(Perc_REC)+2] <- NA
  
  # }
  
  
  PerfOutTime <- rbind(c(mean(AUC_CAR), mean(AUC_fix), mean(AUC_ST), mean(AUC_STTrend), mean(AUC_bayes), mean(AUC), mean(AUC_NNaive)),
                       cbind(REC_CAR, REC_fix, REC_ST, REC_STtrend, REC_bayes, REC, REC_NNaive))
  
  rownames(PerfOutTime) <- c("AUC", paste0("$REC\\_{",
                                           Perc_REC,"} (0-",
                                           round(pmin(rowMeans(ceiling(K*Perc_REC) %*% t(1/sum(Y_fore))),
                                                      1)*100,2),")$"),
                             "Recall", "Precision")
  if(!is.null(AUC_fix)){
  colnames(PerfOutTime) <- c("CAR", "Fix","AST", "ASTtrend", "GLMstan","GLM", "NetNaive")}else{
    colnames(PerfOutTime) <- c("CAR","AST", "ASTtrend", "GLMstan","GLM", "NetNaive")
  }
  
  PerfOutTime <- PerfOutTime*100
  if(!is.null(Theta_glm)){
    return(list(predictionGLM=predictionGLM, 
                Theta_glm=Theta_glm,
                Theta_CARStan=Theta_CARStan, 
                Theta_STfix=Theta_STfix,
                Theta_STcar=Theta_STcar,
                Theta_STtrend= Theta_STtrend,
                predNN=predNN,  PerfOutTime=PerfOutTime))
  }else{
    return(list(predictionGLM=predictionGLM, 
                Theta_CARStan=Theta_CARStan,
                Theta_STfix=Theta_STfix,
                Theta_STcar=Theta_STcar,
                Theta_STtrend= Theta_STtrend,
                predNN=predNN,  PerfOutTime=PerfOutTime))
  }
}