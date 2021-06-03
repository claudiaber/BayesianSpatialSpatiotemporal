MeasureFitted <- function(predictionGLM, predNN,
                          Theta_CARStan, Theta_STcar, Theta_STtrend,
                          Theta_STfix=NULL,
                          Y_matrix,
                          K=nrow(Y_matrix),
                          Perc_REC,Thresh_REC){
REC <- AUC <- NULL
REC_CAR <- AUC_CAR <- NULL
REC_fix <- AUC_fix <- NULL
REC_ST <- AUC_ST <- NULL
REC_STtrend <- AUC_STTrend <- NULL
for(i in 1:(ncol(x1_matrix))){
  print(i)
  print("GLM")
  pred <- ROCR::prediction( predictions = predictionGLM[((i-1)*K+1):(i*K)], labels =Y_matrix[,i] )
  AUC <- c(AUC , ROCR::performance( pred, "auc")@y.values[[1]])
  REC <- cbind(REC ,c(unlist(lapply(Perc_REC, FUN = function(j) {
    sum(Y_matrix[order(predictionGLM[((i-1)*K+1):(i*K)],
                       decreasing = TRUE)[1:ceiling(K*j)],i])/
      sum(Y_matrix[,i])})),
    sum(Y_matrix[predictionGLM[((i-1)*K+1):(i*K)]>Thresh_REC,i])/sum(Y_matrix[,i]),
    sum(Y_matrix[predictionGLM[((i-1)*K+1):(i*K)]>Thresh_REC,i])/sum(predictionGLM[((i-1)*K+1):(i*K)]>Thresh_REC)))
  rm(pred)
  
  print("CAR")
  pred1 <- ROCR::prediction(predictions = Theta_CARStan[((i-1)*K+1):(i*K), "Mean"], labels = Y_matrix[,i])
  AUC_CAR <-c(AUC_CAR,  ROCR::performance( pred1, "auc")@y.values[[1]])
  REC_CAR <- cbind(REC_CAR, c(unlist(lapply(Perc_REC, FUN = function(j) {
    sum(Y_matrix[order(Theta_CARStan[((i-1)*K+1):(i*K), "Mean"],
                       decreasing = TRUE)[1:ceiling(K*j)],i])/
      sum(Y_matrix[,i])})),
    sum(Y_matrix[Theta_CARStan[((i-1)*K+1):(i*K), "Mean"]>Thresh_REC,i])/sum(Y_matrix[,i]),
    sum(Y_matrix[Theta_CARStan[((i-1)*K+1):(i*K), "Mean"]>Thresh_REC,i])/
      sum(Theta_CARStan[((i-1)*K+1):(i*K), "Mean"]>Thresh_REC)))
  rm(pred1)
  
  if(!is.null(Theta_STfix)){
  print("FIX")
  pred1 <- ROCR::prediction(predictions = Theta_STfix[((i-1)*K+1):(i*K), "Mean"], labels = Y_matrix[,i])
  AUC_fix <-c(AUC_fix,  ROCR::performance( pred1, "auc")@y.values[[1]])
  REC_fix <- cbind(REC_fix, c(unlist(lapply(Perc_REC, FUN = function(j) {
    sum(Y_matrix[order(Theta_STfix[((i-1)*K+1):(i*K), "Mean"],
                       decreasing = TRUE)[1:ceiling(K*j)],i])/
      sum(Y_matrix[,i])})),
    sum(Y_matrix[Theta_STfix[((i-1)*K+1):(i*K), "Mean"]>Thresh_REC,i])/sum(Y_matrix[,i]),
    sum(Y_matrix[Theta_STfix[((i-1)*K+1):(i*K), "Mean"]>Thresh_REC,i])/
      sum(Theta_STfix[((i-1)*K+1):(i*K), "Mean"]>Thresh_REC)))
  }
  rm(pred1)
  
  print("AST")
  pred2 <- ROCR::prediction(predictions = Theta_STcar[((i-1)*K+1):(i*K), "Mean"], labels = Y_matrix[,i])
  AUC_ST <-c(AUC_ST,  ROCR::performance( pred2, "auc")@y.values[[1]])
  REC_ST <- cbind(REC_ST, c(unlist(lapply(Perc_REC, FUN = function(j) {
    sum(Y_matrix[order(Theta_STcar[((i-1)*K+1):(i*K), "Mean"],
                       decreasing = TRUE)[1:ceiling(K*j)],i])/
      sum(Y_matrix[,i])})),
    sum(Y_matrix[Theta_STcar[((i-1)*K+1):(i*K), "Mean"]>Thresh_REC,i])/sum(Y_matrix[,i]),
    sum(Y_matrix[Theta_STcar[((i-1)*K+1):(i*K), "Mean"]>Thresh_REC,i])/
      sum(Theta_STcar[((i-1)*K+1):(i*K), "Mean"]>Thresh_REC)))
  rm(pred2)
  
  print("ASTtrend")
  pred3 <- ROCR::prediction(predictions = Theta_STtrend[((i-1)*K+1):(i*K), "Mean"], labels = Y_matrix[,i])
  AUC_STTrend <-c(AUC_STTrend,  ROCR::performance( pred3, "auc")@y.values[[1]])
  REC_STtrend <- cbind(REC_STtrend,c(unlist(lapply(Perc_REC, FUN = function(j) {
    sum(Y_matrix[order(Theta_STtrend[((i-1)*K+1):(i*K), "Mean"],
                       decreasing = TRUE)[1:ceiling(K*j)],i])/
      sum(Y_matrix[,i])})),
    sum(Y_matrix[Theta_STtrend[((i-1)*K+1):(i*K), "Mean"]>Thresh_REC,i])/sum(Y_matrix[,i]),
    sum(Y_matrix[Theta_STtrend[((i-1)*K+1):(i*K), "Mean"]>Thresh_REC,i])/
      sum(Theta_STtrend[((i-1)*K+1):(i*K), "Mean"]>Thresh_REC)))
  rm(pred3)
  
}

print("NetNaive")

AUC_NNaive <- NULL
REC_NNaive <- matrix(ncol=n_tmstp-3, nrow=length(Perc_REC)+2)

for(i in 4:ncol(Y_matrix)){
  print(i)
  
  pred3 <- ROCR::prediction(predictions =predNN[,i-3], labels = Y_matrix[,i])
  # perf3 <- ROCR::performance( pred3, "tpr", "fpr" )
  AUC_NNaive <-c(AUC_NNaive,  ROCR::performance( pred3, "auc")@y.values[[1]])
  for(j in 1:length(Perc_REC)){
    ONES_NN <- order(predNN[,i-3], decreasing = TRUE)[1:ceiling(K*Perc_REC[j])]
    singlePredNN <- rep(0, K)
    singlePredNN[ONES_NN] <- 1
    REC_NNaive[j,i-3] <-sum(Y_matrix[order(singlePredNN, decreasing = TRUE)[1:ceiling(K*Perc_REC[j])],i])/
      sum(Y_matrix[,i])
  }
  REC_NNaive[length(Perc_REC)+1,i-3] <- NA
  REC_NNaive[length(Perc_REC)+2,i-3] <- NA
  rm(pred3)
}

PerfInsample <- rbind(c(AUC_CAR[n_tmstp], mean(AUC_fix), mean(AUC_ST), mean(AUC_STTrend), mean(AUC), mean(AUC_NNaive)),
                      cbind(REC_CAR[,n_tmstp], rowMeans(REC_fix), rowMeans(REC_ST), rowMeans(REC_STtrend), rowMeans(REC), rowMeans(REC_NNaive)))
rownames(PerfInsample) <- c("AUC", paste0("$REC\\_{",
                                          Perc_REC,"} (0-",
                                          round(pmin(rowMeans(ceiling(K*Perc_REC) %*% t(1/colSums(Y_matrix))),
                                                     1)*100,2),")$"),
                            "Recall", "Precision")
colnames(PerfInsample) <- c("CAR", "Fix", "AST", "ASTtrend", "GLM", "NetNaive")

PerfInsample <- PerfInsample*100
return(PerfInsample)
}