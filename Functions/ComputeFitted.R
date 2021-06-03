ComputeFitted <- function(fitcarstan,fitSTcarfix=NULL, fitSTcar, fitSTtrend, fitstan=NULL,
                          x1_matrix, x2_matrix, K=nrow(x1_matrix)){
  
  extcarstan <- rstan::extract(fitcarstan)
  extSTcar <- rstan::extract(fitSTcar)
  extSTfix <- rstan::extract(fitSTcarfix)
  extSTtrend <- rstan::extract(fitSTtrend)
  if(!is.null(fitstan)){
    extglmstan <- rstan::extract(fitstan)
  }
  
  Theta_CARStan <- Theta_STcar <- Theta_STtrend <- Theta_glm <- Theta_STfix <- NULL
  for(i in 1:(ncol(x1_matrix))){
    print(i)
    if(!is.null(fitstan)){
      print("GLMstan")
      Logit_singletime <-  matrix(extglmstan$b_0,
                                  nrow = length(extglmstan$b_0),
                                  ncol = K, byrow = FALSE) +
        extglmstan$b_x1 %*% t(x1_matrix[,i]) + 
        extglmstan$b_x2 %*% t(x2_matrix[,i])
      explogit <- exp(Logit_singletime)
      Theta_samples_Stan_singletime <-  explogit / (1 + explogit)
      Theta_glm <- rbind(Theta_glm, cbind(
        Mean=apply(Theta_samples_Stan_singletime,2,mean),
        t(apply(Theta_samples_Stan_singletime,2,quantile, 
                probs=c(0.025,0.975)))
      ))
      
    }
    
    print("CAR")
    Logit_singletime <-  matrix(extcarstan$b_0,
                                nrow = length(extcarstan$b_0),
                                ncol = K, byrow = FALSE) +
      extcarstan$b_x1 %*% t(x1_matrix[,i]) + 
      extcarstan$b_x2 %*% t(x2_matrix[,i])+ 
      extcarstan$phi 
    explogit <- exp(Logit_singletime)
    Theta_samples_Stan_singletime <-  explogit / (1 + explogit)
    Theta_CARStan <- rbind(Theta_CARStan, cbind(
      Mean=apply(Theta_samples_Stan_singletime,2,mean),
      t(apply(Theta_samples_Stan_singletime,2,quantile, 
              probs=c(0.025,0.975)))
    ))
    
    if(!is.null(fitSTcarfix)){
      print("STfix")
      Logit_singletime <-  matrix(extSTfix$b_0,
                                  nrow = length(extSTfix$b_0),
                                  ncol = K, byrow = FALSE) +
        extSTfix$b_x1 %*% t(x1_matrix[,i]) + 
        extSTfix$b_x2 %*% t(x2_matrix[,i])+ 
        extSTfix$phi 
      explogit <- exp(Logit_singletime)
      Theta_samples_Stan_singletime <-  explogit / (1 + explogit)
      Theta_STfix <- rbind(Theta_STfix, cbind(
        Mean=apply(Theta_samples_Stan_singletime,2,mean),
        t(apply(Theta_samples_Stan_singletime,2,quantile, 
                probs=c(0.025,0.975)))
      ))
    }
    
    
    print("AST")
    Logit_singletime <-  matrix(extSTcar$b_0,
                                nrow = length(extSTcar$b_0),
                                ncol = K, byrow = FALSE) +
      extSTcar$b_x1 %*% t(x1_matrix[,i]) + 
      extSTcar$b_x2 %*% t(x2_matrix[,i])+ 
      extSTcar$phi + 
      matrix(extSTcar$chi[,i], 
             ncol=K,
             nrow=length(extSTcar$chi[,i]))
    explogit <- exp(Logit_singletime)
    Theta_samples_Stan_singletime <- explogit / (1 + explogit)
    Theta_STcar <- rbind(Theta_STcar, cbind(
      Mean=apply(Theta_samples_Stan_singletime,2,mean),
      t(apply(Theta_samples_Stan_singletime,2,quantile, 
              probs=c(0.025,0.975)))
      # ,    Prob =prob[,i]
    ))
    
    print("ASTtrend")
    Logit_singletime <-  matrix(extSTtrend$b_0,
                                nrow = length(extSTtrend$b_0),
                                ncol = K, byrow = FALSE) +
      extSTtrend$b_x1 %*% t(x1_matrix[,i]) + 
      extSTtrend$b_x2 %*% t(x2_matrix[,i])+ 
      matrix( extSTtrend$b_3 *i,
              ncol=K,
              nrow=length(extSTtrend$b_3)) + 
      extSTtrend$phi + 
      matrix(extSTtrend$chi[,i], 
             ncol=K,
             nrow=length(extSTtrend$chi[,i]))
    explogit <- exp(Logit_singletime)
    Theta_samples_Stan_singletime <- explogit / (1 + explogit)
    Theta_STtrend <- rbind(Theta_STtrend, cbind(
      Mean=apply(Theta_samples_Stan_singletime,2,mean),
      t(apply(Theta_samples_Stan_singletime,2,quantile, 
              probs=c(0.025,0.975)))
      # ,    Prob =prob[,i]
    ))
    
  }
  if(!is.null(fitstan)){
  return(list(Theta_CARStan=Theta_CARStan,
              Theta_STcar=Theta_STcar, 
              Theta_STtrend=Theta_STtrend, 
              Theta_glm=Theta_glm, 
              Theta_STfix=Theta_STfix))
  }else{
    return(list(Theta_CARStan=Theta_CARStan,
                Theta_STcar=Theta_STcar, 
                Theta_STtrend=Theta_STtrend,
                Theta_STfix=Theta_STfix))
  }
}
