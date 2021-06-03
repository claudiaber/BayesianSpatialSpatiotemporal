print(Sys.time())
rm(list=ls())
library(MASS)
library(stringr)
library(rstan)
library(dplyr)

source("./Functions/tableSummaryTex.R")
source("./Functions/theta_hat.R")
source("./Functions/ComputeFitted.R")
source("./Functions/MeasureFitted.R")
source("./Functions/ForecastandMeasure.R")

#### Set up a square lattice region ------
griddim <- 50
x.easting <- 1:griddim
x.northing <- 1:griddim
Grid <- expand.grid(x.easting, x.northing)
K <- nrow(Grid)

#### set up distance and neighbourhood (W, based on sharing a common border) matrices ----
distance <- as.matrix(dist(Grid))
W <-array(0, c(K,K))
W[distance==1| distance==(griddim-1)] <-1 	
isSymmetric(distance)

# time dimensions 
n_tmstp <- 10+4+4

#### Generate the covariates and response data
set.seed(125)
x1 <- runif(n = K*n_tmstp, min = 0, max = 1)
x1_matrix <- matrix(x1,ncol=n_tmstp, byrow = FALSE)
x1_test <- c(0,0.2, 0.5, 0.8,1)
set.seed(124)
x2 <- runif(n = K*n_tmstp, min = 0, max = 1)
x2_matrix <- matrix(x2,ncol=n_tmstp, byrow = FALSE)
x2_test <- c(0,0.2, 0.5, 0.8,1)
alpha <- 0.8
tau <- 0.1
D <- diag(apply(W,1,sum))
set.seed(123)
system.time({phi_s <- mvrnorm(mu=rep(0,K), Sigma=solve(tau *(D - alpha *W)))})
# user  system elapsed 
# 43.79    0.11   44.48 

# add temporal dependence ------

ar1_par <- 0
set.seed(128)
phi_t <- rep(0, n_tmstp)

# sum spatial and temporal effects
phi_linear <- expand.grid(phi_s=phi_s, phi_t=phi_t)
phi <- matrix(phi_linear$phi_s + phi_linear$phi_t, ncol = length(phi_t), byrow = FALSE)
rm(phi_linear)

# compute prob 
b0_linear <- 1
b1_linear <-4
b2_linear <- -6
logit <- b0_linear + b1_linear *x1_matrix + b2_linear *x2_matrix + phi
prob <- exp(logit) / (1 + exp(logit))
trials <- rep(1,K)
set.seed(126)
Y <- rbinom(n=K*n_tmstp, size=trials, prob=prob)
Y_matrix <- matrix(Y, ncol=n_tmstp, byrow = FALSE)


Prob_matrix <- matrix(prob, ncol=n_tmstp, byrow = FALSE)
Denistyevolution <- rbind(
  cbind.data.frame(Prob =Prob_matrix[,1], Timestamp="01")
)

plots_ST <- NULL

plots_ST[[1]] <- ggplot(Denistyevolution[,], aes(x=Prob, y=Timestamp)) +
  ggridges::geom_density_ridges(stat="binline",panel_scaling = FALSE, 
                                draw_baseline=FALSE)+
  xlab(expression(theta[t][k]))

plots_ST[[1]] <- plots_ST[[1]]+
  theme(text = element_text(size=20))
pdf(file = "./Figures/ProbabilityTheta_SPACE_bigger.pdf", height = 2)
do.call(gridExtra::grid.arrange, c(plots_ST))
dev.off()

# forecast data
Y_fore_full <- Y_matrix[,11:n_tmstp]
x1_fore_full <- x1_matrix[,11:n_tmstp]
x2_fore_full <- x2_matrix[,11:n_tmstp]

# training data 
n_tmstp <- 10

Y_matrix <- Y_matrix[,1:n_tmstp]
x1_matrix <- x1_matrix[,1:n_tmstp]
x2_matrix <- x2_matrix[,1:n_tmstp]

#test base glm logit --------------
set.seed(127)
fitglm <- glm(c(Y_matrix) ~ c(x1_matrix)+c(x2_matrix),
             family = binomial(link="logit"))
summary(fitglm)

#set params for MCMC -------
Iter = 100 #2500
Chains = 4
Warmup = 50 #1000
Verbose = TRUE
Save_dso=TRUE
Thin = 2

# test stan logit --------------
datalist <- list(N = K*ncol(Y_matrix), Y = c(Y_matrix),
                 x1 = c(x1_matrix), x2 = c(x2_matrix))

set.seed(127)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fitstan <- stan(file = "./STAN/Logit_training_2pred_v2.stan",
                data = datalist, iter = Iter, chains = Chains ,warmup = Warmup, verbose = Verbose,thin = Thin,
                save_dso=Save_dso,
                seed = 127)
print(fitstan)


# test STAN CAR -------------
print(paste0("start CAR at ", Sys.time()))
datalist <- list(N = K, Y = Y_matrix[,10],
                 x1 = x1_matrix[,10], x2 = x2_matrix[,10], 
                 W=as.matrix(W), W_n = sum(W)/2)
require(rstan)
require(dplyr)
set.seed(127)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fitcarstan <- stan(file = "./STAN/CARlogit_training_2pred_v2.stan",
                   data = datalist,iter = Iter, chains = Chains ,warmup = Warmup, verbose = Verbose,thin = Thin,
                   save_dso=Save_dso,
                   seed = 127 )


# test STAN SpaceTime stationary + intercept -------------
print(paste0("start ST at ", Sys.time()))
datalist <- list(N = K, Time=n_tmstp,  Y = Y_matrix,
                 x1 = x1_matrix, x2 = x2_matrix, W=W, W_n = sum(W)/2)
require(rstan)
require(dplyr)
set.seed(127)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fitSTcar <- stan(file = paste0("./STAN/SpaceTimelogit_training_2pred_rhomarginal_v2.stan"),
                   data = datalist, iter = Iter, chains = Chains ,warmup = Warmup, verbose = Verbose,thin = Thin,
                 save_dso=Save_dso,
                   seed = 127 )

# test STAN SpaceTime stationary+ marginal + linear trend-------------
print(paste0("start ST trend at ", Sys.time()))
datalist <- list(N = K, Time=n_tmstp,  Y = Y_matrix,
                 x1 = x1_matrix, x2 = x2_matrix, W=W, W_n = sum(W)/2)
require(rstan)
require(dplyr)
set.seed(127)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fitSTtrend <- stan(file = "./STAN/SpaceTimelogit_training_2pred_rhomargTrend_v2.stan",
                   data = datalist, iter = Iter, chains = Chains ,warmup = Warmup, verbose = Verbose,thin = Thin,
                   save_dso=Save_dso,
                   seed = 127 )

# test STAN SpaceTime with fixed phis -------------
print(paste0("start ST fix phi at ", Sys.time()))
datalist <- list(N = K, Time=n_tmstp,  Y = Y_matrix,
                 x1 = x1_matrix, x2 = x2_matrix, W=W, W_n = sum(W)/2)
require(rstan)
require(dplyr)
# require(coda)
set.seed(127)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fitSTcarfix <- stan(file = paste0("./STAN/SpaceTimelogit_training_2pred_fixphi.stan"),
                    data = datalist, iter = Iter, chains = Chains ,warmup = Warmup, verbose = Verbose,thin = Thin,
                    save_dso=Save_dso,
                    seed = 127 )

save(logit, phi, prob, W, x1_matrix, x2_matrix, Y_matrix, 
     alpha, tau, ar1_par, b0_linear, b1_linear, b2_linear, K, n_tmstp, 
     phi_s, phi_t, 
     fitglm, fitstan, fitcarstan, fitSTcar, fitSTtrend, 
     fitSTcarfix,
     Y_fore_full, x1_fore_full, x2_fore_full,
     Iter ,    Chains,   Warmup,     Verbose,    Save_dso,     Thin,
     file="./Posterior/Simulated_data_models_Space.rdata")
# load("./Posterior/Simulated_data_models_Space.rdata")

# model diagnostics ----------
library(rstan)
fitcarstanSum <- summary(fitcarstan,
                         pars = fitcarstan@sim$pars_oi,
                         probs = c(0.025,
                                   0.975))$summary
Space <- cbind.data.frame(Parameter=c("b_0", "b_x1", "b_x2", "alpha_car"),
                 fitcarstanSum[ c("b_0", "b_x1", "b_x2", "alpha_car"), 
                                c("mean", "2.5%", "97.5%")] )
colnames(Space)[-1] <- paste0(colnames(Space)[-1],"_SPACE")

fitcarfix <- summary(fitSTcarfix,
                         pars = fitSTcarfix@sim$pars_oi,
                         probs = c(0.025, 
                                   0.975))$summary
Timefix <- cbind.data.frame(Parameter=c("b_0", "b_x1", "b_x2", "alpha_car"),
                          fitcarfix[ c("b_0", "b_x1", "b_x2", "alpha_car"), 
                                         c("mean", "2.5%", "97.5%")] )
colnames(Timefix)[-1] <- paste0(colnames(Space)[-1],"_TimeFix")


fitSTcarSum <- summary(fitSTcar,
                         pars = fitSTcar@sim$pars_oi,
                         probs = c(0.025, 
                                   0.975))$summary
SpaceTime <- cbind.data.frame(Parameter=c("b_0", "b_x1", "b_x2", "ar1_par", "alpha_car"),
                              fitSTcarSum[c("b_0", "b_x1", "b_x2", "ar1_par", "alpha_car"),
                                          c("mean", "2.5%", "97.5%")])

colnames(SpaceTime)[-1] <- paste0(colnames(SpaceTime)[-1],"_SPACETime")

fitSTtrendSum <- summary(fitSTtrend,
                         pars = fitSTtrend@sim$pars_oi,
                         probs = c(0.025, 
                                   0.975))$summary
Trend <- cbind.data.frame(Parameter= c("b_0", "b_x1", "b_x2", "b_3", "ar1_par", "alpha_car"),
                          fitSTtrendSum[c("b_0", "b_x1", "b_x2", "b_3", "ar1_par", "alpha_car"),
                                        c("mean", "2.5%", "97.5%")])
colnames(Trend)[-1] <- paste0(colnames(Trend)[-1],"_STTrend")

True <-  cbind.data.frame(Parameter= c("b_0", "b_x1", "b_x2", "b_3", "ar1_par", "alpha_car"),
                          Value= c(b0_linear,b1_linear, b2_linear, 0,ar1_par, alpha))

GLM <- cbind.data.frame(Parameter= c("b_0", "b_x1", "b_x2"),
                 mean=fitglm$coefficients,
                 confint(fitglm))
colnames(GLM)[-1] <- paste0(colnames(GLM)[-1],"_GLM")

fitGLMSum <- summary(fitstan,
                         pars = fitstan@sim$pars_oi,
                         probs = c(0.025, 
                                   0.975))$summary
GLMstan <- cbind.data.frame(Parameter= c("b_0", "b_x1", "b_x2"),
                            fitGLMSum[c("b_0", "b_x1", "b_x2"),
                                        c("mean", "2.5%", "97.5%")])
colnames(GLMstan)[-1] <- paste0(colnames(GLMstan)[-1],"_GLMstan")


Scenario <- merge(True, Space, all = TRUE)
Scenario <- merge(Scenario, Timefix, all = TRUE)
Scenario <- merge(Scenario, SpaceTime, all = TRUE)
Scenario <- merge(Scenario, Trend, all = TRUE)
Scenario <- merge(Scenario, GLMstan, all = TRUE)
Scenario <- merge(Scenario, GLM, all = TRUE)
tableSummaryTex(summarytab = Scenario,round = 2,rownames = NULL,
                label = "tab:summarySpace",
                caption = "prova",
                filename =  "./Tables/simulationstudy_Space.tex")

# Performance in-sample -----
# Fitted values
Perc_REC <- c(0.03, 0.05, 0.1, 0.5)
Thresh_REC <- 0.5

predictionGLM <- predict(fitglm,type = "response")
Fitted <- ComputeFitted(fitcarstan = fitcarstan,fitSTcarfix=fitSTcarfix, fitSTcar = fitSTcar, fitSTtrend = fitSTtrend, 
                        fitstan = fitstan,
              x1_matrix = x1_matrix, x2_matrix = x2_matrix)


save(  Fitted,
       file="./Posterior/Simulated_data_fitted_Space.rdata")
# load("./010_IJF_code_2/Newmodel/Simulated_data_fitted_Space.rdata")

predNN <- W %*% Y_matrix[,-(8:10)]
PerfInsample <- MeasureFitted(predictionGLM = Fitted$Theta_glm, predNN = predNN, 
              Theta_CARStan = Fitted$Theta_CARStan, Theta_STcar = Fitted$Theta_STcar, Theta_STtrend = Fitted$Theta_STtrend,
              Theta_STfix= Fitted$Theta_STfix,
              Y_matrix =Y_matrix , 
              Perc_REC =Perc_REC , Thresh_REC = Thresh_REC)

tableSummaryTex(summarytab = PerfInsample,round = 2 ,rownames = rownames(PerfInsample) ,
                label = "tab:PerfInsampleSpace",
                caption = "prova",
                filename =  "./Posterior/simulationstudy_Space_PerfInsample.tex")
                
# Performance out of time 1-------
Perc_REC <- c(0.03, 0.05, 0.1, 0.5)
Thresh_REC <- 0.5
#update chi
future_dist <- 4# 2018 09 --> 2019 01

x1_fore <- x1_fore_full[,4*1]
x2_fore <- x2_fore_full[, 4*1]
Y_fore <-Y_fore_full[,4*1]

Fore1 <-ForecastandMeasure(future_dist = future_dist, PredictedTime= (ncol(x1_matrix)+future_dist),
                   x1_fore = x1_fore, x2_fore = x2_fore, Y_fore = Y_fore,
                   fitglm = fitglm, fitcarstan = fitcarstan,
                   fitSTcarfix=fitSTcarfix,
                   fitSTcar = fitSTcar, fitSTtrend = fitSTtrend,
                   fitstan = fitstan,
                   W=W,Y_matrix=Y_matrix,
                   Perc_REC = Perc_REC, Thresh_REC = Thresh_REC, 
                   K =K )

tableSummaryTex(summarytab = Fore1$PerfOutTime,round = 2 ,rownames = rownames(Fore1$PerfOutTime) ,
                label = "tab:PerfOutTime_Space",
                caption = "prova",
                filename =  "./Tables/simulationstudy_Space_PerfOutTime.tex")


#Performance out of time 2---------
Y_matrix <- cbind(Y_matrix,Y_fore_full[,1:4])
x1_matrix <- cbind(x1_matrix,x1_fore_full[,1:4])
x2_matrix <- cbind(x2_matrix,x2_fore_full[,1:4])
Y_fore <- Y_fore_full[,4*2]
x1_fore <- x1_fore_full[,4*2]
x2_fore <- x2_fore_full[,4*2]    

# retrain glm
set.seed(127)
fitglm2 <- glm(c(Y_matrix) ~ c(x1_matrix)+c(x2_matrix), 
              family = binomial(link="logit"))



# retrain stan logit --------------
datalist <- list(N = K*ncol(Y_matrix), Y = c(Y_matrix),
                 x1 = c(x1_matrix), x2 = c(x2_matrix))

set.seed(127)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fitstan2 <- stan(file = "./STAN/Logit_training_2pred_v2.stan",
                data = datalist, iter = Iter, chains = Chains ,warmup = Warmup, verbose = Verbose,thin = Thin,
                save_dso=Save_dso,
                seed = 127)
print(fitstan2)


# retrain CAR -------------
print(paste0("Start CAR at ", Sys.time()))
datalist <- list(N = K, Y = Y_matrix[,ncol(Y_matrix)],
                 x1 = x1_matrix[,ncol(Y_matrix)], x2 = x2_matrix[,ncol(Y_matrix)], 
                 W=as.matrix(W), W_n = sum(W)/2)
require(rstan)
require(dplyr)
set.seed(127)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fitcarstan2 <- stan(file = "./STAN/CARlogit_training_2pred_v2.stan",
                   data = datalist, iter = Iter, chains = Chains ,warmup = Warmup, verbose = Verbose,thin = Thin,
                   save_dso=Save_dso,
                   seed = 127 )

# retrain ST fix -----------------
print(paste0("start ST fix at ", Sys.time()))
datalist <- list(N = K, Time=ncol(Y_matrix),  Y = Y_matrix,
                 x1 = x1_matrix, x2 = x2_matrix, W=W, W_n = sum(W)/2)
require(rstan)
require(dplyr)
set.seed(127)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fitSTcarfix2 <- stan(file = paste0("./STAN/SpaceTimelogit_training_2pred_fixphi.stan"),
                  data = datalist, iter = Iter, chains = Chains ,warmup = Warmup, verbose = Verbose,thin = Thin,
                  save_dso=Save_dso,
                  seed = 127 )



# retrain ST -------------------
print(paste0("start ST at ", Sys.time()))
# shuold extract a 
datalist <- list(N = K, Time=ncol(Y_matrix),  Y = Y_matrix,
                 x1 = x1_matrix, x2 = x2_matrix, W=W, W_n = sum(W)/2)
require(rstan)
require(dplyr)
set.seed(127)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fitSTcar2 <- stan(file = paste0("./STAN/SpaceTimelogit_training_2pred_rhomarginal_v2.stan"),
                 data = datalist, iter = Iter, chains = Chains ,warmup = Warmup, verbose = Verbose,thin = Thin,
                 save_dso=Save_dso,
                 seed = 127 )

# retrain STtrend ------------------
print(paste0("start ST trend at ", Sys.time()))

datalist <- list(N = K, Time=ncol(Y_matrix),  Y = Y_matrix,
                 x1 = x1_matrix, x2 = x2_matrix, W=W, W_n = sum(W)/2)
require(rstan)
require(dplyr)
set.seed(127)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fitSTtrend2 <- stan(file = "./STAN/SpaceTimelogit_training_2pred_rhomargTrend_v2.stan",
                   data = datalist, iter = Iter, chains = Chains ,warmup = Warmup, verbose = Verbose,thin = Thin,
                   save_dso=Save_dso,
                   seed = 127 )

save(fitglm2, fitstan2, fitcarstan2, fitSTcar2,fitSTcarfix2, fitSTtrend2, 
     file="./Posterior/Simulated_data_models_retrain_Space.rdata")
# load("./Posterior/Simulated_data_models_retrain_Space.rdata")

Fore2 <-ForecastandMeasure(future_dist = future_dist,PredictedTime= (ncol(x1_matrix)+future_dist),
                           x1_fore = x1_fore, x2_fore = x2_fore, Y_fore = Y_fore,
                           fitglm = fitglm2, fitcarstan = fitcarstan2,
                           fitSTcar = fitSTcar2, fitSTtrend = fitSTtrend2, 
                          fitstan = fitstan2,fitSTcarfix=fitSTcarfix2,
                           W=W,Y_matrix=Y_matrix,
                           Perc_REC = Perc_REC, Thresh_REC = Thresh_REC, K=K)


PerfOutTime_multiple <-(Fore1$PerfOutTime+ Fore2$PerfOutTime)/2

#to be computed
rownames(PerfOutTime_multiple) <- c("AUC", "$REC\\_{0.03} (0-5.84)$",
                                    "$REC\\_{0.05} (0-9.74)$",
                                    "$REC\\_{0.1} (0-19.48)$" ,
                                    "$REC\\_{0.5} (0-97.39)$" ,
                                    "Recall", "Precision" )


tableSummaryTex(summarytab = PerfOutTime_multiple,round = 2 ,rownames = rownames(PerfOutTime_multiple) ,
                label = "tab:PerfOutTime_multiple_Space",
                caption = "prova",
                filename =  "./Tables/simulationstudy_Space_PerfOutTime_multiple.tex")
