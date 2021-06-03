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
# library(shinystan)
# launch_shinystan(fitcarstan)

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
# is.positive.definite(distance)
# is.positive.definite(0.4 * exp(-0.1 * distance))

# time dimensions 
n_tmstp <- 10+4+4

#### Generate the covariates and response data
set.seed(125)
# x1 <- rnorm(K)+1
x1 <- runif(n = K*n_tmstp, min = 0, max = 1)
x1_matrix <- matrix(x1,ncol=n_tmstp, byrow = FALSE)
x1_test <- c(0,0.2, 0.5, 0.8,1)
set.seed(124)
# x2 <- rnorm(K)
x2 <- runif(n = K*n_tmstp, min = 0, max = 1)
x2_matrix <- matrix(x2,ncol=n_tmstp, byrow = FALSE)
x2_test <- c(0,0.2, 0.5, 0.8,1)
# theta <- rnorm(K, sd=0.05)
alpha <- 0.8# 1/median(apply(W,1,sum))
tau <- 0.1
D <- diag(apply(W,1,sum))
# tau *(D - alpha *W)[1:10, 1:10]
# round(solve(tau *(D - alpha *W))[1:10, 1:10],3)
set.seed(123)
# undebug(mvrnorm)
system.time({phi_s <- mvrnorm(mu=rep(0,K), Sigma=solve(tau *(D - alpha *W)))})
# user  system elapsed 
# 43.79    0.11   44.48 
# esempio <- 1946
# Grid[c(esempio, which(W[esempio, ]==1)),]
# phi_s[esempio]
# mean(phi_s[which(W[esempio,]==1)])

# add temporal dependence ------

ar1_par <- 0.7
b3 <- 0.1
set.seed(129)
phi_t <- arima.sim(model = list(ar=ar1_par), n = n_tmstp, n.start = 1000) +
  b3*(1:n_tmstp)


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
  cbind.data.frame(Prob =Prob_matrix[,1], Timestamp="01"),
  cbind.data.frame(Prob =Prob_matrix[,10], Timestamp="10"),
  cbind.data.frame(Prob =Prob_matrix[,14], Timestamp="14")
)
plots_ST <- NULL

plots_ST[[1]] <- ggplot(data.frame("Chi_t_Trend"= phi_t, Timestamp=1:n_tmstp),
                        aes(y=Chi_t_Trend, x=Timestamp)) +
  geom_line() + ylab(expression(chi[t] + beta[3] *t))
plots_ST[[2]] <- ggplot(Denistyevolution[,], aes(x=Prob, y=Timestamp)) +
  ggridges::geom_density_ridges(stat="binline",panel_scaling = FALSE, 
                                draw_baseline=FALSE) +
  xlab(expression(theta[t][k]))
pdf(file = "./010_IJF_code_2/ProbabilityTheta_STtrend_20210228.pdf")
do.call(gridExtra::grid.arrange, c(plots_ST,
                                   list(nrow=2,
                                        layout_matrix=matrix(c(1,2,2,2), ncol=1 ))
))
dev.off()



plots_ST[[1]] <- plots_ST[[1]]+
  theme(text = element_text(size=20))
plots_ST[[2]] <- plots_ST[[2]] +
  theme(text = element_text(size=20))
pdf(file = "./010_IJF_code_2/ProbabilityTheta_STtrend_bigger.pdf")
do.call(gridExtra::grid.arrange, c(plots_ST,
                                   list(nrow=2,
                                        layout_matrix=matrix(c(1,2,2,2), ncol=1 ))
))
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
fitglm <- glm(c(Y_matrix) ~ c(x1_matrix)+c(x2_matrix), #data = cbind.data.frame(Y, x1,x2), 
             family = binomial(link="logit"))
summary(fitglm)
#set params for MCMC -------
Iter = 2500
Chains = 4
Warmup = 1000
Verbose = TRUE
Save_dso=TRUE
Thin = 2

# test stan logit --------------
datalist <- list(N = K*ncol(Y_matrix), Y = c(Y_matrix),
                 x1 = c(x1_matrix), x2 = c(x2_matrix))

set.seed(127)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fitstan <- stan(file = "./010_IJF_code_2/Logit_training_2pred_v2.stan",
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
# require(coda)
set.seed(127)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fitcarstan <- stan(file = "./010_IJF_code_2/CARlogit_training_2pred_v2.stan",
                   data = datalist,iter = Iter, chains = Chains ,warmup = Warmup, verbose = Verbose,thin = Thin,
                   save_dso=Save_dso,
                   seed = 127 )


# test STAN SpaceTime stationary + intercept -------------
print(paste0("start ST at ", Sys.time()))
# shuold extract a 
datalist <- list(N = K, Time=n_tmstp,  Y = Y_matrix,
                 x1 = x1_matrix, x2 = x2_matrix, W=W, W_n = sum(W)/2)
# rstan::stan(file = "./simulated data.stan",data = datalist, iter = 1000, chains = 4,verbose = TRUE )
require(rstan)
require(dplyr)
# require(coda)
set.seed(127)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fitSTcar <- stan(file = paste0("./010_IJF_code_2/SpaceTimelogit_training_2pred_rhomarginal_v2.stan"),
                 data = datalist, iter = Iter, chains = Chains ,warmup = Warmup, verbose = Verbose,thin = Thin,
                 save_dso=Save_dso,
                 seed = 127 )

# test STAN SpaceTime stationary+ marginal + linear trend-------------
print(paste0("start ST trend at ", Sys.time()))

datalist <- list(N = K, Time=n_tmstp,  Y = Y_matrix,
                 x1 = x1_matrix, x2 = x2_matrix, W=W, W_n = sum(W)/2)
# rstan::stan(file = "./simulated data.stan",data = datalist, iter = 1000, chains = 4,verbose = TRUE )
require(rstan)
require(dplyr)
# require(coda)
set.seed(127)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fitSTtrend <- stan(file = "./010_IJF_code_2/SpaceTimelogit_training_2pred_rhomargTrend_v2.stan",
                   data = datalist, iter = Iter, chains = Chains ,warmup = Warmup, verbose = Verbose,thin = Thin,
                   save_dso=Save_dso,
                   seed = 127 )

save(logit, phi, prob, W, x1_matrix, x2_matrix, Y_matrix, 
     alpha, tau, ar1_par, b0_linear, b1_linear, b2_linear, b3, K, n_tmstp, 
     phi_s, phi_t, 
     fitglm, fitstan, fitcarstan, fitSTcar, fitSTtrend, 
     Y_fore_full, x1_fore_full, x2_fore_full,
     Iter ,    Chains,   Warmup,     Verbose,    Save_dso,     Thin,
     file="./010_IJF_code_2/Simulated_data_models_STrend.rdata")

# load("./010_IJF_code_2/Simulated_data_models_STrend.rdata")
# model diagnostics ----------
# library(rstan)
fitcarstanSum <- summary(fitcarstan,
                         pars = fitcarstan@sim$pars_oi,
                         probs = c(0.025, #0.25, 0.5, 
                                   # 0.75,
                                   0.975))$summary
Space <- cbind.data.frame(Parameter=c("b_0", "b_x1", "b_x2", "alpha_car"),
                          fitcarstanSum[ c("b_0", "b_x1", "b_x2", "alpha_car"), 
                                         c("mean", "2.5%", "97.5%")] )
colnames(Space)[-1] <- paste0(colnames(Space)[-1],"_SPACE")

fitSTcarSum <- summary(fitSTcar,
                       pars = fitSTcar@sim$pars_oi,
                       probs = c(0.025, #0.25, 0.5, 
                                 # 0.75,
                                 0.975))$summary
SpaceTime <- cbind.data.frame(Parameter=c("b_0", "b_x1", "b_x2", "ar1_par", "alpha_car"),
                              fitSTcarSum[c("b_0", "b_x1", "b_x2", "ar1_par", "alpha_car"),
                                          c("mean", "2.5%", "97.5%")])

colnames(SpaceTime)[-1] <- paste0(colnames(SpaceTime)[-1],"_SPACETime")

fitSTtrendSum <- summary(fitSTtrend,
                         pars = fitSTtrend@sim$pars_oi,
                         probs = c(0.025, #0.25, 0.5, 
                                   # 0.75,
                                   0.975))$summary
Trend <- cbind.data.frame(Parameter= c("b_0", "b_x1", "b_x2", "b_3", "ar1_par", "alpha_car"),
                          fitSTtrendSum[c("b_0", "b_x1", "b_x2", "b_3", "ar1_par", "alpha_car"),
                                        c("mean", "2.5%", "97.5%")])
colnames(Trend)[-1] <- paste0(colnames(Trend)[-1],"_STTrend")

True <-  cbind.data.frame(Parameter= c("b_0", "b_x1", "b_x2", "b_3", "ar1_par", "alpha_car"),
                          Value= c(b0_linear,b1_linear, b2_linear, b3 ,ar1_par, alpha))

GLM <- cbind.data.frame(Parameter= c("b_0", "b_x1", "b_x2"),
                        mean=fitglm$coefficients,
                        confint(fitglm))
colnames(GLM)[-1] <- paste0(colnames(GLM)[-1],"_GLM")


fitGLMSum <- summary(fitstan,
                     pars = fitstan@sim$pars_oi,
                     probs = c(0.025, #0.25, 0.5, 
                               # 0.75,
                               0.975))$summary
GLMstan <- cbind.data.frame(Parameter= c("b_0", "b_x1", "b_x2"),
                            fitGLMSum[c("b_0", "b_x1", "b_x2"),
                                      c("mean", "2.5%", "97.5%")])
colnames(GLMstan)[-1] <- paste0(colnames(GLMstan)[-1],"_GLMstan")


Scenario <- merge(True, Space, all = TRUE)
Scenario <- merge(Scenario, SpaceTime, all = TRUE)
Scenario <- merge(Scenario, Trend, all = TRUE)
Scenario <- merge(Scenario, GLMstan, all = TRUE)
Scenario <- merge(Scenario, GLM, all = TRUE)
tableSummaryTex(summarytab = Scenario,round = 2,rownames = NULL,
                label = "tab:summarySTrend",
                caption = "prova",
                filename =  "./010_IJF_code_2/simulationstudy_STtrend.tex")

# Performance in-sample -----
# Fitted values
Perc_REC <- c(0.03, 0.05, 0.1, 0.5)
Thresh_REC <- 0.5

predictionGLM <- predict(fitglm,type = "response")
Fitted <- ComputeFitted(fitcarstan = fitcarstan, fitSTcar = fitSTcar, fitSTtrend = fitSTtrend, 
                        fitstan = fitstan,
                        x1_matrix = x1_matrix, x2_matrix = x2_matrix)
save(  Fitted,
       file="./010_IJF_code_2/Simulated_data_fitted_STtrend.rdata")
# load("./010_IJF_code_2/Simulated_data_fitted_STtrend.rdata")
predNN <- W %*% Y_matrix[,-(8:10)]
PerfInsample <- MeasureFitted(predictionGLM = Fitted$Theta_glm, predNN = predNN,
                              Theta_CARStan = Fitted$Theta_CARStan, Theta_STcar = Fitted$Theta_STcar, Theta_STtrend = Fitted$Theta_STtrend,
                              Y_matrix =Y_matrix , 
                              Perc_REC =Perc_REC , Thresh_REC = Thresh_REC)

tableSummaryTex(summarytab = PerfInsample,round = 2 ,rownames = rownames(PerfInsample) ,
                label = "tab:PerfInsampleSTtrend",
                caption = "prova",
                filename =  "./010_IJF_code_2/simulationstudy_STtrend_PerfInsample.tex")

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
                           fitSTcar = fitSTcar, fitSTtrend = fitSTtrend, 
                           fitstan = fitstan,
                           W=W,Y_matrix=Y_matrix,
                           Perc_REC = Perc_REC, Thresh_REC = Thresh_REC, 
                           K =K )

tableSummaryTex(summarytab = Fore1$PerfOutTime,round = 2 ,rownames = rownames(Fore1$PerfOutTime) ,
                label = "tab:PerfOutTime_STtrend",
                caption = "prova",
                filename =  "./010_IJF_code_2/simulationstudy_STtrend_PerfOutTime.tex")


#Performance out of time 2---------
Y_matrix <- cbind(Y_matrix,Y_fore_full[,1:4])
x1_matrix <- cbind(x1_matrix,x1_fore_full[,1:4])
x2_matrix <- cbind(x2_matrix,x2_fore_full[,1:4])
Y_fore <- Y_fore_full[,4*2]
x1_fore <- x1_fore_full[,4*2]
x2_fore <- x2_fore_full[,4*2]    

# retrain glm
set.seed(127)
fitglm2 <- glm(c(Y_matrix) ~ c(x1_matrix)+c(x2_matrix), #data = cbind.data.frame(Y, x1,x2), 
               family = binomial(link="logit"))

# retrain stan logit --------------
datalist <- list(N = K*ncol(Y_matrix), Y = c(Y_matrix),
                 x1 = c(x1_matrix), x2 = c(x2_matrix))

set.seed(127)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fitstan2 <- stan(file = "./010_IJF_code_2/Logit_training_2pred_v2.stan",
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
# require(coda)
set.seed(127)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fitcarstan2 <- stan(file = "./010_IJF_code_2/CARlogit_training_2pred_v2.stan",
                    data = datalist, iter = Iter, chains = Chains ,warmup = Warmup, verbose = Verbose,thin = Thin,
                    save_dso=Save_dso,
                    seed = 127 )


# retrain ST -------------------
print(paste0("start ST at ", Sys.time()))
# shuold extract a 
datalist <- list(N = K, Time=ncol(Y_matrix),  Y = Y_matrix,
                 x1 = x1_matrix, x2 = x2_matrix, W=W, W_n = sum(W)/2)
# rstan::stan(file = "./simulated data.stan",data = datalist, iter = 1000, chains = 4,verbose = TRUE )
require(rstan)
require(dplyr)
# require(coda)
set.seed(127)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fitSTcar2 <- stan(file = paste0("./010_IJF_code_2/SpaceTimelogit_training_2pred_rhomarginal_v2.stan"),
                  data = datalist, iter = Iter, chains = Chains ,warmup = Warmup, verbose = Verbose,thin = Thin,
                  save_dso=Save_dso,
                  seed = 127 )

# retrain STtrend ------------------
print(paste0("start ST trend at ", Sys.time()))

datalist <- list(N = K, Time=ncol(Y_matrix),  Y = Y_matrix,
                 x1 = x1_matrix, x2 = x2_matrix, W=W, W_n = sum(W)/2)
# rstan::stan(file = "./simulated data.stan",data = datalist, iter = 1000, chains = 4,verbose = TRUE )
require(rstan)
require(dplyr)
# require(coda)
set.seed(127)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fitSTtrend2 <- stan(file = "./010_IJF_code_2/SpaceTimelogit_training_2pred_rhomargTrend_v2.stan",
                    data = datalist, iter = Iter, chains = Chains ,warmup = Warmup, verbose = Verbose,thin = Thin,
                    save_dso=Save_dso,
                    seed = 127 )

save(fitglm2, fitstan2, fitcarstan2, fitSTcar2, fitSTtrend2, 
     file="./010_IJF_code_2/Simulated_data_models_retrain_STtrend.rdata")
# load("./010_IJF_code_2/Simulated_data_models_retrain_Space.rdata")

Fore2 <-ForecastandMeasure(future_dist = future_dist,PredictedTime= (ncol(x1_matrix)+future_dist),
                           x1_fore = x1_fore, x2_fore = x2_fore, Y_fore = Y_fore,
                           fitglm = fitglm2, fitcarstan = fitcarstan2,
                           fitSTcar = fitSTcar2, fitSTtrend = fitSTtrend2,  
                           fitstan = fitstan2,
                           W=W,Y_matrix=Y_matrix,
                           Perc_REC = Perc_REC, Thresh_REC = Thresh_REC, K=K)


PerfOutTime_multiple <-(Fore1$PerfOutTime+ Fore2$PerfOutTime)/2
rownames(PerfOutTime_multiple) <- c("AUC", "$REC\\_{0.03} (0-4.69)$",
                                    "$REC\\_{0.05} (0-7.82)$",
                                    "$REC\\_{0.1} (0-15.62)$" ,
                                    "$REC\\_{0.5} (0-78.11)$" , 
                                    "Recall", "Precision" )

tableSummaryTex(summarytab = PerfOutTime_multiple,round = 2 ,rownames = rownames(PerfOutTime_multiple) ,
                label = "tab:PerfOutTime_multiple_STtrend",
                caption = "prova",
                filename =  "./010_IJF_code_2/simulationstudy_STtrend_PerfOutTime_multiple.tex")
