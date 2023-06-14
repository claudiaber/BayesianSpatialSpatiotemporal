# BayesianSpatialSpatiotemporal
In this repo you may find the code to reproduce results on simulated data presented in this paper 
Claudia Berloco, Raffaele Argiento, Silvia Montagna,
Forecasting short-term defaults of firms in a commercial network via Bayesian spatial and spatio-temporal methods,
International Journal of Forecasting,
2022,
ISSN 0169-2070,
https://doi.org/10.1016/j.ijforecast.2022.05.003.
(https://www.sciencedirect.com/science/article/pii/S0169207022000632)

The three scripts in the root excute the similar pipelines and test three different models on 3 different datasets:
1) data with spatial autocorrelation (".._Space.R"), 
2) data with spatio-temporal autocorrelation ("... _ST.R")
3) data with spatio-temporal autocorrelation and a trend("..._STtrend.R)

For each dataset three models are fitted and tested out-of-time on two future timestamps. Three models are defined and excuted in STAN:
1) a spatial model
2) a spatio temporal model 
3) a spatio temporal model with a trend component 

In each pipeline, the code does the following:
a) import functions from folder /Functions;
b) simulate data;
c) fit a Bayesian model in STAN (STAN code is in the /STAN folder);
d) evaluate the model insample, out of time on the first timestamp;
e) compute metrics;
f) retrain and evaluate on the second timestamp.
