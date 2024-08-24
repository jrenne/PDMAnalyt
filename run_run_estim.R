
# Optimization set up:
maxit.nlminb <- 50
maxit.NlMd   <- 2000
nb_loops     <- 2

nb_attemps <- 20

best <- 100000

for(iii in 1:nb_attemps){
  source("estimation/run_estim.R")
  if(res.estim["value"]<best){
    best_param <- param
    best <- res.estim["value"]
  }
}


