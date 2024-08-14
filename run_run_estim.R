
nb_attemps <- 20

best <- 100000

for(iii in 1:nb_attemps){
  source("estimation/run_estim.R")
  if(res.estim["value"]<best){
    best_param <- param
  }
}


