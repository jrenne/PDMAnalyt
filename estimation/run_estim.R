
print("---- New model estimation (macro block) ----")

if(indic_use_last_res){
  load(file=paste("results/",file_with_input_param,sep=""))
}else{
  source("estimation/set_ini_model.R")
}

Model_ini <- Model

param <- model2param(Model_ini)
# Randomize initial values:
param <- param + random_factor * rnorm(param)

for(ii in 1:nb_loops){
  for(algo in c("Nelder-Mead","nlminb")){
    res.estim <- optimx(param,
                        #compute_distance,
                        #compute_loglik,
                        compute_total_distance,
                        targets = targets,
                        Model_ini = Model_ini,
                        method = algo,
                        control=
                          list(trace=0,
                               maxit=ifelse(algo=="Nelder-Mead",
                                            maxit.NlMd,maxit.nlminb),
                               kkt = FALSE))
    param     <- c(as.matrix(res.estim)[1:length(param)])
    print(paste("Log-lik value at end of iteration: ",-round(res.estim$value,2),sep=""))
  }
}

compute_total_distance(param,targets,Model_ini)

Model <- make_model(param,Model_ini)

source("outputs/make_chart_fit.R")

