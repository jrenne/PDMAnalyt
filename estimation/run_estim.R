
print("---- New model estimation (macro block) ----")

# Resize dataset:
indic_fst <- which(format(DATA$date,"%Y")==start_year)
DATA <- DATA[indic_fst:dim(DATA)[1],]

# Targetted moments:
targets <- list(
  target_10_nom = mean(DATA$SVENY10,na.rm=TRUE)/100,
  target_01_nom = mean(DATA$SVENY01,na.rm=TRUE)/100,
  target_10_rea = mean(DATA$TIPSY10,na.rm=TRUE)/100,
  target_02_rea = mean(DATA$TIPSY02,na.rm=TRUE)/100,
  target_slop_nom = mean(DATA$SVENY10,na.rm=TRUE)/100 - 
    mean(DATA$SVENY01,na.rm=TRUE)/100,
  target_slop_rea = mean(DATA$TIPSY10,na.rm=TRUE)/100 -
    mean(DATA$TIPSY02,na.rm=TRUE)/100,
  target_avg_Pi = mean(DATA$pi,na.rm=TRUE),
  target_avg_Dy = mean(DATA$dy,na.rm=TRUE),
  target_std_10_nom = sd(DATA$SVENY10,na.rm=TRUE)/100,
  target_std_10_rea = sd(DATA$TIPSY10,na.rm=TRUE)/100,
  IRP10 = 0 # 10-year inflation risk premium
  )

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

