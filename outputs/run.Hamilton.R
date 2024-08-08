

path.results <- paste("results/res_estim_",ctry,".Rdat",sep="")
load(file=path.results)

# Construct estimated model:
model.est     <- Theta2Model(full_theta_est,model_ini_sol,bounds)
model.sol.est <- solve_model(model.est,indic_compute_beta = 1)

n_w <- dim(model.sol.est$Phi)[1]
n_x <- dim(model.sol.est$Phi_x)[1]

T     <- dim(DATA$CDS)[1]
all.X <- matrix(0,T,n_x)
all.X[,n_w+1] <- DATA$r
all.X[,n_w+2] <- DATA$d
all.X[,n_w+3] <- c(DATA$d[1],DATA$d[1:(T-1)])

load(file=paste("results/res_simul_",ctry,suffix,".Rdat",sep=""))

all.sim.ell1 <- NULL

for(count in 1:length(all.w)){
  
  full_theta_perturb <- all.theta[[count]]
  
  model.est <- Theta2Model(full_theta_perturb,model_ini_sol,bounds)
  model.sol.est <- solve_model(model.est,indic_compute_beta = 1)
  
  all.X[,1:n_w] <- all.w[[count]]
  
  ell1 <- compute.debt.limits(model.sol.est,all.X,Proba=.01,h=4)/model.sol.est$freq
  if(sum(ell1>=5)==0){
    all.sim.ell1 <- cbind(all.sim.ell1,ell1)
  }
}

