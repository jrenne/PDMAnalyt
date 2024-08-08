# ==============================================================================
# Script that runs parallel simulations for countries in list.of.ctries
# ==============================================================================

cl <- makeCluster(number.of.cores)
registerDoParallel(cl)

save.image("results/toto.Rdata")
clusterEvalQ(cl,load("results/toto.Rdata"))

clusterEvalQ(cl,library(optimx))
clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(expm))
clusterEvalQ(cl,library(seasonal))
clusterEvalQ(cl,library(stringr))
clusterEvalQ(cl,library(zoo))
clusterEvalQ(cl,library(mFilter))
clusterEvalQ(cl,library(splines))
clusterEvalQ(cl,library(Rcpp))
clusterEvalQ(cl,library(RcppEigen))
clusterEvalQ(cl,library(pracma))

foreach(i = 1:length(list.of.ctries),
        .combine = "rbind",
        .noexport = c("compute_cds_cpp",
                      "compute_condit_Exp_cpp",
                      "solve_model",
                      "EKF_filter_cpp"))  %dopar% {
                        
                        sourceCpp("procedures/pricing_cpp.cpp") # load C++ Kalman Filter
                        
                        ctry <- list.of.ctries[i]
                        
                        # Initialize random number generator:
                        set.seed(123)

                        path.results <- paste("results/res_estim_",ctry,suffix,".Rdat",sep="")
                        load(file=path.results)
                        # Construct estimated model:
                        model.est     <- Theta2Model(full_theta_est,model_ini_sol,bounds)
                        model.sol.est <- solve_model(model.est,indic_compute_beta = 1)
                        resEKF_est    <- EKF_smoother(model.sol.est,list.stdv,DATA)
                        logl0         <- resEKF_est$loglik

                        all.w     <- list()
                        all.theta <- list()
                        
                        # Load covariance matrix of parameter estimates:
                        load(file=paste("results/res_MatCov_",ctry,suffix,".Rdat",sep=""))
                        CovMat <- res.stdv$Mat.var.cov
                        
                        indic.zero <- (apply(CovMat,1,sum)==0)
                        if(length(indic.zero)>0){
                          CovMat[indic.zero,indic.zero] <- .001
                        }
                        Sigma  <- t(chol(CovMat))
                        
                        count <- 0
                        
                        while(count < N.replic){
                          full_theta_perturb <- full_theta_est
                          
                          list.perturb.theta <- res.stdv$List.param.estim.2keep
                          
                          perturb.theta <- Sigma %*% rnorm(length(list.perturb.theta))
                          full_theta_perturb[list.perturb.theta,] <-
                            full_theta_perturb[list.perturb.theta,] + perturb.theta
                          
                          model.est     <- Theta2Model(full_theta_perturb,model.sol.est,bounds)
                          model.sol.est <- solve_model(model.est,indic_compute_beta = 1)
                          
                          resEKF_est <- EKF_smoother(model.sol.est,list.stdv,DATA)
                          
                          logl1 <- resEKF_est$loglik
                          print(logl1)
                          
                          if(!is.na(logl1)){
                            if(abs(logl0-logl1)<max.diff.logl){# otherwise implausible parameters
                              count <- count + 1
                              all.w[[count]] <- t(resEKF_est$W.smoothed)
                              all.theta[[count]] <- full_theta_perturb
                            }
                          }
                        }
                        
                        # Save results:
                        save(all.w,all.theta,
                             file=paste("results/res_simul_",ctry,suffix,".Rdat",sep=""))
                        
                        success <- 1
                      }

stopCluster(cl)


