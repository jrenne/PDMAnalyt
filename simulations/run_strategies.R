# ==============================================================================
# Analysis of public debt management strategies
# ==============================================================================


outputs <- c("mean_d","stdv_d","DaR95","mean_rr","stdv_rr",
             "stdv_Delta_d","avg_PD[maxH]","avg_spreads[maxH]")
# Convert the names of output into "expressions" for charts labels & titles:
outputs4charts <- c(expression(paste(E(d),sep="")),
                    expression(paste(sqrt(V(d)),sep="")),
                    expression(paste(q[95](d),sep="")),
                    expression(paste(E(r),sep="")),
                    expression(paste(sqrt(V(r)),sep="")),
                    expression(paste(sqrt(Delta(d)),sep="")),
                    expression(paste(E(PD[10]),sep="")),
                    expression(paste(E(spd[10]),sep="")))
# For Latex table:
latex.column.names <- c("$\\mathbb{E}(d)$","$\\sqrt{\\mathbb{V}(d)}$",
                        "$q_{95}(d)$",
                        "$\\mathbb{E}(r)$","$\\sqrt{\\mathbb{V}(r)}$",
                        "$\\sqrt{\\mathbb{V}(\\Delta d)}$",
                        "$\\mathbb{E}(PD)$","$\\mathbb{E}(spd)$")


grids <- make_grid(nb_grid = nb_grid,
                   min_d  = min_d,
                   max_d  = max_d,
                   min_rr = min_rr,
                   max_rr = max_rr,
                   sigma_eps = Model$sigma_eps,
                   all_quantiles_eps = c(-2,-1,1,2))

M <- NULL
parameters <- NULL

for(chi in values_of_chi){
  for(kappa_pi in values_of_kappa_pi){
    for(kappa_y in values_of_kappa_y){
      
      print(paste("chi = ",chi," , kappa_pi = ",kappa_pi,
                  " , kappa_y = ",kappa_y,sep=""))
      Model_i <- Model
      Model_i$kappa_pi <- kappa_pi
      Model_i$kappa_y  <- kappa_y
      Model_i$chi      <- chi
      
      # ============================================
      # To model liquidity risks:
      Model_i$delta <- Model$delta - .0 * kappa_pi
      # ============================================
      
      Model_solved_i <- solve_ToyModel(Model_i,grids,
                                       nb_iter = nb_iter,
                                       nb_iter_sdf = nb_iter_sdf)
      
      # -----------------------------------
      p <- compute_uncond_distri(Model_solved_i$indicators_x,
                                 Model_solved_i$Probas,1000)
      distri_d  <- compute_distri_x(grids$all_d,Model_solved_i$d,p)
      plot(grids$all_d,distri_d,type="l")
      distri_rr  <- compute_distri_x(grids$all_rr,Model_solved_i$rr,p)
      plot(grids$all_rr,distri_rr,type="l")
      # -----------------------------------
      
      strat_i <- run_strategy(Model_solved_i,maxH=10,
                              nb_iter_sdf = nb_iter_sdf)
      
      thisLine <- NULL
      for(output in outputs){
        eval(parse(text = gsub(" ","",paste("aux <- strat_i$",output,sep=""))))
        thisLine <- c(thisLine,aux)
      }
      
      M <- rbind(M,thisLine)
      parameters <- rbind(parameters,
                          c(chi,kappa_pi,kappa_y))
    }
  }
}

colnames(parameters) <- c("chi","kappa_pi","kappa_y")
colnames(M)          <- outputs

plot(M[,"mean_d"],M[,"stdv_d"])
plot(M[,"mean_rr"],M[,"avg_PD[maxH]"])
plot(M[,"mean_rr"],M[,"DaR95"])

plot(M[(parameters[,"kappa_y"]==0),"mean_d"],M[(parameters[,"kappa_y"]==0),"stdv_d"])
plot(M[(parameters[,"kappa_y"]==0)&(parameters[,"chi"]==0.9),"mean_rr"],
     M[(parameters[,"kappa_y"]==0)&(parameters[,"chi"]==0.9),"avg_PD[maxH]"])

# Save results:
save(parameters,M,Model,grids,file="results/results_strategies.Rda")






