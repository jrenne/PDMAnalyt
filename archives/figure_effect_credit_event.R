# ------------------------------------------------------------------------------
# Prepare plot showing effect of credit risk on nominal yield curves
# (in the context of the demand/supply exercise)
# ------------------------------------------------------------------------------

FILE = paste("figures/Figure_yds_curve_wCredit_DemaSupp.pdf",sep="")
pdf(file=FILE, pointsize=10, width=9, height=5)

par(mfrow=c(1,2))
par(plt=c(.12,.95,.1,.85))

for(regime in c("Demand","Supply")){
  
  eval(parse(text = gsub(" ","",paste("Model <- Model_",regime,sep=""))))
  
  # Nominal bonds:
  Model$kappa_pi <- 0
  Model$kappa_y  <- 0
  
  # Without credit risk: -------------------------------------------------------
  
  RES <- prepare_returns_yds(Model,maxH)
  
  ylim <- c(min(RES$avg_annual_nominal_returns)-.005,
            max(RES$avg_annual_nominal_returns)+.01)
  
  plot(1:maxH,rep(0,maxH),las=1,
       xlab="",ylab="",type="l",ylim=ylim,lwd=2,col="white",
       main=paste(ifelse(regime=="Demand","(a) ","(b) "),regime,"-driven economy",sep=""))
  
  lines(1:maxH,c(RES$avg_annual_nominal_returns),col="blue",lwd=2)
  
  
  # With credit risk: ----------------------------------------------------------
  
  # First case, nu_y = 0
  Model$nu_y <- 0
  Model_solved <- solve_ToyModel(Model,
                                 grids,nb_iter,
                                 nb_iter_sdf)
  res_bond_prices <- compute_bond_prices(Model_solved,maxH=10,nb_iter_sdf)
  p <- compute_uncond_distri(Model_solved$indicators_x,
                             Model_solved$Probas,nb_iter = 2000)
  avg_yds <- t(p) %*% res_bond_prices$all_rth
  lines(1:maxH,c(avg_yds),col="blue",lwd=2,lty=2)
  
  # Second case, nu_y < 0
  Model$nu_y <- - 0.1
  Model_solved <- solve_ToyModel(Model,
                                 grids,nb_iter,
                                 nb_iter_sdf)
  res_bond_prices <- compute_bond_prices(Model_solved,maxH=10,nb_iter_sdf)
  p <- compute_uncond_distri(Model_solved$indicators_x,
                             Model_solved$Probas,nb_iter = 2000)
  avg_yds <- t(p) %*% res_bond_prices$all_rth
  lines(1:maxH,c(avg_yds),col="blue",lwd=2,lty=3)
  
  grid()
  
  if(regime=="Supply"){
    legend("topright", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend
           c("no credit risk",
             "with credit risk, no macro impact of default",
             "with credit risk, with macro impact of default"),
           lty=c(1,2,3), # gives the legend appropriate symbols (lines)
           lwd=c(2,2,2), # line width
           col=c("blue","blue","blue"),
           bg="white",
           #pch=c(3,NaN,NaN),
           seg.len = 3,
           cex=1)
  }
}

dev.off()

