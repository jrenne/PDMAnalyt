

# Create stylized models

low_pi <- -.01
med_pi <- .03
hig_pi <- .07

low_y <- -.01
med_y <- .02
hig_y <- .05

mu_y  <- matrix(c(low_y,med_y,hig_y),ncol=1)

# Demand:
mu_pi_demand  <- matrix(c(low_pi,med_pi,hig_pi),ncol=1)
# Supply:
mu_pi_supply  <- matrix(c(hig_pi,med_pi,low_pi),ncol=1)

rho <- .8
Omega <- diag(rep(rho,3))
Omega[1,2] <- 1 - rho
Omega[3,2] <- 1 - rho
Omega[2,1] <- (1 - rho)/2
Omega[2,3] <- (1 - rho)/2

Model <- list(  mu_pi = NaN,
                mu_y = mu_y,
                Omega = Omega,
                gamma = 10,
                delta = .99,
                chi = .8,
                mu_eta = 1*mu_y,
                beta = .1,
                alpha = .3,
                d_star =.8,
                s_star =.05,
                RR= .5,
                nu_pi = 0,
                nu_y = 0)

Model_Demand <- Model
Model_Demand$mu_pi <- mu_pi_demand

Model_Supply <- Model
Model_Supply$mu_pi <- mu_pi_supply

FILE = paste("figures/Figure_expected_returns_DemaSupp.pdf",sep="")
pdf(file=FILE, pointsize=10, width=9, height=5)

par(mfrow=c(1,2))
par(plt=c(.12,.95,.1,.85))

for(regime in c("Demand","Supply")){
  
  eval(parse(text = gsub(" ","",paste("Model <- Model_",regime,sep=""))))
  
  RES <- prepare_returns_yds(Model,maxH)
  
  plot(1:maxH,rep(0,maxH),las=1,
       xlab="",ylab="",type="l",ylim=c(0.0,.11),lwd=2,col="white",
       main=paste(ifelse(regime=="Demand","(a) ","(b) "),regime,"-driven economy",sep=""))
  
  lower_bound <- RES$avg_annual_nominal_returns - 1.0 * RES$Std_annual_nominal_returns
  upper_bound <- RES$avg_annual_nominal_returns + 1.0 * RES$Std_annual_nominal_returns
  polygon(c(1:maxH,rev(1:maxH)),c(lower_bound,rev(upper_bound)),border = NaN,
          col="#66AAAA44")
  lines(1:maxH,c(RES$avg_annual_nominal_returns),col="blue",lwd=2)
  
  lines(1:maxH,RES$avg_annual_TIPS_returns,col="red",lwd=2)
  
  lower_bound <- RES$avg_annual_GDPLB_returns - 1.0 * RES$Std_annual_GDPLB_returns
  upper_bound <- RES$avg_annual_GDPLB_returns + 1.0 * RES$Std_annual_GDPLB_returns
  polygon(c(1:maxH,rev(1:maxH)),c(lower_bound,rev(upper_bound)),border = NaN,
          col="#AAAAAA44")
  lines(1:maxH,RES$avg_annual_GDPLB_returns,col="black",lwd=2,lty=2)
  
  grid()
  
  if(regime=="Supply"){
    legend("topright", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend
           c("Nominal bonds","TIPS","GDP-LB"),
           lty=c(1,1,2), # gives the legend appropriate symbols (lines)
           lwd=c(2,2,2), # line width
           col=c("blue","red","black"),
           bg="white",
           #pch=c(3,NaN,NaN),
           seg.len = 3,
           cex=1,title = "Annualized expected nominal returns")
    
  }
}

dev.off()

# Debt issuance strategies:
Model_Demand$chi    <- .8
Model_Demand$mu_eta <- 1*Model_Demand$mu_y
Model_Demand$beta   <- .1
Model_Demand$alpha  <- .3
Model_Demand$d_star <- .8
Model_Demand$s_star <- .05
Model_Demand$RR     <- .5
Model_Demand$nu_pi  <- 0
Model_Demand$nu_y   <- 0

Models <- prepare_and_solve_3(Model_Supply,
                              all_d,all_rr,all_eps,proba_eps)

strat_nominal <- run_strategy(Models$Model_solved_nominal,maxH=10)
strat_TIPS    <- run_strategy(Models$Model_solved_TIPS,maxH=10)
strat_GDPLB   <- run_strategy(Models$Model_solved_GDPLB,maxH=10)

M <- rbind(c(strat_nominal$mean_d,strat_TIPS$mean_d,strat_GDPLB$mean_d),
           c(strat_nominal$stdv_d,strat_TIPS$stdv_d,strat_GDPLB$stdv_d),
           c(strat_nominal$mean_rr,strat_TIPS$mean_rr,strat_GDPLB$mean_rr),
           c(strat_nominal$stdv_rr,strat_TIPS$stdv_rr,strat_GDPLB$stdv_rr),
           c(strat_nominal$mean_Delta_d,strat_TIPS$mean_Delta_d,strat_GDPLB$mean_Delta_d),
           c(strat_nominal$stdv_Delta_d,strat_TIPS$stdv_Delta_d,strat_GDPLB$stdv_Delta_d))

colnames(M) <- c("Nominal","TIPS","GDPLB")
rownames(M) <- c("Mean d","Stdv d",
                 "Mean r","Stdv r",
                 "Mean D(d)","Stdv D(d)")

print(M)

