

# FILE = paste("figures/Figure_Toy_compstat.pdf",sep="")
# #pdf(file=FILE, pointsize=8, width=7, height=4)
# pdf(file=FILE, pointsize=9, width=5, height=3)

par(mfrow=c(1,1))
par(plt=c(.13,.95,.2,.95))

# res0_notRcpp <- solve_ToyModel_notRcpp(all_d,all_rr,all_eps,proba_eps,
#                                        Model,nb_iter = nb_iter)

res0 <- solve_ToyModel(all_d,all_rr,all_eps,proba_eps,
                       Model,nb_iter = nb_iter)


PD <- compute_proba_def(maxH=maxH,
                        indicators_x = res0$indicators_x,
                        all_proba_def = res0$all_proba_def,
                        Probas = res0$Probas)

res_LTprices <- compute_LTRF_bond_prices(Model,maxH)
avgLTyields <- c(t(res0$stat_distri) %*% res_LTprices$all_LT_rth)
plot(avgLTyields)
print(paste("Slope of yield curve: ",
            round(10000*(avgLTyields[length(avgLTyields)]-avgLTyields[1]),2),
            " basis points",sep=""))
print(cbind(Model$mu_pi,Model$mu_y))

stop()

# PD_notRcpp <- compute_proba_def_notRcpp(maxH=10,
#                                         indicators_x = res0$indicators_x,
#                                         all_proba_def = res0$all_proba_def,
#                                         Probas = res0$Probas)

par(mfrow=c(1,2))
v <- 15
v <- which(all_d==d_star)
Ts <- which((res0$d==res0$d_1)&(res0$d==all_d[v])&(res0$rr==all_rr[v]))
t <- Ts[1]

plot(PD[t,])


stop()

# # ============================================================================
# # Figure illustrating the effect of RR on q and q0:
# # ============================================================================
# 
# all.values.d <- c(.6,1)
# 
# all.RR <- seq(0,1,by=.05)
# 
# all.q  <- matrix(NaN,length(all.RR),length(all.values.d))
# all.q0 <- matrix(NaN,length(all.RR),length(all.values.d))
# 
# indic.value.d <- 0
# for(value.d in all.values.d){
#   indic.value.d <- indic.value.d + 1
#   
#   value.d_1 <- value.d
#   all_x <- cbind(res0$d,res0$d_1,res0$rr)
#   # Determine position of closest state:
#   indics <- which(
#     ((all_x[,1]-value.d)^2==min((all_x[,1]-value.d)^2))&
#       ((all_x[,2]-value.d_1)^2==min((all_x[,2]-value.d_1)^2)&
#          (all_x[,3]-value.rr)^2==min((all_x[,3]-value.rr)^2)))[1]
#   
#   indic.RR <- 0
#   for(RR in all.RR){
#     indic.RR <-   indic.RR + 1
#     Model_modif <- Model
#     Model_modif$RR <- RR
#     res_modif <- solve_ToyModel(all_d,all_rr,all_eps,
#                                 Model_modif,nb_iter)
#     
#     all.q[indic.RR,indic.value.d]   <- 100*res_modif$q[indics]
#     all.q0[indic.RR,indic.value.d]  <- 100*res_modif$q0[indics]
#   }
#   
# }
# 
# names.4.legend <- NULL
# for(i in 1:length(all.values.d)){
#   if(i==1){
#     plot(all.RR,all.q[,i],type="l",lwd=2,ylim=c(-1,5),
#          xlab="Recovery rate (RR)",ylab="yields, in percent",las=1,
#          main="Panel (a)")
#   }else{
#     lines(all.RR,all.q[,i],lty=i,lwd=2,col="black")
#   }
#   lines(all.RR,all.q0[,i],lty=i,lwd=2,col="grey")
#   abline(v=Model$RR,col="grey",lty=3)
#   abline(h=100*res0$q[indics],col="grey",lty=3)
#   names.4.legend <- c(names.4.legend,
#                       toString(round(all.values.d[i],2)))
# }
# 
# arrows(.7,0,Model$RR,-1,length = .1)
# text(.7,0,expression(paste("RR = exp(-",gamma,"b)",sep="")),pos=4)
# 
# legend("topright",
#        names.4.legend,
#        lty=c(1,2,3), # gives the legend appropriate symbols (lines)       
#        lwd=2, # line width
#        col=c("black","black"), # gives the legend lines the correct color and width
#        pt.bg=c(1,1),
#        pch = NaN,#symbols,
#        pt.cex = c(1,1),
#        bg="white",
#        seg.len = 3,
#        cex=1,
#        title="Values of debt-to-GDP:"
# )
# 
# legend("bottomleft",
#        c("Sovereign bond yield","Risk-free yield"),
#        lty=c(1), # gives the legend appropriate symbols (lines)       
#        lwd=4, # line width
#        col=c("black","grey"), # gives the legend lines the correct color and width
#        bg="white",
#        seg.len = 1,
#        cex=1,
#        title="Type of yield:"
# )


# ==============================================================================
# Figure illustrating the effect of d on q and q0:
# ==============================================================================

all_x <- cbind(res0$d,res0$d_1,res0$rr)

all.values.RR <- c(.3,.6) # for Panel (b)

value.rr     <- mean(res0$q)*mean(res0$d)

# Determine position of closest state:
indics <- which(
  (all_x[,1]==all_x[,2])&
    (all_x[,3]-value.rr)^2 == min((all_x[,3]-value.rr)^2))

all_q  <- 100*res0$q[indics,]
all_q0 <- 100*res0$q0[indics,]

names.4.legend <- toString(round(Model$RR,2))

eval(parse(text = gsub(" ","",
                       paste("names.4.legend <- expression(paste('exp(-',gamma,b[y],') = ",toString(round(Model$RR,2)),"',sep=''))")
)))


for(i in 1:length(all.values.RR)){
  Model_modif <- Model
  Model_modif$RR <- all.values.RR[i]
  res_modif <- solve_ToyModel(all_d,all_rr,all_eps,
                              Model_modif,nb_iter)
  all_q  <- cbind(all_q,100*res_modif$q[indics,])
  all_q0 <- cbind(all_q0,100*res_modif$q0[indics,])
  names.4.legend <- c(names.4.legend,
                      toString(round(Model_modif$RR,2)))
}

plot(all_d,all_q[,1],lwd=2,las=1,
     ylim=c(-1,5),type="l",
     xlab="Debt-to-GDP ratio",ylab="yields, in percent")

lines(all_d,all_q0[,1],col="grey",lwd=2)
for(i in 1:length(all.values.RR)){
  lines(all_d,all_q[,1+i],col="black",lwd=2,lty=i+1)
  lines(all_d,all_q0[,1+i],col="grey",lwd=2,lty=i+1)
}

grid()

#abline(v=all.values.d[1],col="grey",lty=3)
#abline(v=all.values.d[2],col="grey",lty=3)

legend("bottomleft",
       names.4.legend,
       lty=c(1,2,3), # gives the legend appropriate symbols (lines)
       lwd=2, # line width
       col=c("black","black"), # gives the legend lines the correct color and width
       pt.bg=c(1,1),
       pch = NaN,#symbols,
       pt.cex = c(1,1),
       bg="white",
       seg.len = 3,
       cex=1,
       title="Values of RR:"
)

legend("topleft",
       c("Sovereign bond yield","Risk-free yield"),
       lty=c(1), # gives the legend appropriate symbols (lines)
       lwd=4, # line width
       col=c("black","grey"), # gives the legend lines the correct color and width
       bg="white",
       seg.len = 1,
       cex=1,
       title="Type of yield:"
)

dev.off()





# ==============================================================================
# Figure illustrating the concept of credit risk premium:
# ==============================================================================

FILE = paste("figures/Figure_Toy_CDSs.pdf",sep="")
pdf(file=FILE, pointsize=10, width=8, height=4)

res0 <- solve_ToyModel(all_d,all_rr,all_eps,
                       Model,nb_iter)

# Compute bond prices and PDs:
res_PDs <- compute_PDs_bond_prices(all_d,all_rr,all_eps,
                                   Model,nb_iter,
                                   maxH)

# Use Shadow-rate approximation:
X <- cbind(res0$d,res0$d_1,res0$rr,0)
res_PDs_SR <- compute_PDs_bond_prices_SRapprox(X,Model,maxH)

all.values.d <- c(.3,.6,.9)
value.rr <- 0.02
par(mfrow=c(1,length(all.values.d)))

for(value.d in all.values.d){
  
  value.d_1 <- value.d
  
  # Determine position of closest state:
  indics <- which(
    ((all_x[,1]-value.d)^2==min((all_x[,1]-value.d)^2))&
      ((all_x[,2]-value.d_1)^2==min((all_x[,2]-value.d_1)^2)&
         (all_x[,3]-value.rr)^2==min((all_x[,3]-value.rr)^2)))[1]
  
  digitCDS    <- 1/(1-Model$RR)*(res_PDs$B0th-res_PDs$Bth)/res_PDs$B0th
  digitCDS_SR <- 1/(1-Model$RR)*(res_PDs_SR$B0th-res_PDs_SR$Bth)/res_PDs_SR$B0th
  
  plot(1:maxH,100*res_PDs$all.prob.def[indics,],type="l",lwd=2,las=1,
       ylim=c(0,18),
       xlab="Maturity",ylab="Probabilities, in percent",
       main=paste("Debt-to-GDP: ",toString(100*value.d),"%",sep=""))
  lines(1:maxH,100*digitCDS[indics,],lwd=2,lty=3)
  points(1:maxH,100*res_PDs_SR$all.prob.def[indics,],pch=6,col="dark grey",lwd=2)
  points(1:maxH,100*digitCDS_SR[indics,],pch=2,col="dark grey",lwd=2)
  
  grid()
  
  if(value.d==all.values.d[1]){
    legend("topleft",
           c("Default Proba.","Digital CDS"),
           lty=c(1,3), # gives the legend appropriate symbols (lines)       
           lwd=3, # line width
           col=c("black","black"), # gives the legend lines the correct color and width
           bg="white",
           seg.len = 2,
           cex=1,
           title="Numerical solution:"
    )
    legend("topright",
           c("Default Proba.","Digital CDS"),
           lty=NaN, # gives the legend appropriate symbols (lines)       
           lwd=2, # line width
           col=c("dark grey"), # gives the legend lines the correct color and width
           pch=c(6,2),
           bg="white",
           seg.len = 2,
           cex=1,
           title="Approximation:"
    )
    
  }
}

dev.off()




# ==============================================================================
# Impulse response functions:
# ==============================================================================


X <- matrix(0,2,4) # faked values, just to run next line
res_PDs_SR <- compute_PDs_bond_prices_SRapprox(X,Model,maxH)

Mu_x    <- res_PDs_SR$Mu_x
Phi_x   <- res_PDs_SR$Phi_x
Sigma_x <- res_PDs_SR$Sigma_x

EX <- solve(diag(4)-Phi_x) %*% Mu_x

maxH.IRF <- 26

shock <- -1
x <- Sigma_x * shock
all_x <- x
for(i in 1:(maxH.IRF-1)){
  x <- Phi_x %*% x
  all_x <- cbind(all_x,x)
}

# Compute proba of default for different starting values of debt


FILE = paste("figures/Figure_Toy_IRFs.pdf",sep="")
pdf(file=FILE, pointsize=10, width=7, height=4)
par(mfrow=c(2,2))
par(plt=c(.16,.95,.15,.8))

# Debt:
plot(0:(maxH.IRF-1),100*all_x[1,],type="l",lwd=2,las=1,col="dark grey",
     xlab="Time after shock",ylab="Debt, in percent of GDP",
     main="(a) Debt-to-GDP")
grid()

# Surplus:
plot(0:(maxH.IRF-1),100*all_x[4,],type="l",lwd=2,las=1,col="dark grey",
     xlab="Time after shock",ylab="Surplus, in percent of GDP",
     main="(b) Surplus-to-GDP")
grid()

# Proba of default and spreads:

maturity <- 5

PDs4values_d <- matrix(NaN,maxH.IRF,length(all.values.d))
spd4values_d <- matrix(NaN,maxH.IRF,length(all.values.d))

iii <- 0
for(value.d in all.values.d){
  iii <- iii+1
  
  # baseline:
  x_baseline <- EX
  x_baseline[1:2] <- value.d
  x_shock <- x_baseline + shock*Sigma_x
  all_x_baseline <- x_baseline
  all_x_shock    <- x_shock
  for(i in 1:(maxH.IRF-1)){
    x_baseline <- Mu_x + Phi_x %*% x_baseline
    x_shock    <- Mu_x + Phi_x %*% x_shock
    all_x_baseline <- cbind(all_x_baseline,x_baseline)
    all_x_shock    <- cbind(all_x_shock,x_shock)
  }
  
  res_baseline <- compute_PDs_bond_prices_SRapprox(t(all_x_baseline),Model,maturity+5)
  res_shock    <- compute_PDs_bond_prices_SRapprox(t(all_x_shock),Model,maturity+5)
  
  PDs4values_d[,iii] <- res_shock$all.prob.def[,maturity] - 
    res_baseline$all.prob.def[,maturity]
  spd4values_d[,iii] <- (res_shock$yth[,maturity]-res_shock$y0th[,maturity]) - 
    (res_baseline$yth[,maturity]-res_baseline$y0th[,maturity])
  
}


PDs4values_d <- PDs4values_d*100
spd4values_d <- spd4values_d*10000

# PDs:
names.4.legend <- NULL
for(i in 1:length(all.values.d)){
  if(i==1){
    plot(0:(maxH.IRF-1),PDs4values_d[,i],lty=i,type="l",lwd=2,las=1,
         ylim=c(0,max(PDs4values_d)),
         xlab="Time after shock",ylab="Probability, in percent",
         main=paste("(c) Default probability, ",toString(maturity),"-year horizon",sep=""))
  }else{
    lines(0:(maxH.IRF-1),PDs4values_d[,i],lty=i,lwd=2)
  }
  names.4.legend <- c(names.4.legend,
                      toString(round(all.values.d[i],2)))
}
grid()

legend("topright",
       names.4.legend,
       lty=c(1,2,3), # gives the legend appropriate symbols (lines)       
       lwd=2, # line width
       col=c("black","black"), # gives the legend lines the correct color and width
       pt.bg=c(1,1),
       pch = NaN,#symbols,
       pt.cex = c(1,1),
       bg="white",
       seg.len = 3,
       cex=1,
       title="Initial debt-to-GDP:"
)


# Spreads:
for(i in 1:length(all.values.d)){
  if(i==1){
    plot(0:(maxH.IRF-1),spd4values_d[,i],lty=i,type="l",lwd=2,las=1,
         ylim=c(0,max(spd4values_d)),
         xlab="Time after shock",ylab="Spreads, in basis points",
         main=paste("(d) Bond spread, ",toString(maturity),"-year maturity",sep=""))
  }else{
    lines(0:(maxH.IRF-1),spd4values_d[,i],lty=i,lwd=2)
  }
}
grid()


dev.off()



# 
# 
# # ==============================================================================
# # Figure illustrating the effect of d on q and q0 and on PDs:
# # ==============================================================================
# 
# 
# FILE = paste("figures/Figure_Toy_compstat2.pdf",sep="")
# #pdf(file=FILE, pointsize=8, width=7, height=4)
# pdf(file=FILE, pointsize=8, width=7, height=3)
# 
# par(mfrow=c(1,2))
# par(plt=c(.11,.95,.2,.85))
# 
# # Solve baseline model:
# res0 <- solve_ToyModel(all_d,all_rr,all_eps,
#                        Model,nb_iter)
# 
# # Compute probabilities od default:
# res_PDs <- compute_PDs_bond_prices(all_d,all_rr,all_eps,
#                                    Model,nb_iter,maxH)
# 
# value.rr     <- mean(res0$q)*mean(res0$d)
# 
# all_x <- cbind(res0$d,res0$d_1,res0$rr)
# 
# all.values.RR <- c(.3,.6)
# 
# value.rr     <- mean(res0$q)*mean(res0$d)
# 
# # Determine position of closest state:
# indics <- which(
#   (all_x[,1]==all_x[,2])&
#     (all_x[,3]-value.rr)^2==min((all_x[,3]-value.rr)^2))
# 
# all_q  <- 100*res0$q[indics,]
# all_q0 <- 100*res0$q0[indics,]
# all_PD <- 100*res_PDs$all.prob.def[indics,maxH]
# 
# names.4.legend <- toString(round(Model$RR,2))
# 
# eval(parse(text = gsub(" ","",
#                        paste("names.4.legend <- expression(paste('exp(-',gamma,' b) = ",toString(round(Model$RR,2)),"',sep=''))")
# )))
# 
# 
# for(i in 1:length(all.values.RR)){
#   Model_modif <- Model
#   Model_modif$RR <- all.values.RR[i]
#   res_modif <- solve_ToyModel(all_d,all_rr,all_eps,
#                               Model_modif,nb_iter)
#   all_q  <- cbind(all_q,100*res_modif$q[indics,])
#   all_q0 <- cbind(all_q0,100*res_modif$q0[indics,])
#   
#   res_PDs_modif <- compute_PDs_bond_prices(all_d,all_rr,all_eps,
#                                            Model_modif,nb_iter,maxH)
#   all_PD <- cbind(all_PD,100*res_PDs_modif$all.prob.def[indics,maxH])
#   
#   names.4.legend <- c(names.4.legend,
#                       toString(round(Model_modif$RR,2)))
# }
# 
# plot(all_d,all_q[,1],lwd=2,las=1,
#      ylim=c(-1,5),type="l",
#      main = "(a) Perpetuity's yield-to-maturity",
#      xlab="Debt-to-GDP ratio",ylab="yields, in percent")
# 
# lines(all_d,all_q0[,1],col="grey",lwd=2)
# for(i in 1:length(all.values.RR)){
#   lines(all_d,all_q[,1+i],col="black",lwd=2,lty=i+1)
#   #lines(all_d,all_q0[,1+i],col="grey",lwd=2,lty=i+1)
# }
# 
# grid()
# 
# #abline(v=all.values.d[1],col="grey",lty=3)
# #abline(v=all.values.d[2],col="grey",lty=3)
# 
# legend("bottomleft",
#        names.4.legend,
#        lty=c(1,2,3), # gives the legend appropriate symbols (lines)
#        lwd=2, # line width
#        col=c("black","black"), # gives the legend lines the correct color and width
#        pt.bg=c(1,1),
#        pch = NaN,#symbols,
#        pt.cex = c(1,1),
#        bg="white",
#        seg.len = 3,
#        cex=1,
#        title="Values of RR:"
# )
# 
# legend("topleft",
#        c("Sovereign bond yield","Risk-free yield"),
#        lty=c(1), # gives the legend appropriate symbols (lines)
#        lwd=4, # line width
#        col=c("black","grey"), # gives the legend lines the correct color and width
#        bg="white",
#        seg.len = 1,
#        cex=1,
#        title="Type of yield:"
# )
# 
# plot(all_d,all_PD[,1],lwd=2,las=1,
#      ylim=c(-1,20),type="l",
#      main = "(b) 10-year probability of default",
#      xlab="Debt-to-GDP ratio",ylab="yields, in percent")
# for(i in 1:length(all.values.RR)){
#   lines(all_d,all_PD[,1+i],col="black",lwd=2,lty=i+1)
#   #lines(all_d,all_q0[,1+i],col="grey",lwd=2,lty=i+1)
# }
# 
# grid()
# 
# legend("topleft",
#        names.4.legend,
#        lty=c(1,2,3), # gives the legend appropriate symbols (lines)
#        lwd=2, # line width
#        col=c("black","black"), # gives the legend lines the correct color and width
#        pt.bg=c(1,1),
#        pch = NaN,#symbols,
#        pt.cex = c(1,1),
#        bg="white",
#        seg.len = 3,
#        cex=1,
#        title="Values of RR:"
# )
# 
# dev.off()


