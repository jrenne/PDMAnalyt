# =================================================================
# Produces figure with fiscal limits
# =================================================================


case.suffix  <- c("","","_noCDS")
case.noCDS   <- c(0,1,1)
case.col     <- c("black","red","blue")


case.suffix  <- c("","")
case.noCDS   <- c(0,1)
case.col     <- c("black","red")


nb.ctries <- length(list.of.ctries)

FILE = "figures/FigureFL_noCDS.pdf"

pdf(file=FILE,pointsize=11,width=8, height=6)

par(mfrow=c(2,2)) 

par(plt=c(.1,.95,.2,.8))

counter.ctry <- 0

for(ctry in list.of.ctries){
  counter.ctry <- counter.ctry + 1
  
  print(paste("==========",ctry,"=========="))
  
  for(case in 1:length(case.suffix)){

    suffix       <- case.suffix[case]
    indic.no.CDS <- case.noCDS[case]
    col          <- case.col[case]
    
    path.results <- paste("results/res_estim_",ctry,suffix,".Rdat",sep="")
    load(file=path.results)
    
    if(indic.no.CDS==1){
      DATA$CDS <- NaN * DATA$CDS
    }
    
    vec.dates <- DATA$DATE
    
    # Construct estimated model:
    model.est     <- Theta2Model(full_theta_est,model_ini_sol,bounds)
    model.sol.est <- solve_model(model.est,indic_compute_beta = 1)
    
    n_w <- dim(model.sol.est$Phi)[1]
    n_x <- dim(model.sol.est$Phi_x)[1]
    
    resEKF_est <- EKF_smoother(model.sol.est,list.stdv,DATA)
    
    logl0 <- resEKF_est$loglik
    
    T     <- dim(DATA$CDS)[1]
    all.X <- matrix(0,T,n_x)
    all.X[,1:n_w] <- t(resEKF_est$W.smoothed)
    all.X[,n_w+1] <- DATA$r
    all.X[,n_w+2] <- DATA$d
    all.X[,n_w+3] <- c(DATA$d[1],DATA$d[1:(T-1)])
    
    elld <- DATA$d/4
    
    ELL1 <- compute.debt.limits(model.sol.est,all.X,Proba=.01,h=4,MIN.d=-4)/model.sol.est$freq

    # Compute std dev of FL estimates: -----------------------------------------
    # Compute derivatives wrt w:
    Epsilon <- .001
    all.d.ell1 <- NULL
    for(jjj in 1:n_w){
      all.X.pert <- all.X
      all.X.pert[,jjj] <- all.X.pert[,jjj] + Epsilon
      ell1.pert <-
        compute.debt.limits(model.sol.est,all.X.pert,Proba=.01,h=4,MIN.d=-4)/model.sol.est$freq
      d.ell1 <- (ell1.pert - ELL1)/Epsilon
      all.d.ell1 <- cbind(all.d.ell1,d.ell1)
    }
    # Compute variance of ell1 estimates:
    stdv.ell1 <- matrix(0,dim(all.X)[1],1)
    for(t in 1:dim(all.X)[1]){
      aux <- matrix(all.d.ell1[t,],ncol=1)
      var.ell1.t <- t(aux) %*% resEKF_est$P.smoothed[[t]] %*% aux
      stdv.ell1[t] <- sqrt(var.ell1.t)
    }
    # --------------------------------------------------------------------------
    
    # min.y <- -.3
    # if(counter.ctry==1){
    #   max.y <- max(min.y + 1.5*(max(ELL1-elld) - min.y),.5)
    # }else{
    #   max.y <- max(min.y + 1.2*(max(ELL1-elld) - min.y),.5)
    # }
    
    min.y <- 0
    
    if(ctry=="US"){max.y <- 3}    
    if(ctry=="UK"){max.y <- 2}    
    if(ctry=="CA"){max.y <- 2}    
    if(ctry=="JP"){max.y <- 4}    
    if(ctry=="BR"){max.y <- 1.5}    
    if(ctry=="RU"){max.y <- .4}    
    if(ctry=="IN"){max.y <- 1.6}    
    if(ctry=="CN"){max.y <- 1.6}    
    
    if(case == 1){
      plot(vec.dates,elld,type="l",
           las=1,
           ylim=c(min.y,max.y),
           xlab="",ylab="",
           #col=ColS,
           col="grey",lwd=2,lty=1,
           cex.main=parcex, cex.lab=parcex, cex.axis=parcex)
      title(main=paste(DATA$name.of.country$long.name,sep=""),
            cex.main=parcex)
    }
    lines(vec.dates,ELL1,lwd=2,lty=1,
          col=col)

    # Add confidence intervals:
    lines(vec.dates,ELL1-qnorm(.995)*stdv.ell1,lty=2,col=col,lwd=1)
    lines(vec.dates,ELL1+qnorm(.995)*stdv.ell1,lty=2,col=col,lwd=1)
    
  }
  
  if(ctry == "JP"){
    legend("bottomleft", 
           c("Debt-to-GDP",
             "",
             "Fiscal limit (baseline estim.)",
             "Fiscal limit (baseline estim., no CDS data)"),
           lty=c(1,1,1,1), # gives the legend appropriate symbols (lines)       
           lwd=c(2,2,2,2), # line width
           col=c("grey","white","black","red"),
           bg="white",
           seg.len = 2,
           cex=1)
  }
  
  
  
}

dev.off()


