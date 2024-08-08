# ==============================================================================
# Produces figure with fiscal limits
# ==============================================================================


ColS <- ifelse(indic.colors,rgb(0,0,1),rgb(0.5,0.5,0.5))


color.area.0 <- "#11111122" # color used for the first area (weak default risk)
color.area.1 <- "#11111144" # color used for the first area (mild default risk)
color.area.2 <- "#11111166" # color used for the first area (higher default risk)

      
    
# ==============================================================================
indic_compute_ell <- 0
# ==============================================================================

# First step, construct names of files

# compute number of files:
nb.ctries <- length(list.of.ctries)
nb.files  <- nb.ctries/(2*max.nb.rows.per.figure) # because two columns
nb.files <- ifelse(nb.files==round(nb.files),nb.files,trunc(nb.files)+1)

list.names.of.file <- NULL
list.nb.ctries <- list()
ctries.already.done <- 0
for(i in 1:nb.files){
  ctries.on.chart.i <- (ctries.already.done+1):
    (ctries.already.done+min(nb.ctries-ctries.already.done,2*max.nb.rows.per.figure))
  ctries.already.done <- ctries.already.done + length(ctries.on.chart.i)
  list.nb.ctries[[i]] <- ctries.on.chart.i
  ctries.string <- list.of.ctries[ctries.on.chart.i[1]]
  for(j in 2:length(ctries.on.chart.i)){
    ctries.string <- paste(ctries.string,list.of.ctries[ctries.on.chart.i[j]],  sep="-")
  }
  list.names.of.file <- c(list.names.of.file,paste("FigureFL", ctries.string,  sep="_"))
}


for(j in 1:length(list.names.of.file)){
  
  list_of_countries <- list.nb.ctries[[j]]
  nb.ctries.on.chart.j <- length(list_of_countries)
  
  FILE = paste("figures/",list.names.of.file[j],suffix,".pdf",sep="")
  
  pdf(file=FILE,pointsize=11,width=8, height=1.3*nb.ctries.on.chart.j)
  
  
  par(mfrow=c(round(nb.ctries.on.chart.j/2),2)) 
  
  par(plt=c(.1,.95,.2,.8))
  
  counter.ctry <- 0
  
  for(ctry in list.of.ctries[list_of_countries]){
    counter.ctry <- counter.ctry + 1
    
    path.results <- paste("results/res_estim_",ctry,suffix,".Rdat",sep="")
    load(file=path.results)
    
    vec.dates <- DATA$DATE
    
    # Construct estimated model:
    model.est <- Theta2Model(full_theta_est,model_ini_sol,bounds)
    model.sol.est <- solve_model(model.est,indic_compute_beta = 1)
    
    n_w <- dim(model.sol.est$Phi)[1]
    n_x <- dim(model.sol.est$Phi_x)[1]
    
    resEKF_est <- EKF_smoother(model.sol.est,list.stdv,DATA)
    
    T     <- dim(DATA$CDS)[1]
    all.X <- matrix(0,T,n_x)
    all.X[,1:n_w] <- t(resEKF_est$W.smoothed)
    all.X[,n_w+1] <- DATA$r
    all.X[,n_w+2] <- DATA$d
    all.X[,n_w+3] <- c(DATA$d[1],DATA$d[1:(T-1)])
    
    elld <- DATA$d/4
    
    # Compute fiscal limits:
    ell1 <-
      compute.debt.limits(model.sol.est,all.X,Proba=.01,h=4)/model.sol.est$freq
    ell2 <-
      compute.debt.limits(model.sol.est,all.X,Proba=.05,h=4)/model.sol.est$freq
    
    # Compute RN fiscal limits:
    ell1.RN <-
      compute.debt.limits(model.sol.est,all.X,Proba=.01,h=4,indic.P = 0)/model.sol.est$freq
    ell2.RN <-
      compute.debt.limits(model.sol.est,all.X,Proba=.05,h=4,indic.P = 0)/model.sol.est$freq
    
    min.y <- .5*min(elld,ell1,ell2)
    if(counter.ctry==1){
      max.y <- min.y + 1.5*(max(elld,ell1) - min.y)
    }else{
      max.y <- min.y + 1.3*(max(elld,ell1) - min.y)
    }
    
    plot(vec.dates,elld,type="l",
         las=1,lwd=2,
         ylim=c(min.y,max.y),
         xlab="",ylab="",col=ColS,
         cex.main=parcex, cex.lab=parcex, cex.axis=parcex)
    title(main=paste(DATA$name.of.country$long.name,sep=""),
          cex.main=parcex)
    
    if(indic_compute_ell==1){
      # ========================================================================
      # Compute alternative ell:
      s_star <- model.sol.est$mu_star +
        t(resEKF_est$W.updated) %*% (model.sol.est$sigma_star)
      s_noBetaD <- t(resEKF_est$W.updated) %*% (model.sol.est$sigma_s)
      d_test <- (model.sol.est$d_star + (s_star-s_noBetaD)/model.sol.est$beta)/4
      # Compute Ell:
      model.test <- model.sol.est
      model.test$iota <- 10000
      #model.test$delta <- .985
      res.ell <- compute_AB_ell(model.test,100)
      ell <- exp(all.X[,1:n_w] %*% res.ell$a.l + res.ell$b.l)/model.sol.est$freq
      lines(vec.dates,ell,col="green",lwd=2,lty=2)
      lines(vec.dates,d_test,col="green",lwd=2,lty=1)
      # ========================================================================
    }
    
    polygon(c(vec.dates,rev(vec.dates)),c(ell1,rev(ell2)),
            col=color.area.0,border=NA)
    polygon(c(vec.dates,rev(vec.dates)),c(ell2,rep(1000,length(elld))),
            col=color.area.2,border=NA)
    
    #main = paste(model.sol$name.of.country$long.name,sep=""))
    lines(vec.dates,ell1,lty=3,col="black",lwd=3)
    lines(vec.dates,ell1.RN,lty=2,col="red",lwd=1)
    # lines(vec.dates,exp.adj.Ell.2,lty=3,col="red",lwd=3)
    abline(h=0,col="grey")
    
    list.of.dates <- DATA$list.of.dates
    for(dddate in list.of.dates){
      abline(v=as.Date(dddate,"%d/%m/%Y"),col="black",lwd=1)
    }
    
    if(counter.ctry==1){
      legend("topleft", 
             c("Debt-to-GDP",
               "Fiscal limit",
               "Risk-neutral fiscal limit"),
             lty=c(1,3,2), # gives the legend appropriate symbols (lines)       
             lwd=c(3,3,1), # line width
             col=c(ColS,"black","red"),
             bg="white",
             seg.len = 2,
             cex=1)
    }
    
  }
  
  dev.off()

}





