# =================================================================
# Produces figure with fiscal limits
# =================================================================



ColS <- ifelse(indic.colors,rgb(0,0,1),rgb(0.5,0.5,0.5))

color.area.0 <- "#11111122" # color used for the first area (weak default risk)
color.area.1 <- "#11111144" # color used for the first area (mild default risk)
color.area.2 <- "#11111166" # color used for the first area (higher default risk)


# Matrix where FL and FS of different countries will be stored:
all.FS <- data.frame(DATE=as.Date("2000-01-01"))
# Matrix where w of different countries will be stored:
all.W  <- data.frame(DATE=as.Date("2000-01-01"))

# First step, construct names of files

# compute number of files:
nb.ctries <- length(list.of.ctries)
nb.files  <- nb.ctries/(2*max.nb.rows.per.figure) # because two columns
nb.files <- ifelse(nb.files==round(nb.files),nb.files,trunc(nb.files)+1)

list.names.of.file <- NULL
list.nb.ctries <- list()
ctries.already.done <- 0
for(i in 1:nb.files){
  ctries.on.chart.i <- (ctries.already.done+1):(ctries.already.done+min(nb.ctries-ctries.already.done,2*max.nb.rows.per.figure))
  ctries.already.done <- ctries.already.done + length(ctries.on.chart.i)
  list.nb.ctries[[i]] <- ctries.on.chart.i
  ctries.string <- list.of.ctries[ctries.on.chart.i[1]]
  for(j in 2:length(ctries.on.chart.i)){
    ctries.string <- paste(ctries.string,list.of.ctries[ctries.on.chart.i[j]],  sep="-")
  }
  list.names.of.file <- c(list.names.of.file,paste("FigureFS", ctries.string,  sep="_"))
}


for(j in 1:length(list.names.of.file)){
  
  list_of_countries <- list.nb.ctries[[j]]
  nb.ctries.on.chart.j <- length(list_of_countries)
  
  FILE = paste("figures/",list.names.of.file[j],suffix,
               ifelse(indic.no.CDS,"_noCDS",""),
               ".pdf",sep="")
  
  pdf(file=FILE,pointsize=11,width=8, height=1.3*nb.ctries.on.chart.j)

  par(mfrow=c(round(nb.ctries.on.chart.j/2),2)) 
  
  par(plt=c(.1,.95,.2,.8))
  
  counter.ctry <- 0
  
  for(ctry in list.of.ctries[list_of_countries]){
    counter.ctry <- counter.ctry + 1
    
    print(paste("==========",ctry,"=========="))
    
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
    
    ELL1 <- compute.debt.limits(model.sol.est,all.X,Proba=.01,h=4)/model.sol.est$freq
    ELL2 <- compute.debt.limits(model.sol.est,all.X,Proba=.05,h=4)/model.sol.est$freq
    
    # Compute std dev of FL estimates: -----------------------------------------
    # Compute derivatives wrt w:
    Epsilon <- .001
    all.d.ell1 <- NULL
    for(jjj in 1:n_w){
      all.X.pert <- all.X
      all.X.pert[,jjj] <- all.X.pert[,jjj] + Epsilon
      ell1.pert <-
        compute.debt.limits(model.sol.est,all.X.pert,Proba=.01,h=4)/model.sol.est$freq
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
    
    if(indic.comput.stdv.param==1){
      source("outputs/run.Hamilton.R") # Compute ell1 when param drawn in asympt distri of param
      stdv.ell1.Hamilton <- apply(all.sim.ell1,1,sd)
      # Compute total stdv for ell1:
      stdv.ell1 <- sqrt(stdv.ell1^2 + stdv.ell1.Hamilton^2)
    }
    
    min.y <- -.3
    if(counter.ctry==1){
      max.y <- max(min.y + 1.5*(max(ELL1-elld) - min.y),.5)
    }else{
      max.y <- max(min.y + 1.5*(max(ELL1-elld) - min.y),.5)
    }
    
    
    
    # === Prepare outputs (csv files for FS and w factors) =====================
    # Store FL and FS results:
    eval(parse(text = gsub(" ","",
                           paste("FS <- data.frame(DATE=DATA$DATE,FS_",ctry,"=ELL1-elld,
                                 FL_",ctry,"=ELL1)",
                                 sep=""))))
    all.FS <- merge(FS,all.FS,by="DATE",all = TRUE)

    # Store w results:
    names_w <- paste("w",1:n_w,"_",ctry,sep="")
    for(iii in 1:n_w){
      eval(parse(text = gsub(" ","",
                             paste("w <- data.frame(DATE=DATA$DATE,",names_w[iii],"=all.X[,iii])",
                                   sep=""))))
      all.W <- merge(w,all.W,by="DATE",all = TRUE)
    }
    # ==========================================================================
    
    
    
    plot(vec.dates,ELL1-elld,type="l",
         las=1,lwd=3,lty=3,
         ylim=c(min.y,max.y),
         xlab="",ylab="",
         #col=ColS,
         col="red",
         cex.main=parcex, cex.lab=parcex, cex.axis=parcex)
    title(main=paste(DATA$name.of.country$long.name,sep=""),
          cex.main=parcex)
    
    polygon(c(vec.dates,rev(vec.dates)),c(ELL1-elld,rev(ELL2-elld)),
            col=color.area.0,border=NA)
    polygon(c(vec.dates,rev(vec.dates)),c(ELL2-elld,rep(1000,length(elld))),
            col=color.area.2,border=NA)
    
    abline(h=0,col="grey")

    list.of.dates <- DATA$list.of.dates
    for(dddate in list.of.dates){
      abline(v=as.Date(dddate,"%d/%m/%Y"),col="black",lwd=1)
    }
    
    # Add confidence intervals:
    lines(vec.dates,ELL1-elld-1.64*stdv.ell1,lty=1,col="red",lwd=1)
    lines(vec.dates,ELL1-elld+1.64*stdv.ell1,lty=1,col="red",lwd=1)
    
  }
  
  dev.off()
  
}

if(suffix==""){# Save FS and w factors
  # Only in the baseline estimation case
  write.csv(all.FS,"results/FS.csv")
  write.csv(all.W,"results/W.csv")
}


