# =================================================================
# Produces figure with yields
# =================================================================

ColQ <- ifelse(indic.colors,rgb(0,0,1),rgb(0.5,0.5,0.5))


indic.plot.P <- 0

# First step, construct names of files

# compute number of files:
nb.ctries <- length(list.of.ctries)
nb.files  <- nb.ctries/max.nb.rows.per.figure
nb.files <- ifelse(nb.files==round(nb.files),nb.files,round(nb.files)+1)

list.names.of.file <- NULL
list.nb.ctries <- list()
ctries.already.done <- 0
for(i in 1:nb.files){
  ctries.on.chart.i <- (ctries.already.done+1):(ctries.already.done+min(nb.ctries-ctries.already.done,max.nb.rows.per.figure))
  ctries.already.done <- ctries.already.done + length(ctries.on.chart.i)
  list.nb.ctries[[i]] <- ctries.on.chart.i
  ctries.string <- list.of.ctries[ctries.on.chart.i[1]]
  for(j in 2:length(ctries.on.chart.i)){
    ctries.string <- paste(ctries.string,list.of.ctries[ctries.on.chart.i[j]],  sep="-")
  }
  list.names.of.file <- c(list.names.of.file,paste("FigureCDS", ctries.string,  sep="_"))
}


for(j in 1:length(list.names.of.file)){
  
  list_of_countries <- list.nb.ctries[[j]]
  nb.ctries.on.chart.j <- length(list_of_countries)
  
  FILE = paste("figures/",list.names.of.file[j],suffix,
               ".pdf",sep="")
  
  pdf(file=FILE,pointsize=11,width=9, height=2*nb.ctries.on.chart.j)
  
  ColS <- ifelse(indic.colors,rgb(0,0,1),rgb(0.5,0.5,0.5))
  
  par(mfrow=c(nb.ctries.on.chart.j,nb.of.CDS.maturities))
  
  par(plt=c(.15,.95,.25,.8))
  
  counter.ctry <- 0
  
  for(ctry in list.of.ctries[list_of_countries]){
    counter.ctry <- counter.ctry + 1
    
    path.results <- paste("results/res_estim_",ctry,suffix,".Rdat",sep="")
    load(file=path.results)
    
    vec.dates <- DATA$DATE
    
    # Construct estimated model:
    model.est     <- Theta2Model(full_theta_est,model_ini_sol,bounds)
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

    aux <- compute.CDS(model.sol.est,all.X,
                       max.h = max(DATA$maturity.CDS),
                       indic.cpp = 1)
    modeled.CDS.Q <- aux$CDS.spds
    
    # Under P:
    model_new_P <- model.sol.est
    model_new_P$gamma <- 0.00001
    model_new_sol_P <- solve_model(model_new_P,indic_compute_beta = 1)
    aux <- compute.CDS(model_new_sol_P,all.X,
                       max.h = max(DATA$maturity.CDS),
                       indic.cpp = 1)
    modeled.CDS.P <- aux$CDS.spds
    
    for(i.maturity in 1:nb.of.CDS.maturities){
      modeled.cds.q <- pmax(40000*modeled.CDS.Q[,DATA$maturity.CDS[i.maturity]],0)
      modeled.cds.p <- pmax(40000*modeled.CDS.P[,DATA$maturity.CDS[i.maturity]],0)
      observd.cds.q <- pmax(40000*DATA$CDS[,i.maturity],0)
      
      min.y <- 0
      max.y <- 1.3*max(modeled.cds.q,
                       #modeled.cds.p,
                       observd.cds.q,na.rm = TRUE)
      
      par(mfg=c(counter.ctry,i.maturity))

      plot(vec.dates,modeled.cds.q,#type="l",
           las=1,lwd=1,pch=3,
           ylim=c(min.y,max.y),
           xlab="",ylab="",col="white",
           main = paste(DATA$name.of.country$short.name,
                        ", Maturity: ",toString(DATA$maturity.CDS[i.maturity]/4),
                        ifelse(DATA$maturity.CDS[i.maturity]<=4," year"," years"),sep=""),
           cex.main=parcex, cex.lab=parcex, cex.axis=parcex)
      lines(vec.dates,modeled.cds.q,lty=1,col=ColQ,lwd=3)
      points(vec.dates,observd.cds.q,lty=1,lwd=1,pch=3,col="black")
      
      abline(h=0,col="grey",lty=1)
      
      if(indic.plot.P){
        lines(vec.dates,modeled.cds.p,lty=1,col="black",lwd=1)
      }
      
      if((counter.ctry==1)*(i.maturity==1)){
        legend("topright", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend
               c("Obs. CDS",
                 "Model"),
               lty=c(NaN,1,3), # gives the legend appropriate symbols (lines)
               lwd=c(1,3,3), # line width
               col=c("black",ColQ),
               bg="white",
               pch=c(3,NaN,NaN),
               seg.len = 1,
               ncol=1,
               cex=parcexleg)
      }
      
    }
  }
  
  dev.off()
  
}





