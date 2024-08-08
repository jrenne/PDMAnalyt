# =================================================================
# Produces figure with yields
# =================================================================

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
  list.names.of.file <- c(list.names.of.file,paste("FigureYds", ctries.string,  sep="_"))
}


for(j in 1:length(list.names.of.file)){
  
  list_of_countries <- list.nb.ctries[[j]]
  nb.ctries.on.chart.j <- length(list_of_countries)
  
  FILE = paste("figures/",list.names.of.file[j],suffix,
               ".pdf",sep="")
  
  pdf(file=FILE,pointsize=11,width=9, height=2*nb.ctries.on.chart.j)

  ColS <- ifelse(indic.colors,rgb(0,0,1),rgb(0.5,0.5,0.5))
  
  par(mfrow=c(nb.ctries.on.chart.j,max.nb.of.yds))
  
  par(plt=c(.15,.95,.25,.77))
  
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
    
    yds <- all.X[,1:n_w] %*% model.sol.est$a.nom + 
      matrix(1,T,1) %*% model.sol.est$b.nom
    
    for(i.horizon in 1:length(DATA$maturity.yields)){
      
      Y      <- 400*yds[,DATA$maturity.yields[i.horizon]]
      Y.data <- 400*DATA$yields[,i.horizon]
      
      min.y <- ifelse(min(Y,Y.data)<0,1.4*min(Y,Y.data),0.5*min(Y,Y.data))
      max.y <- max(1.4*max(Y,Y.data),6)
      
      par(mfg=c(counter.ctry,i.horizon))
      
      plot(vec.dates,Y,type="l",
           las=1,lwd=3,
           ylim=c(min.y,max.y),
           xlab="",ylab="",col=ColS,  #rgb(.5,.5,.5)
           main = paste(DATA$name.of.country$short.name,
                        ", Maturity: ",toString(DATA$maturity.yields[i.horizon]/4),
                        ifelse(DATA$maturity.yields[i.horizon]<=4," year"," years"),sep=""),
           cex.main=parcex, cex.lab=parcex, cex.axis=parcex)
      
      points(vec.dates,Y.data,col="black",lwd=1,pch=3) 
      
      abline(h=0,col="grey",lty=1)
      
      if((counter.ctry==1)*(i.horizon==1)){
        legend("topright", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend 
               c("Model-implied yields","Observed yields"),
               lty=c(1,NaN), # gives the legend appropriate symbols (lines)       
               lwd=c(3,1), # line width
               col=c(ColS,"black"), #rgb(.5,.5,.5)
               bg="white",
               pch=c(NaN,3),
               seg.len = 2,
               cex=parcexleg)
      }
    }
  }
  
  dev.off()
  
}





