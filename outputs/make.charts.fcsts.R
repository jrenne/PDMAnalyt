# =================================================================
# Produces figure with forecasts
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
  list.names.of.file <- c(list.names.of.file,paste("FigureFcsts",
                                                   variable, ctries.string,  sep="_"))
}


for(j in 1:length(list.names.of.file)){
  
  list_of_countries <- list.nb.ctries[[j]]
  nb.ctries.on.chart.j <- length(list_of_countries)
  
  FILE = paste("figures/",list.names.of.file[j],suffix,
               ".pdf",sep="")
  
  pdf(file=FILE,pointsize=11,width=9, height=2*nb.ctries.on.chart.j)
  
  ColS <- ifelse(indic.colors,rgb(0,0,1),rgb(0.5,0.5,0.5))
  
  par(mfrow=c(nb.ctries.on.chart.j,max.nb.of.fcsts))
  
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
    
    nb_obs_macro <- 4
    nb_mat_CDS   <- length(DATA$maturity.CDS)
    nb_mat_yds   <- length(DATA$maturity.yields)
    nb_horiz_fcsts <- length(DATA$horizon.fcsts)
    aux <- nb_obs_macro + nb_mat_yds + nb_mat_CDS
    debt.fcsts <- t(resEKF_est$Obs.predicted[(aux+1):(aux+nb_horiz_fcsts),])
    aux <- aux + nb_horiz_fcsts
    inflation.fcsts <- t(resEKF_est$Obs.predicted[(aux+1):(aux+nb_horiz_fcsts),])
    aux <- aux + nb_horiz_fcsts
    growth.fcsts <- t(resEKF_est$Obs.predicted[(aux+1):(aux+nb_horiz_fcsts),])
    
    for(i.horizon in 1:nb_horiz_fcsts){
      
      if(variable=="debt"){
        Y      <- debt.fcsts[,i.horizon]/model.sol.est$freq
        Y.data <- DATA$debt.fcsts[,i.horizon]/model.sol.est$freq
      }
      if(variable=="infl"){
        Y      <- inflation.fcsts[,i.horizon]
        Y.data <- DATA$infl.fcsts[,i.horizon]
      }
      if(variable=="growth"){
        Y      <- growth.fcsts[,i.horizon]
        Y.data <- DATA$growth.fcsts[,i.horizon]
      }
      
      min.y <- ifelse(min(Y,Y.data,na.rm = TRUE)<0,
                      1.4*min(Y,Y.data,na.rm = TRUE),0.5*min(Y,Y.data,na.rm = TRUE))
      max.y <- 1.5*max(Y,Y.data,na.rm = TRUE)
      
      par(mfg=c(counter.ctry,i.horizon))
      
      plot(vec.dates,Y,type="l",
           las=1,lwd=3,
           ylim=c(min.y,max.y),
           xlab="",ylab="",col=ColS,  #rgb(.5,.5,.5)
           main = paste(DATA$name.of.country$short.name,
                        ", Horizon: ",toString(DATA$horizon.fcsts[i.horizon]/4),
                        ifelse(DATA$horizon.fcsts[i.horizon]<=4," year"," years"),sep=""),
           cex.main=parcex, cex.lab=parcex, cex.axis=parcex)
      
      points(vec.dates,Y.data,col="black",lwd=1,pch=3) 
      
      abline(h=0,col="grey",lty=1)
      
    }
  }
  
  dev.off()
  
}





