# =================================================================
# Produces figure with estimated factors w_t
# =================================================================


nb.stdv  <- 2 # for Conf. Intervals
n_w2plot <- 4 # the first n_w2plot w factors willbe plotted


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
  list.names.of.file <- c(list.names.of.file,paste("FigureFactors", ctries.string,  sep="_"))
}

for(j in 1:length(list.names.of.file)){
  
  list_of_countries <- list.nb.ctries[[j]]
  nb.ctries.on.chart.j <- length(list_of_countries)
  

  FILE = paste("figures/",list.names.of.file[j],suffix,".pdf",sep="")
  
  pdf(file=FILE,pointsize=11,width=9, height=2*nb.ctries.on.chart.j)
  
  ColS <- ifelse(indic.colors,rgb(0,0,1),rgb(0.5,0.5,0.5))
  
  par(mfrow=c(nb.ctries.on.chart.j,n_w2plot))
  
  par(plt=c(.15,.95,.25,.77))
  
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
    
    X <- t(resEKF_est$W.smoothed)
    #X <- t(resEKF_est$W.updated)
    
    for(i.factor in 1:n_w2plot){
      
      Y <- X[,i.factor]
      stdv.Y <- sqrt(
          sapply(resEKF_est$P.smoothed,function(x){x[i.factor,i.factor]})
          )
      
      min.y <- min(Y - nb.stdv*stdv.Y)
      max.y <- max(Y + nb.stdv*stdv.Y)

      par(mfg=c(counter.ctry,i.factor))
      
      Title <- expression(paste(DATA$name.of.country$short.name,", ","sigma[c]",sep=""))
      
      if(i.factor == 1){
        eval(parse(text = gsub(" "," ",paste("Title <- expression(paste('",DATA$name.of.country$long.name,", ', ",
                                             "sigma[y]","))",sep=""))))
      }else if(i.factor==2){
        eval(parse(text = gsub(" "," ",paste("Title <- expression(paste('",DATA$name.of.country$long.name,", ', ",
                                             "sigma[pi]","))",sep=""))))
      }else if(i.factor==3){
        eval(parse(text = gsub(" "," ",paste("Title <- expression(paste('",DATA$name.of.country$long.name,", ', ",
                                             "sigma[s]","))",sep=""))))
      }else if(i.factor==4){
        eval(parse(text = gsub(" "," ",paste("Title <- expression(paste('",DATA$name.of.country$long.name,", ', ",
                                             "sigma['*']","))",sep=""))))
      }
      
      plot(vec.dates,Y,type="l",
           las=1,lwd=3,
           ylim=c(min.y,max.y),
           xlab="",ylab="",col=ColS,  #rgb(.5,.5,.5)
           main = Title,
           cex.main=parcex, cex.lab=parcex, cex.axis=parcex)
      
      lwer.bd <- Y - nb.stdv*stdv.Y
      uper.bd <- Y + nb.stdv*stdv.Y
      
      color.area <- "#11111144" # color used for the first area (mild default risk)
      polygon(c(vec.dates,rev(vec.dates)),c(lwer.bd,rev(uper.bd)),
              col=color.area,border=NA)

      abline(h=0,col="grey",lty=1)
      
    }
  }
  
  dev.off()
  
}

