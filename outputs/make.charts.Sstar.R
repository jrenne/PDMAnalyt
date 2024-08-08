# =================================================================
# Produces figure with maximum budget surplus
# =================================================================


ColS <- ifelse(indic.colors,rgb(0,0,1),rgb(0.5,0.5,0.5))

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
  list.names.of.file <- c(list.names.of.file,paste("FigureSstar", ctries.string,  sep="_"))
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
    
    # beta.d.d.star <- model.sol.est$beta * (DATA$d - model.sol.est$d_star)
    beta.d.d.star <- model.sol.est$beta * (all.X[,n_w+3] - model.sol.est$d_star)

    s_star <- model.sol.est$mu_star + t(resEKF_est$W.updated) %*% model.sol.est$sigma_star
    
    check_s <- beta.d.d.star + t(resEKF_est$W.updated) %*% model.sol.est$sigma_s
    
    stdv_s_star <- lapply(resEKF_est$P.updated,
                         function(x){
                           sqrt(t(model.sol.est$sigma_star) %*% x %*% 
                                  model.sol.est$sigma_star)})
    stdv_s_star <- as.numeric(stdv_s_star)

    min.y <- 1.2*min(DATA$s,s_star-2*stdv_s_star,beta.d.d.star)
    max.y <- 1.2*max(DATA$s,s_star+2*stdv_s_star,beta.d.d.star)
    
    plot(vec.dates,s_star,type="l",
         las=1,lwd=2,
         ylim=c(min.y,max.y),
         xlab="",ylab="",col="black",
         cex.main=parcex, cex.lab=parcex, cex.axis=parcex)
    title(main=paste(DATA$name.of.country$long.name,sep=""),
          cex.main=parcex)
    
    lwer.bd <- s_star - 2*stdv_s_star
    uper.bd <- s_star + 2*stdv_s_star
    
    color.area <- "#11111144" # color used for the first area (mild default risk)
      polygon(c(vec.dates,rev(vec.dates)),c(lwer.bd,rev(uper.bd)),
              col=color.area,border=NA)
      
      
    #main = paste(model.sol$name.of.country$long.name,sep=""))
    lines(vec.dates,s_star,lty=1,col="black",lwd=2)
    lines(vec.dates,DATA$s,lty=1,col="blue",lwd=2)
    lines(vec.dates,beta.d.d.star,lty=3,col="blue",lwd=1)
    abline(h=0,col="grey")
    
    
    if(counter.ctry==1){
      legend("bottomleft", 
             c("s*",
               "budget surplus",
               expression(paste(beta,"(d-d*)",sep=""))),
             lty=c(1,1,3), # gives the legend appropriate symbols (lines)       
             lwd=c(2,2,2), # line width
             col=c("black","blue","blue"),
             bg="white",
             seg.len = 2,
             cex=parcexleg)
    }
  }
  
  dev.off()

}





