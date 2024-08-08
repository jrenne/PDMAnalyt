# ==============================================================================
# Produces figure with fiscal limits
# ==============================================================================


# Define groups of estimations:
group.estimations <- list()

group.estimations[[1]] <- list(
  charact.model.large.alpha,
  charact.model.small.alpha,
  charact.model.noCDS
)

group.estimations[[2]] <- list(
  charact.model.low_bc,
  charact.model.smallGamma,
  charact.model.largeGamma,
  charact.model.smallGamma.largeRR,
  charact.model.largeGamma.smallRR
)

group.estimations[[3]] <- list(
  charact.model.low_bc,
  charact.model.large.alpha,
  charact.model.small.alpha,
  charact.model.smallGamma,
  charact.model.largeGamma,
  charact.model.smallGamma.largeRR,
  charact.model.largeGamma.smallRR,
  charact.model.noCDS
)


ColS <- ifelse(indic.colors,rgb(0,0,1),rgb(0.5,0.5,0.5))


for(count.group in 1:length(group.estimations)){
  
  all.charact.model <- group.estimations[[count.group]]
  
  if(all.charact.model[[1]]$suffix!=""){
    new.all.charact.model <- list()
    new.all.charact.model[[1]] <- charact.model.baseline
    for(i in 1:length(all.charact.model)){
      new.all.charact.model[[i+1]] <- all.charact.model[[i]]
    }
  }
  all.charact.model <- new.all.charact.model
  
  
  # Latex table:
  Tab.robust.beta   <- matrix(NaN,length(all.charact.model),length(list.of.ctries))
  Tab.robust.loglik <- matrix(NaN,length(all.charact.model),length(list.of.ctries))
  
  # First step, construct names of files
  
  # compute number of files:
  nb.ctries <- length(list.of.ctries)
  
  
  FILE = paste("figures/FigureFL_robustness","_group",count.group,".pdf",sep="")
  
  pdf(file=FILE,pointsize=11,width=9, height=7)
  
  par(mfrow=c(2,3)) 
  
  par(plt=c(.15,.95,.2,.85))
  
  counter.ctry <- 0
  
  for(ctry in list.of.ctries){
    counter.ctry <- counter.ctry + 1
    
    print(paste("======== ",ctry," ========",sep=""))
    
    for(case in 1:length(all.charact.model)){
      
      suffix <- all.charact.model[[case]]$suffix
      
      path.results <- paste("results/res_estim_",ctry,suffix,".Rdat",sep="")
      load(file=path.results)
      
      vec.dates <- DATA$DATE
      
      # Construct estimated model:
      model.est <- Theta2Model(full_theta_est,model_ini_sol,bounds)
      model.sol.est <- solve_model(model.est,indic_compute_beta = 1)
      
      n_w <- dim(model.sol.est$Phi)[1]
      n_x <- dim(model.sol.est$Phi_x)[1]
      
      resEKF_est <- EKF_smoother(model.sol.est,list.stdv,DATA)
      print(paste(suffix,": loglik = ",toString(round(resEKF_est$loglik,2)),
                  ", beta: ",round(model.sol.est$beta,4),
                  sep=""))
      
      # Store results:
      Tab.robust.beta[case,counter.ctry]   <- model.sol.est$beta
      Tab.robust.loglik[case,counter.ctry] <- resEKF_est$loglik
      
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
      
      min.y <- 0
      if(ctry=="US"){max.y <- 3}    
      if(ctry=="UK"){max.y <- 2}    
      if(ctry=="CA"){max.y <- 3}    
      if(ctry=="JP"){max.y <- 4}    
      if(ctry=="BR"){max.y <- 1.5}    
      if(ctry=="RU"){max.y <- .4}    
      if(ctry=="IN"){max.y <- 3}    
      if(ctry=="CN"){
        if(sum(unlist(lapply(all.charact.model,function(x){x$suffix}))=="_largeAlpha")>0){
          max.y <- 4
        }else{
          max.y <- 2
        }
      }    
      
      if(case==1){
        plot(vec.dates,elld,type="l",
             las=1,
             ylim=c(min.y,max.y),
             xlab="",ylab="",
             col = "grey",
             lwd = 2,
             lty = 1,
             cex.main=parcex, cex.lab=parcex, cex.axis=parcex)
        title(main=paste(DATA$name.of.country$long.name,sep=""),
              cex.main=parcex)
      }
      
      abline(h=0,lty=3)
      
      lines(vec.dates,ell1,
            col = all.charact.model[[case]]$col,
            lwd = all.charact.model[[case]]$lwd,
            lty = all.charact.model[[case]]$lty)
      
    }
    
  }
  
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("topleft", 
         c("Debt-to-GDP ratio","",
           paste("FL - ",unlist(lapply(all.charact.model,function(x){x$name.4.chart})),sep="")),
         lty=c(1,1,unlist(lapply(all.charact.model,function(x){x$lty}))), # gives the legend appropriate symbols (lines)       
         lwd=c(2,1,unlist(lapply(all.charact.model,function(x){x$lwd}))), # line width
         col=c("grey","white",unlist(lapply(all.charact.model,function(x){x$col}))),
         bg="white",
         seg.len = 3,
         cex=1.2)
  
  
  dev.off()
  
  
}





# Prepare Latex tables:

first.line <- ""
for(j in 1:length(list.of.ctries)){
  first.line <- paste(first.line,"&",list.of.ctries[j])
}
first.line <- paste(first.line,"\\\\")


Latex.tab.robust.beta   <- rbind(first.line,"\\hline")
Latex.tab.robust.loglik <- rbind(first.line,"\\hline")

for(i in 1:dim(Tab.robust.loglik)[1]){
  
  this.line.tab.robust.beta   <- all.charact.model[[i]]$name
  this.line.tab.robust.loglik <- all.charact.model[[i]]$name
  
  for(j in 1:dim(Tab.robust.loglik)[2]){
    this.line.tab.robust.beta <- paste(this.line.tab.robust.beta,"&",
                                       ifelse(Tab.robust.beta[i,j]<.01,"{\\bf \\color{red}",""),
                                       round.fixed.length(Tab.robust.beta[i,j],4),
                                       ifelse(Tab.robust.beta[i,j]<.01,"}",""),
                                       sep="")
    this.line.tab.robust.loglik <- paste(this.line.tab.robust.loglik,"&",
                                         round.fixed.length(Tab.robust.loglik[i,j],0),sep="")
  }
  
  this.line.tab.robust.beta   <- paste(this.line.tab.robust.beta,"\\\\",sep="")
  this.line.tab.robust.loglik <- paste(this.line.tab.robust.loglik,"\\\\",sep="")
  
  Latex.tab.robust.beta   <- rbind(Latex.tab.robust.beta,
                                   this.line.tab.robust.beta)
  Latex.tab.robust.loglik <- rbind(Latex.tab.robust.loglik,
                                   this.line.tab.robust.loglik)
}


Latex.tab.robust.beta <- rbind(
  "\\begin{table}[h!]",
  "\\caption{Estimates of $\\beta$}",
  "\\label{tab:robustBeta}",
  "\\begin{center}",
  "\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}ccccccccccccccccccc}",
  "\\\\",
  "\\hline",
  Latex.tab.robust.beta,
  "\\hline",
  "\\end{tabular*}",
  "\\end{center}",
  "\\begin{footnotesize}",
  "\\parbox{\\linewidth}{",
  "Note: This table reports the estimates of parameter $\\beta$ obtained while imposing different types of restrictions during the estimation.",
  "{\\bf Low b\\_y}: $b_y$ is set to 10\\% (versus 20\\% in the baseline case);",
  "{\\bf Large alpha}: $\\alpha$ is set to 10 (in the baseline case, it is estimated, but smaller than 2);",
  "{\\bf Small alpha}: $\\alpha$ is set to 0.01;",
  "{\\bf Small gamma}: $\\gamma$ is set to 3 (versus 4 in the baseline case);",
  "{\\bf Large gamma}: $\\gamma$ is set to 5 (versus 4 in the baseline case);",
  "{\\bf Small gamma, large RR}: $\\gamma$ is set to 3 (versus 4 in the baseline case) and $b_y$ is adjusted to give the same $RR$ as in the baseline case;",
  "{\\bf Large gamma, small RR}: $\\gamma$ is set to 5 (versus 3 in the baseline case) and $b_y$ is adjusted to give the same $RR$ as in the baseline case;",
  "{\\bf No CDS data}: no CDS data are used in the estimation approach (i.e., there is no measurement equations involving CDS spreads).",
  "}",
  "\\end{footnotesize}",
  "\\end{table}"
)
write(Latex.tab.robust.beta,"tables/Tab.robust.beta.tex")

Latex.tab.robust.loglik <- rbind(
  "\\begin{table}[h!]",
  "\\caption{Maximum log-likelihoods}",
  "\\label{tab:robustLoglik}",
  "\\begin{center}",
  "\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}ccccccccccccccccccc}",
  "\\\\",
  "\\hline",
  Latex.tab.robust.loglik,
  "\\hline",
  "\\end{tabular*}",
  "\\end{center}",
  "\\begin{footnotesize}",
  "\\parbox{\\linewidth}{",
  "Note: This table reports the estimates of parameter $\\beta$ obtained while imposing different types of restrictions during the estimation.",
  "{\\bf Low b\\_y}: $b_y$ is set to 10\\% (versus 20\\% in the baseline case);",
  "{\\bf Large alpha}: $\\alpha$ is set to 10 (in the baseline case, it is estimated, but smaller than 2);",
  "{\\bf Small alpha}: $\\alpha$ is set to 0.01;",
  "{\\bf Small gamma}: $\\gamma$ is set to 3 (versus 4 in the baseline case);",
  "{\\bf Large gamma}: $\\gamma$ is set to 5 (versus 4 in the baseline case);",
  "{\\bf Small gamma, large RR}: $\\gamma$ is set to 3 (versus 4 in the baseline case) and $b_y$ is adjusted to give the same $RR$ as in the baseline case;",
  "{\\bf Large gamma, small RR}: $\\gamma$ is set to 5 (versus 3 in the baseline case) and $b_y$ is adjusted to give the same $RR$ as in the baseline case;",
  "{\\bf No CDS data}: no CDS data are used in the estimation approach (i.e., there is no measurement equations involving CDS spreads).",
  "}",
  "\\end{footnotesize}",
  "\\end{table}"
)
write(Latex.tab.robust.loglik,"tables/Tab.robust.loglik.tex")
