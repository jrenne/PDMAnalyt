
round.fixed.length <- function(x,nb.dec){
  format.nb <- paste("%.",nb.dec,"f",sep="")
  output <- paste("",sprintf(format.nb,x),"",sep="")
  return(output)
}


sim.mod.ext <- function(mod.solution ,nb.sim, x.0){
  
  n.w <- dim(mod.solution$Phi)[1]
  n.x <- dim(mod.solution$Phi_x)[1]
  n_eta <- length(mod.solution$sigma_c)
  
  eta.sim <- matrix(rnorm(nb.sim * n_eta), n_eta, nb.sim)
  epsilon.sim <- matrix(rnorm(nb.sim * n.w), n.w, nb.sim)
  epsilon.s.sim <- rnorm(nb.sim)
  
  epsilon.x <- rbind(epsilon.sim,
                     epsilon.s.sim,
                     eta.sim)
  
  x.sim <- matrix(0,n.x,nb.sim)
  x.sim[,1] <- x.0
  for(t in 2:nb.sim){
    x.sim[,t] <- mod.solution$Mu_x + mod.solution$Phi_x %*% x.sim[,(t-1)] + mod.solution$Sigma_x %*% epsilon.x[,t]   
  }
  
  # Consumption
  delta_c.sim.lowfreq <- c(mod.solution$mu_c + t(mod.solution$Lambda_c)%*%x.sim[1:n.w,])
  delta_c.sim         <- c(delta_c.sim.lowfreq + t(mod.solution$sigma_c)%*%eta.sim)
  #Inflation
  pi.sim <- c(mod.solution$mu_pi + t(mod.solution$Lambda_pi)%*%x.sim[1:n.w,] + t(mod.solution$sigma_pi)%*%eta.sim)
  # Output
  delta_y.sim.lowfreq <- c(mod.solution$mu_y + t(mod.solution$Lambda_y)%*%x.sim[1:n.w,])
  delta_y.sim         <- c(delta_y.sim.lowfreq + t(mod.solution$sigma_y)%*%eta.sim)
  # mu_star
  delta_s_star.sim <- c(mod.solution$mu_s.star + t(mod.solution$Lambda_s)%*%x.sim[1:n.w,])
  
  return(list(
    x.sim       = x.sim,
    delta_c.sim = delta_c.sim,
    delta_y.sim = delta_y.sim,
    delta_c.sim.lowfreq = delta_c.sim.lowfreq,
    delta_y.sim.lowfreq = delta_y.sim.lowfreq,
    delta_s_star.sim = delta_s_star.sim,
    pi.sim      = pi.sim
  ))
}




prepare.vector.of.param.4.table <- function(model.sol){
  
  param.notation <- c(
    "$RR$ (eq.\\,\\ref{eq:_CDS})",
    "$\\chi$ (Sub.\\,\\ref{sub:debtaccum})",
    "$\\gamma$ (eq.\\,\\ref{eq:sdf})",
    "$b_c = b_y$ (eq.\\,\\ref{eq:cypis_affine})",
    "$-b_\\pi$ (eq.\\,\\ref{eq:cypis_affine})",
    "$\\exp(\\bar{d})$ (eq.\\,\\ref{eq:sd})",
    "$\\alpha$ (eq.\\,\\ref{eq:probadefF})",
    "$\\sigma_\\nu$ (eq.\\,\\ref{eq:probadefF})",
    "$h^*$ (eq.\\,\\ref{eq:qhstar})",
    "$\\mu_{s^*}$ (eq.\\,\\ref{eq:ell})",
    paste("\\hline $sd_t$ (eq.\\,\\ref{eq:sd})  \\\\ \\hline ",
          "$\\gamma_d$",sep=""),
    "$\\overline{sd}$ (eq.\\,\\ref{eq:notations_rd_sd})")
  
  param.values <- c(model.sol$RR,
                    model.sol$Chi,
                    model.sol$gamma,
                    model.sol$b_c,
                    -model.sol$b_pi,
                    exp(model.sol$d_bar)/4,
                    model.sol$alpha,
                    model.sol$sigma.nu,
                    model.sol$h_star/4,
                    model.sol$mu_s.star,
                    model.sol$gamma_d,
                    model.sol$sd_bar)
  
  n_w <- dim(model.sol$Phi)[1]
  
  count.variable <- 0
  for(indic.variable in c("s","c","y","pi")){
    count.variable <- count.variable + 1
    
    if(indic.variable != "s"){
      param.notation <- c(param.notation,
                          paste(
                            "\\hline $",
                            ifelse(indic.variable=="pi","\\pi",indic.variable),
                            "_t$ (eq.\\,\\ref{eq:cypis_affine})   \\\\ \\hline ",
                            "$\\mu_",
                            ifelse(indic.variable=="pi","\\pi",indic.variable),
                            "$",sep=""))
      eval(parse(text = gsub(" ","",paste("param.aux <- model.sol$mu_",
                                          indic.variable,
                                          sep=""))))
      param.values <- c(param.values,
                        param.aux)
    }
    
    # Show only relevant value of Lambda_y:
    if(indic.variable!="pi"){
      if(indic.variable %in% c("y","c")){
        i <- 1
      }else( # budget surplus
        i <- 3
      )
      param.notation <- c(param.notation,
                          paste("$\\Lambda_{",ifelse(indic.variable=="pi",
                                                     "\\pi",indic.variable),",",
                                toString(i),"}$",sep=""))
      eval(parse(text = gsub(" ","",paste("param.aux <- model.sol$Lambda_",
                                          indic.variable,
                                          "[",toString(i),"]",
                                          sep=""))))
      param.values <- c(param.values,
                        param.aux)
    }
    
    param.notation <- c(param.notation,
                        paste("$\\lVert\\sigma_",ifelse(indic.variable=="pi","\\pi",indic.variable),"\\rVert$",sep=""))
    eval(parse(text = gsub(" ","",paste("param.aux <- sqrt(sum(model.sol$sigma_",
                                        indic.variable,"^2))",
                                        sep=""))))
    param.values <- c(param.values,
                      param.aux)
  }
  
  indices.2.report <- rbind(c(1,1),
                            c(1,2),
                            c(2,2),
                            c(3,3),
                            c(3,4),
                            c(4,4))
  
  for(i in 1:dim(indices.2.report)[1]){
    if(i == 1){
      param.notation <- c(param.notation,
                          paste("\\hline $w_t$ (eqs.\\,\\ref{eq:w} and \\ref{eq:Phi})  \\\\ \\hline
                                  $\\Phi_{",
                                toString(indices.2.report[i,1]),",",
                                toString(indices.2.report[i,2]),"}$",sep=""))
    }else{
      param.notation <- c(param.notation,
                          paste("$\\Phi_{",
                                toString(indices.2.report[i,1]),",",
                                toString(indices.2.report[i,2]),"}$",sep="") )
    }
    eval(parse(text = gsub(" ","",paste("param.aux <- model.sol$Phi",
                                        "[",toString(indices.2.report[i,1]),",",
                                        toString(indices.2.report[i,2]),"]",
                                        sep=""))))
    param.values <- c(param.values,
                      param.aux)
  }
  
  return(list(
    param.notation = param.notation,
    param.values   = param.values
  ))
}
