

maxH <- 10
nb_iter <- 30 # to solve models

# Create stylized models

low_pi <- -.01
med_pi <- .03
hig_pi <- .07

low_y <- -.01
med_y <- .02
hig_y <- .05

mu_y  <- matrix(c(low_y,hig_y,med_y,hig_y,low_y),ncol=1)

# Demand:
mu_pi_demand  <- matrix(c(low_pi,hig_pi,med_pi,hig_pi,low_pi),ncol=1)
# Supply:
mu_pi_supply  <- matrix(c(hig_pi,low_pi,med_pi,low_pi,hig_pi),ncol=1)

rho1 <- .8
rho2 <- .6
Omega <- diag(rep(rho2,5))
Omega[1,2] <- 1 - rho2
Omega[2,3] <- 1 - rho2
Omega[3,1] <- (1 - rho1)/2
Omega[3,3] <- rho1
Omega[3,4] <- (1 - rho1)/2
Omega[4,5] <- 1 - rho2
Omega[5,3] <- 1 - rho2


Model <- list(  mu_pi = NaN,
                mu_y = mu_y,
                Omega = Omega,
                gamma = 10,
                delta = .99,
                chi = .8,
                mu_eta = 1*mu_y,
                beta = .1,
                alpha = .3,
                d_star =.8,
                s_star =.05,
                RR= .5,
                nu_pi = 0,
                nu_y = 0,
                sigma_nu = .1)

Model_Demand <- Model
Model_Demand$mu_pi <- mu_pi_demand

Model_Supply <- Model
Model_Supply$mu_pi <- mu_pi_supply




# ------------------------------------------------------------------------------
# Prepare Latex table showing specifications
# ------------------------------------------------------------------------------

make.entry <- function(x,format.nb,dollar=1){
  if(dollar==1){
    output <- paste("$",sprintf(format.nb,x),"$",sep="")
  }else{
    output <- sprintf(format.nb,x)
  }
  return(output)
}
format.nb0 <- paste("%.",0,"f",sep="")
format.nb1 <- paste("%.",1,"f",sep="")
format.nb2 <- paste("%.",2,"f",sep="")
format.nb3 <- paste("%.",3,"f",sep="")

nb_m <- dim(Omega)[1]

columns <- NULL
for(i in 1:6){
  columns <- paste(columns,"r",sep="")
}
for(i in 1:nb_m){
  columns <- paste(columns,"c",sep="")
}

latex_table <- rbind("\\begin{table}[ph!]",
                     "\\caption{Stylized models: parameterizations}",
                     "\\label{tab:DemSup_parameters}",
                     paste("\\begin{tabular*}{\\textwidth}{c@{\\extracolsep{\\fill}}",columns,"}",sep=""),
                     "\\hline",
                     paste("Regime&\\multicolumn{2}{c}{$\\mu_\\pi$}&&$\\mu_y$&&\\multicolumn{",nb_m,"}{c}{$\\Omega$}\\\\",sep=""),
                     "&D&S\\\\",
                     "\\hline")

for(i in 1:nb_m){
  pi_i <- NULL
  for(j in 1:nb_m){
    pi_i <- paste(pi_i,"&",make.entry(Model$Omega[i,j],format.nb3),sep="")
  }
  this_line <- paste(i,
                     "&",make.entry(Model_Demand$mu_pi[i],format.nb3),
                     "&",make.entry(Model_Supply$mu_pi[i],format.nb3),
                     "&&",make.entry(Model$mu_y[i],format.nb3),
                     "&",pi_i,
                     "\\\\",sep="")
  latex_table <- rbind(latex_table,
                       this_line)
}

latex_table <- rbind(latex_table,
                     "\\hline",
                     "\\end{tabular*}",
                     "\\begin{footnotesize}",
                     paste("\\parbox{\\linewidth}{\\textit{Notes}: This table shows the parameterizations of the stylized demand/supply models. ",
                           "We also have $\\alpha=",make.entry(Model$alpha,format.nb = format.nb2,dollar=0),"$, ",
                           "$\\beta=",make.entry(Model$beta,format.nb = format.nb2,dollar=0),"$, ",
                           "$\\gamma=",make.entry(Model$gamma,format.nb = format.nb2,dollar=0),"$, ",
                           "$\\beta=",make.entry(Model$beta,format.nb = format.nb2,dollar=0),"$, ",
                           "$d^*=",make.entry(Model$d_star,format.nb = format.nb2,dollar=0),"$, ",
                           "$s^*=",make.entry(Model$s_star,format.nb = format.nb2,dollar=0),"$, ",
                           "$RR=",make.entry(Model$RR,format.nb = format.nb2,dollar=0),"$.}",
                           sep=""),
                     "\\end{footnotesize}",
                     "\\end{table}")

name.of.file <- "table_param_DemSup"
latex.file <- paste(name.of.file,".txt", sep="")
write(latex_table, paste("tables/",latex.file,sep=""))






# ------------------------------------------------------------------------------
# Prepare plot showing average yield curves in the two models
# ------------------------------------------------------------------------------


FILE = paste("figures/Figure_expected_returns_DemaSupp.pdf",sep="")
pdf(file=FILE, pointsize=10, width=9, height=5)

par(mfrow=c(1,2))
par(plt=c(.12,.95,.1,.85))

for(regime in c("Demand","Supply")){
  
  eval(parse(text = gsub(" ","",paste("Model <- Model_",regime,sep=""))))
  
  RES <- prepare_returns_yds(Model,maxH)
  
  plot(1:maxH,rep(0,maxH),las=1,
       xlab="",ylab="",type="l",ylim=c(0.0,.11),lwd=2,col="white",
       main=paste(ifelse(regime=="Demand","(a) ","(b) "),regime,"-driven economy",sep=""))
  
  lower_bound <- RES$avg_annual_nominal_returns - 1.0 * RES$Std_annual_nominal_returns
  upper_bound <- RES$avg_annual_nominal_returns + 1.0 * RES$Std_annual_nominal_returns
  polygon(c(1:maxH,rev(1:maxH)),c(lower_bound,rev(upper_bound)),border = NaN,
          col="#66AAAA44")
  lines(1:maxH,c(RES$avg_annual_nominal_returns),col="blue",lwd=2)
  
  lines(1:maxH,RES$avg_annual_TIPS_returns,col="red",lwd=2)
  
  lower_bound <- RES$avg_annual_GDPLB_returns - 1.0 * RES$Std_annual_GDPLB_returns
  upper_bound <- RES$avg_annual_GDPLB_returns + 1.0 * RES$Std_annual_GDPLB_returns
  polygon(c(1:maxH,rev(1:maxH)),c(lower_bound,rev(upper_bound)),border = NaN,
          col="#AAAAAA44")
  lines(1:maxH,RES$avg_annual_GDPLB_returns,col="black",lwd=2,lty=2)
  
  grid()
  
  if(regime=="Supply"){
    legend("topright", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend
           c("Nominal bonds","TIPS","GDP-LB"),
           lty=c(1,1,2), # gives the legend appropriate symbols (lines)
           lwd=c(2,2,2), # line width
           col=c("blue","red","black"),
           bg="white",
           #pch=c(3,NaN,NaN),
           seg.len = 3,
           cex=1,title = "Annualized expected nominal returns")
    
  }
}

dev.off()



stop()


# Debt issuance strategies: ----------------------------------------------------
# Prepare table showing performances of debt issuance strategies
#    for different values of chi
# ------------------------------------------------------------------------------

values_of_chi <- c(.2,.9)

grids <- make_grid(nb_grid = 27,
                   max_d = 1.5,
                   max_rr=.12,
                   sigma_eps = .03,
                   all_quantiles_eps = c(-2,-1,1,2))

extensions       <- c("_nominal","_TIPS","_GDPLB")
names_extensions <- c("Nominal","TIPS","GDP-LBs")

outputs <- c("mean_d","stdv_d","mean_rr","stdv_rr",
             "stdv_Delta_d","avg_PD[maxH]","avg_spreads[maxH]")
latex.column.names <- c("$\\mathbb{E}(d)$","$\\sqrt{\\mathbb{V}ar(d)}$",
                        "$\\mathbb{E}(r)$","$\\sqrt{\\mathbb{V}ar(r)}$",
                        "$\\sqrt{\\mathbb{V}ar(\\Delta d)}$",
                        "$\\mathbb{E}(PD)$","$\\mathbb{E}(spd)$")


columns <- "c" # empty column
for(i in 1:length(latex.column.names)){
  columns <- paste(columns,"c",sep="")
}
column_names <- NULL
for(i in 1:length(latex.column.names)){
  column_names <- paste(column_names,"&",latex.column.names[i],sep="")
}

latex_table <- rbind("\\begin{table}[ph!]",
                     "\\caption{Performances of debt issuance strategies --- stylized version of the model}",
                     "\\label{tab:DemSup}",
                     paste("\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}",columns,"}",sep=""),
                     "\\hline",
                     paste("&",column_names,"\\\\",sep=""),
                     "\\hline")

for(chi in values_of_chi){
  
  latex_table <- rbind(latex_table,
                       "\\\\",
                       "\\hline",
                       paste("\\multicolumn{",length(outputs)+2,"}{l}{Coupon decay rate $\\chi=",chi,"$}\\\\"))
    
  Model_Demand$chi    <- chi
  Model_Supply$chi    <- chi
  
  for(regime in c("Demand","Supply")){
    
    latex_table <- rbind(latex_table,"\\\\",
                         paste("&&\\multicolumn{",length(outputs),"}{l}{",
                               regime,"-driven economy ($\\chi=",chi,"$)}\\\\",sep=""),
                         "\\hline")
    
    eval(parse(text = gsub(" ","",paste("Models <- prepare_and_solve_3(Model_",
                                        regime,",grids,nb_iter)",sep=""))))

    strat_nominal <- run_strategy(Models$Model_solved_nominal,maxH)
    strat_TIPS    <- run_strategy(Models$Model_solved_TIPS,maxH)
    strat_GDPLB   <- run_strategy(Models$Model_solved_GDPLB,maxH)
    
    M <- rbind(c(strat_nominal$mean_d,strat_TIPS$mean_d,strat_GDPLB$mean_d),
               c(strat_nominal$stdv_d,strat_TIPS$stdv_d,strat_GDPLB$stdv_d),
               c(strat_nominal$mean_rr,strat_TIPS$mean_rr,strat_GDPLB$mean_rr),
               c(strat_nominal$stdv_rr,strat_TIPS$stdv_rr,strat_GDPLB$stdv_rr),
               c(strat_nominal$stdv_Delta_d,strat_TIPS$stdv_Delta_d,strat_GDPLB$stdv_Delta_d),
               c(strat_nominal$avg_PD[maxH],strat_TIPS$avg_PD[maxH],strat_GDPLB$avg_PD[maxH]),
               c(strat_nominal$avg_spreads[maxH],strat_TIPS$avg_spreads[maxH],strat_GDPLB$avg_spreads[maxH]))
    colnames(M) <- names_extensions
    rownames(M) <- latex.column.names
    print(M)
    
    count_extension <- 0
    for(extension in extensions){
      count_extension <- count_extension + 1
      this.line <- paste(names_extensions[count_extension],"&")
      for(output in outputs){
        eval(parse(text = gsub(" ","",paste("this.line <- paste(this.line,'&',",
                                            "make.entry(strat",extension,"$",output,",format.nb2),sep='')",sep=""))))
      }
      this.line <- paste(this.line,"\\\\")
      latex_table <- rbind(latex_table,this.line)
    }
    
  }
}

latex_table <- rbind(latex_table,
                     "\\\\",
                     "\\hline",
                     "\\end{tabular*}",
                     "\\begin{footnotesize}",
                     "\\parbox{\\linewidth}{\\textit{Notes}: This table shows performance metrics associated with three different debt issuance strategies; each strategy consists in issuing a given type of perpetuities: a nominal perpetuity ($\\kappa_\\pi=0$ and $\\kappa_y=0$), an inflation-indexed perpetuity nominal ($\\kappa_\\pi=1$ and $\\kappa_y=0$), and a GDP-indexed perpetuity nominal ($\\kappa_\\pi=1$ and $\\kappa_y=1$). We consider two different values of $\\chi$ (the higher $\\chi$, the higher the average debt maturity). '$d$' denotes the debt-to-GDP ratio. '$r$' denotes the debt service, including debt indexation (in percent of GDP) '$PD$' stands for '10-year probability of default' (expressed in percent); '$spd$' stands for '10-year credit spread' (expressed in basis point).}",
                     "\\end{footnotesize}",
                     "\\end{table}")

name.of.file <- "table_DemSup"
latex.file <- paste(name.of.file,".txt", sep="")
write(latex_table, paste("tables/",latex.file,sep=""))



