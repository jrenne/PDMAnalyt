# ==============================================================================
# This script runs the demand/supply exercise
# ==============================================================================

# prepare additional figure showing effect of credit risk on nominal yields:
indic_addit_figure <- 0

# ==============================================================================
# ==============================================================================
elasticity_of_surpluses <- 1
#abs_nu_y <- .1
abs_nu_y <- 0
# ==============================================================================
# ==============================================================================

maxH <- 10
nb_iter <- 30 # to solve models
nb_iter_sdf <- 10 # to solve SDF
nb_grid <- 25

# Create stylized models

low_pi <- .0
med_pi <- .03
hig_pi <- .06

low_y <- .0
med_y <- .02
hig_y <- .04

# low_pi <- -.01
# med_pi <- .03
# hig_pi <- .07
# 
# low_y <- -.01
# med_y <- .02
# hig_y <- .05

# low_pi <- .03
# med_pi <- .03
# hig_pi <- .03
# 
# low_y <- .02
# med_y <- .02
# hig_y <- .02

mu_y  <- matrix(c(low_y,med_y,hig_y),ncol=1)

# Demand:
mu_pi_demand  <- matrix(c(low_pi,med_pi,hig_pi),ncol=1)
# Supply:
mu_pi_supply  <- matrix(c(hig_pi,med_pi,low_pi),ncol=1)

rho <- .8
Omega <- diag(rep(rho,3))
Omega[1,2] <- 1 - rho
Omega[3,2] <- 1 - rho
Omega[2,1] <- (1 - rho)/2
Omega[2,3] <- (1 - rho)/2

Model <- list(mu_pi = mu_pi_demand,
              mu_y = mu_y,
              Omega = Omega,
              gamma = 10,
              delta = .99,
              chi = .7,
              mu_eta = 0 * mu_y,
              beta = .1,
              alpha = .1,
              d_star = 1,
              s_star = NaN,
              RR = .5,
              nu_pi = 0,
              nu_y = -.1,
              sigma_eps = .02,
              sigma_nu = .1,
              kappa_pi = 0,
              kappa_y = 0)
# ==============================================================================
Model$nu_y <- - abs_nu_y
# ==============================================================================

res_aux <- compute_determ_steady_state(Model,
                                       indic_d_bar_from_s_star=0,
                                       d_bar = .8)
Model$s_star <- res_aux$s_star

Model$mu_eta <- elasticity_of_surpluses * Model$mu_y

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
                           "$\\gamma=",make.entry(Model$gamma,format.nb = format.nb0,dollar=0),"$, ",
                           "$\\beta=",make.entry(Model$beta,format.nb = format.nb2,dollar=0),"$, ",
                           "$d^*=",make.entry(Model$d_star,format.nb = format.nb2,dollar=0),"$, ",
                           "$s^*=",make.entry(Model$s_star,format.nb = format.nb2,dollar=0),"$, ",
                           "$\\sigma_\\nu=",make.entry(Model$sigma_nu,format.nb = format.nb2,dollar=0),"$, ",
                           # "$\\mu_\\eta=",elasticity_of_surpluses,"\\times \\mu_y$",
                           "$\\nu_y=",make.entry(Model$nu_y,format.nb = format.nb2,dollar=0),"$, ",
                           "$\\nu_\\pi=",make.entry(Model$nu_pi,format.nb = format.nb2,dollar=0),"$, ",
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
       xlab="",ylab="",type="l",ylim=c(0.02,.10),lwd=2,col="white",
       main=paste(ifelse(regime=="Demand","(a) ","(b) "),regime,"-driven economy",sep=""))
  
  lower_bound <- RES$avg_annual_nominal_returns - 1.0 * RES$Std_annual_nominal_returns
  upper_bound <- RES$avg_annual_nominal_returns + 1.0 * RES$Std_annual_nominal_returns
  polygon(c(1:maxH,rev(1:maxH)),c(lower_bound,rev(upper_bound)),border = NaN,
          col="#66AAAA44")
  lines(1:maxH,c(RES$avg_annual_nominal_returns),col="blue",lwd=2)
  
  lines(1:maxH,RES$avg_annual_ILB_returns,col="red",lwd=2)
  
  lower_bound <- RES$avg_annual_GDPLB_returns - 1.0 * RES$Std_annual_GDPLB_returns
  upper_bound <- RES$avg_annual_GDPLB_returns + 1.0 * RES$Std_annual_GDPLB_returns
  polygon(c(1:maxH,rev(1:maxH)),c(lower_bound,rev(upper_bound)),border = NaN,
          col="#AAAAAA44")
  lines(1:maxH,RES$avg_annual_GDPLB_returns,col="dark grey",lwd=2,lty=1)
  
  # Add term structure of nominal yields with credit risk:
  Model$nu_y <- 0 # to be sure
  Model$kappa_pi <- 0 # issuance of nominal perpet.
  Model$kappa_y  <- 0 # issuance of nominal perpet.
  Model_solved <- solve_ToyModel(Model,
                                 grids,nb_iter,
                                 nb_iter_sdf)
  Model_solved$Model$kappa_pi <- 0 # consider nominal ZC bonds.
  Model_solved$Model$kappa_y  <- 0 # consider nominal ZC bonds.
  res_bond_prices <- compute_bond_prices(Model_solved,maxH=10,nb_iter_sdf)
  p <- compute_uncond_distri(Model_solved$indicators_x,
                             Model_solved$Probas,nb_iter = 2000)
  avg_yds <- t(p) %*% res_bond_prices$all_rth
  lines(1:maxH,c(avg_yds),col="blue",lwd=2,lty=2)

  grid()
  
  if(regime=="Supply"){
    legend("topright", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend
           c("Nominal bonds (w/o credit risk)","ILB (w/o credit risk)","GDP-LB (w/o credit risk)","Nominal bonds, with credit risk"),
           lty=c(1,1,1,2), # gives the legend appropriate symbols (lines)
           lwd=c(2,2,2,2), # line width
           col=c("blue","red","dark grey","blue"),
           bg="white",
           #pch=c(3,NaN,NaN),
           seg.len = 3,
           cex=1,title = "Annualized expected nominal returns")
    
  }
}

dev.off()




# Debt issuance strategies: ----------------------------------------------------
# Prepare table showing performances of debt issuance strategies
#    for different values of chi
# ------------------------------------------------------------------------------

values_of_chi <- c(.2,.9)

grids <- make_grid(nb_grid = nb_grid,
                   min_d = 0.2,
                   max_d = 1.6,
                   min_rr = 0,
                   max_rr=.15,
                   sigma_eps = Model$sigma_eps,
                   all_quantiles_eps = c(-2,-1,1,2))

extensions       <- c("_nominal","_ILB","_GDPLB")
names_extensions <- c("Nominal","ILB","GDP-LB")

outputs <- c("mean_d","stdv_d","DaR95","mean_rr","stdv_rr",
             "stdv_Delta_d","avg_PD[maxH]","avg_spreads[maxH]")
latex.column.names <- c("$\\mathbb{E}(d)$","$\\sqrt{\\mathbb{V}(d)}$",
                        "$q_{95}(d)$",
                        "$\\mathbb{E}(r)$","$\\sqrt{\\mathbb{V}(r)}$",
                        "$\\sqrt{\\mathbb{V}(\\Delta d)}$",
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
                     paste("\\caption{Performances of debt issuance strategies in stylized versions of the model, $\\mu_\\eta=",
                           elasticity_of_surpluses,"\\times \\mu_y$ and ",
                           "$\\nu_y = ",ifelse(abs_nu_y<0,"-",""),abs_nu_y,"$",
                           "}",sep=""),
                     paste("\\label{tab:DemSup_elast",
                           elasticity_of_surpluses,
                           "_nu",abs_nu_y,
                           "}",sep=""),
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
                                        regime,",grids,nb_iter,nb_iter_sdf)",sep=""))))
    
    strat_nominal <- run_strategy(Models$Model_solved_nominal,
                                  maxH,
                                  nb_iter_sdf = nb_iter_sdf)
    strat_ILB     <- run_strategy(Models$Model_solved_ILB,
                                  maxH,
                                  nb_iter_sdf = nb_iter_sdf)
    strat_GDPLB   <- run_strategy(Models$Model_solved_GDPLB,
                                  maxH,
                                  nb_iter_sdf = nb_iter_sdf)
    
    M <- rbind(c(strat_nominal$mean_d,strat_ILB$mean_d,strat_GDPLB$mean_d),
               c(strat_nominal$stdv_d,strat_ILB$stdv_d,strat_GDPLB$stdv_d),
               c(strat_nominal$DaR95,strat_ILB$DaR95,strat_GDPLB$DaR95),
               c(strat_nominal$mean_rr,strat_ILB$mean_rr,strat_GDPLB$mean_rr),
               c(strat_nominal$stdv_rr,strat_ILB$stdv_rr,strat_GDPLB$stdv_rr),
               c(strat_nominal$stdv_Delta_d,strat_ILB$stdv_Delta_d,strat_GDPLB$stdv_Delta_d),
               c(strat_nominal$avg_PD[maxH],strat_ILB$avg_PD[maxH],strat_GDPLB$avg_PD[maxH]),
               c(strat_nominal$avg_spreads[maxH],strat_ILB$avg_spreads[maxH],strat_GDPLB$avg_spreads[maxH]))
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
                     "\\parbox{\\linewidth}{\\textit{Notes}: This table shows performance metrics associated with three different debt issuance strategies; each strategy consists in issuing a given type of perpetuities: a nominal perpetuity ($\\kappa_\\pi=0$ and $\\kappa_y=0$), an inflation-indexed perpetuity nominal ($\\kappa_\\pi=1$ and $\\kappa_y=0$), and a GDP-indexed perpetuity nominal ($\\kappa_\\pi=1$ and $\\kappa_y=1$). We consider two different values of $\\chi$ (the higher $\\chi$, the higher the average debt maturity). '$d$' denotes the debt-to-GDP ratio. '$r$' denotes the debt service, including debt indexation (in percent of GDP). '$\\sqrt{\\mathbb{V}(x)}$' corresponds to the standard deviation of variable $x$; '$PD$' stands for '10-year probability of default' (expressed in percent); '$spd$' stands for '10-year credit spread' (expressed in basis point), '$q_{95}(d)$' is the $95^{th}$ percentile of the debt-to-GDP distribution.}",
                     "\\end{footnotesize}",
                     "\\end{table}")

name.of.file <- "table_DemSup"
latex.file <- paste(name.of.file,
                    "_elastsurplus",elasticity_of_surpluses,
                    "_nu",abs_nu_y,
                    ".txt", sep="")
write(latex_table, paste("tables/",latex.file,sep=""))




# ------------------------------------------------------------------------------
# Prepare plot showing effect of credit risk on nominal yield curves
# ------------------------------------------------------------------------------

if(indic_addit_figure==1){
  FILE = paste("figures/Figure_yds_curve_wCredit_DemaSupp.pdf",sep="")
  pdf(file=FILE, pointsize=10, width=9, height=5)
  
  par(mfrow=c(1,2))
  par(plt=c(.12,.95,.1,.85))
  
  for(regime in c("Demand","Supply")){
    
    eval(parse(text = gsub(" ","",paste("Model <- Model_",regime,sep=""))))
    
    # Nominal bonds:
    Model$kappa_pi <- 0
    Model$kappa_y  <- 0
    
    # Without credit risk: -------------------------------------------------------
    
    RES <- prepare_returns_yds(Model,maxH)
    
    ylim <- c(min(RES$avg_annual_nominal_returns)-.005,
              max(RES$avg_annual_nominal_returns)+.01)
    
    plot(1:maxH,rep(0,maxH),las=1,
         xlab="",ylab="",type="l",ylim=ylim,lwd=2,col="white",
         main=paste(ifelse(regime=="Demand","(a) ","(b) "),regime,"-driven economy",sep=""))
    
    lines(1:maxH,c(RES$avg_annual_nominal_returns),col="blue",lwd=2)
    
    
    # With credit risk: ----------------------------------------------------------
    
    # First case, nu_y = 0
    Model$nu_y <- 0
    Model_solved <- solve_ToyModel(Model,
                                   grids,nb_iter,
                                   nb_iter_sdf)
    res_bond_prices <- compute_bond_prices(Model_solved,maxH=10,nb_iter_sdf)
    p <- compute_uncond_distri(Model_solved$indicators_x,
                               Model_solved$Probas,nb_iter = 2000)
    avg_yds <- t(p) %*% res_bond_prices$all_rth
    lines(1:maxH,c(avg_yds),col="blue",lwd=2,lty=2)
    
    # Second case, nu_y < 0
    Model$nu_y <- - 0.1
    Model_solved <- solve_ToyModel(Model,
                                   grids,nb_iter,
                                   nb_iter_sdf)
    res_bond_prices <- compute_bond_prices(Model_solved,maxH=10,nb_iter_sdf)
    p <- compute_uncond_distri(Model_solved$indicators_x,
                               Model_solved$Probas,nb_iter = 2000)
    avg_yds <- t(p) %*% res_bond_prices$all_rth
    lines(1:maxH,c(avg_yds),col="blue",lwd=2,lty=3)
    
    grid()
    
    if(regime=="Supply"){
      legend("topright", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend
             c("no credit risk",
               "with credit risk, no macro impact of default",
               "with credit risk, with macro impact of default"),
             lty=c(1,2,3), # gives the legend appropriate symbols (lines)
             lwd=c(2,2,2), # line width
             col=c("blue","blue","blue"),
             bg="white",
             #pch=c(3,NaN,NaN),
             seg.len = 3,
             cex=1)
    }
  }
  
  dev.off()
}



FILE = paste("figures/Figure_nu_effect_DemaSupp.pdf",sep="")
pdf(file=FILE, pointsize=12, width=9, height=6)

par(mfrow=c(2,2))
#par(plt=c(.12,.95,.1,.85))
par(plt = c(.15,.95,.15,.8))

for(regime in c("Demand","Supply")){
  
  eval(parse(text = gsub(" ","",paste("Model <- Model_",regime,sep=""))))
  
  # The government issues nominal bonds:
  Model$kappa_pi <- 0
  Model$kappa_y  <- 0
  
  for(nu_y in c(0,-.1)){
    Model$nu_y <- nu_y
    Model_solved <- solve_ToyModel(Model,
                                   grids,nb_iter,
                                   nb_iter_sdf)
    
    indic <- which(Model_solved$d==Model_solved$d_1)
    
    eval(parse(text = gsub(" "," ",paste("main.t <- expression(paste('",regime,"-driven economy, ',nu[y],' = ',",nu_y,",sep=''))",sep=""))))

    plot(Model_solved$d[indic],Model_solved$q0[indic],
         col="#AAAAAA22",pch=19,
         ylim=c(0,.10),
         main=main.t,
         xlab="Debt-to-GDP",ylab="Yield-to-maturity",las=1)
    points(Model_solved$d[indic],Model_solved$q[indic],col="#5555AA22",pch=17)
    
    grid()
    
    if((regime=="Supply")&(nu_y==0)){
      legend("bottomleft", # places a legend at the appropriate place c("Health","Defense"), # puts text in the legend
             c("Risk-free perpetuity",
               "Defaultable perpetuity"),
             lty=c(NaN,NaN), # gives the legend appropriate symbols (lines)
             lwd=1, # line width
             col=c("#AAAAAA77","#5555AA99"),
             bg="white",
             pch=c(19,17),
             seg.len = 2,
             cex=1)
    }
  }
}

dev.off()

