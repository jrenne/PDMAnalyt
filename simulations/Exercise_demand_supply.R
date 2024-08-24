# ==============================================================================
# This script runs the demand/supply exercise
# ==============================================================================

# Create stylized models

low_pi <- .0
med_pi <- .03
hig_pi <- .06

low_y <- .0
med_y <- .02
hig_y <- .04

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




# Debt issuance strategies: ----------------------------------------------------
# Prepare table showing performances of debt issuance strategies
#    for different values of chi
# ------------------------------------------------------------------------------

grids <- make_grid(nb_grid = nb_grid,
                   min_d  = min_d,
                   max_d  = max_d,
                   min_rr = min_rr,
                   max_rr = max_rr,
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

