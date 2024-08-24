
# ------------------------------------------------------------------------------
# Prepare Latex table showing specifications (Demand / Supply exercise)
# ------------------------------------------------------------------------------

nb_m <- dim(Model$Omega)[1]

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
                           # "$\\nu_y=",make.entry(Model$nu_y,format.nb = format.nb2,dollar=0),"$, ",
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
# (Demand / Supply exercise)
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

