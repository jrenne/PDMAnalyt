
# ==============================================================================
# Prepare table showing model parameterization
# ==============================================================================

nb_m <- dim(Model$Omega)[1]

columns <- NULL
for(i in 1:5){
  columns <- paste(columns,"r",sep="")
}
for(i in 1:nb_m){
  columns <- paste(columns,"c",sep="")
}

latex_table <- rbind("\\begin{table}[ph!]",
                     "\\caption{Model parameterization}",
                     "\\label{tab:param}",
                     paste("\\begin{tabular*}{\\textwidth}{c@{\\extracolsep{\\fill}}",columns,"}",sep=""),
                     "\\hline",
                     paste("Regime&$\\mu_\\pi$&&$\\mu_y$&&\\multicolumn{",nb_m,"}{c}{$\\Omega$}\\\\",sep=""),
                     "\\hline")

for(i in 1:nb_m){
  pi_i <- NULL
  for(j in 1:nb_m){
    pi_i <- paste(pi_i,"&",make.entry(Model$Omega[i,j],format.nb3),sep="")
  }
  this_line <- paste(i,"&",make.entry(Model$mu_pi[i],format.nb3),
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
                     "\\parbox{\\linewidth}{\\textit{Notes}: This table shows the model parameterization of the baseline model.",
                     paste("We also have: $\\alpha=",Model$alpha,"$,",
                           "$\\beta=",make.entry(Model$beta,format.nb = format.nb2,dollar=0),"$, ",
                           "$\\gamma=",make.entry(Model$gamma,format.nb = format.nb0,dollar=0),"$, ",
                           "$\\delta=",make.entry(Model$delta,format.nb = format.nb2,dollar=0),"$, ",
                           "$d^*=",make.entry(Model$d_star,format.nb = format.nb2,dollar=0),"$, ",
                           "$s^*=",make.entry(Model$s_star,format.nb = format.nb3,dollar=0),"$, ",
                           "$\\nu_y=",make.entry(Model$nu_y,format.nb = format.nb3,dollar=0),"$, ",
                           "$\\nu_\\pi=",make.entry(Model$nu_pi,format.nb = format.nb3,dollar=0),"$, ",
                           "$\\mu_\\eta=",Model$mu_eta[1]/Model$mu_y[1],"\\times \\mu_y$, ",
                           "$RR=",make.entry(Model$RR,format.nb = format.nb2,dollar=0),"$."),
                     "}",
                     "\\end{footnotesize}",
                     "\\end{table}")

name.of.file <- "table_param"
latex.file <- paste(name.of.file,".txt", sep="")
write(latex_table, paste("tables/",latex.file,sep=""))


