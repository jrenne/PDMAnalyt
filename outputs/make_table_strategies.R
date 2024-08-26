

min_chi <- min(parameters[,"chi"])
max_chi <- max(parameters[,"chi"])

min_kappa_pi <- min(parameters[,"kappa_pi"])
max_kappa_pi <- max(parameters[,"kappa_pi"])

min_kappa_y <- min(parameters[,"kappa_y"])
max_kappa_y <- max(parameters[,"kappa_y"])

matrix_values_4table <- NULL
matrix_values_4table <- rbind(matrix_values_4table,
                              c(min_chi,0,0))
matrix_values_4table <- rbind(matrix_values_4table,
                              c(min_chi,0,max_kappa_y))
matrix_values_4table <- rbind(matrix_values_4table,
                              c(min_chi,max_kappa_pi,0))
matrix_values_4table <- rbind(matrix_values_4table,
                              c(min_chi,max_kappa_pi,max_kappa_y))
matrix_values_4table <- rbind(matrix_values_4table,
                              c(max_chi,0,0))
matrix_values_4table <- rbind(matrix_values_4table,
                              c(max_chi,0,max_kappa_y))
matrix_values_4table <- rbind(matrix_values_4table,
                              c(max_chi,max_kappa_pi,0))
matrix_values_4table <- rbind(matrix_values_4table,
                              c(max_chi,max_kappa_pi,max_kappa_y))

matrix_values_4table <- rbind(matrix_values_4table,
                               parameters[which.min(M[,"stdv_d"]),])
matrix_values_4table <- rbind(matrix_values_4table,
                              parameters[which.min(M[,"DaR95"]),])
matrix_values_4table <- rbind(matrix_values_4table,
                               parameters[which.min(M[,"avg_PD[maxH]"]),])

colnames(matrix_values_4table) <- c("chi","kappa_pi","kappa_y")

# names_cases <- NULL
# for(i in 1:dim(matrix_values_4table)[1]){
#   names_cases <- c(names_cases,
#                    paste(
#                      "$\\chi=",make.entry(matrix_values_4table[i,"chi"],format.nb1,dollar=0),"$, ",
#                      "$\\kappa_\\pi=",make.entry(matrix_values_4table[i,"kappa_pi"],format.nb1,dollar=0),"$, ",
#                      "$\\kappa_y=",make.entry(matrix_values_4table[i,"kappa_y"],format.nb1,dollar=0),"$",
#                      sep=""))
# }
names_cases <- NULL
for(i in 1:dim(matrix_values_4table)[1]){
  names_cases <- c(names_cases,
                   paste(
                     "$(",
                     make.entry(matrix_values_4table[i,"chi"],format.nb1,dollar=0),",",
                     make.entry(matrix_values_4table[i,"kappa_pi"],format.nb1,dollar=0),",",
                     make.entry(matrix_values_4table[i,"kappa_y"],format.nb1,dollar=0),")$",
                     sep=""))
}

columns <- "c" # empty column
for(i in 1:length(latex.column.names)){
  columns <- paste(columns,"c",sep="")
}
column_names <- "$(\\chi,\\kappa_\\pi,\\kappa_y)$&"
for(i in 1:length(latex.column.names)){
  column_names <- paste(column_names,"&",latex.column.names[i],sep="")
}

latex_table <- rbind("\\begin{table}[ph!]",
                     paste("\\caption{Performances of debt issuance strategies in the calibrated model",
                           #", $\\mu_\\eta=",elasticity_of_surpluses,"\\times \\mu_y$",
                           "}",sep=""),
                     paste("\\label{tab:Strategies",
                           #"_elast",elasticity_of_surpluses,
                           "}",sep=""),
                     paste("\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}",columns,"}",sep=""),
                     "\\hline",
                     paste(column_names,"\\\\",sep=""),
                     "\\hline")

count_case <- 0
for(case in 1:dim(matrix_values_4table)[1]){
  count_case <- count_case + 1
  indic_in_M <- which((matrix_values_4table[case,"chi"]==parameters[,"chi"])&
                        (matrix_values_4table[case,"kappa_pi"]==parameters[,"kappa_pi"])&
                        (matrix_values_4table[case,"kappa_y"]==parameters[,"kappa_y"]))
  this.line <- paste(names_cases[count_case],"&")
  for(output in outputs){
    this.line <- paste(this.line,"&",make.entry(M[indic_in_M,output],format.nb2))
  }
  this.line <- paste(this.line,"\\\\")
  if(case == dim(matrix_values_4table)[1] - 2){
    latex_table <- rbind(latex_table,
                         "\\hline",
                         this.line)
  }else{
    latex_table <- rbind(latex_table,this.line)
  }
}

latex_table <- rbind(latex_table,
                     "\\\\",
                     "\\hline",
                     "\\end{tabular*}",
                     "\\begin{footnotesize}",
                     "\\parbox{\\linewidth}{\\textit{Notes}: This table shows performance metrics associated with different debt issuance strategies characterized by the issuance of perpetuities of different durations (captured by the coupon decay rate $\\chi$), a coefficient of indexation to inflation $\\kappa_\\pi$ and a coefficient of indexation to real GDP $\\kappa_y$. The model is the one whose parameterization is reported in Table~\\ref{tab:param}. '$d$' denotes the debt-to-GDP ratio. '$r$' denotes the debt service, including debt indexation (in percent of GDP). '$\\sqrt{\\mathbb{V}(x)}$' corresponds to the standard deviation of variable $x$; '$PD$' stands for '10-year probability of default' (expressed in percent); '$spd$' stands for '10-year credit spread' (expressed in basis point), '$q_{95}(d)$' is the $95^{th}$ percentile of the debt-to-GDP distribution. The last three rows show the performances of the strategies implying the lowest $\\sqrt{\\mathbb{V}(d)}$, $q_{95}(d)$, and $\\mathbb{E}(PD)$, respectively.}",
                     "\\end{footnotesize}",
                     "\\end{table}")




name.of.file <- "table_strategies"
latex.file <- paste(name.of.file,
                    #"_elastsurplus",elasticity_of_surpluses,
                    ".txt", sep="")
write(latex_table, paste("tables/",latex.file,sep=""))


