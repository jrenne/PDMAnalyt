

# ==============================================================================
# Prepare table showing moment matching results
# ==============================================================================

latex_table <- rbind("\\begin{table}[ph!]",
                     "\\caption{Model-implied versus targeted moments}",
                     "\\label{tab:moment_matching}",
                     "\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}rr}",
                     "\\hline",
                     "Moment&Model&Target\\\\",
                     "\\hline")

table4Latex <- compare_model_target_moments(Model,targets)$table_compare_moments
Rownames <- rownames(table4Latex)

for(i in 1:dim(table4Latex)[1]){
  this_line <- paste(Rownames[i],"&",make.entry(table4Latex[i,1],format.nb3),
                     "&",make.entry(table4Latex[i,2],format.nb3),
                     "\\\\",sep="")
  latex_table <- rbind(latex_table,
                       this_line)
}

latex_table <- rbind(latex_table,
                     "\\hline",
                     "\\end{tabular*}",
                     "\\begin{footnotesize}",
                     "\\parbox{\\linewidth}{\\textit{Notes}: This table compares model-implied with targeted moments. The distance between these moments is part of the loss function that is minimized to estimate the components of $\\mu_\\pi$, $\\mu_y$, and $\\Omega$. See Subsection~\\ref{sub:calibration} for more details.}",
                     "\\end{footnotesize}",
                     "\\end{table}")

name.of.file <- "table_moment_matching"
latex.file <- paste(name.of.file,".txt", sep="")
write(latex_table, paste("tables/",latex.file,sep=""))


