# =================================================================
# Codes generating illustrative charts.
# =================================================================

proba.def.Q.quarters.fction.of.FS <- function(values.fs,
                                              alpha,
                                              sigma.nu,
                                              nb.Q){
  # values.fs are the values of the (log) fiscal space
  E.exp.mAlpha.max <- pnorm(values.fs/sigma.nu) + exp(alpha*values.fs + alpha^2 * sigma.nu^2/2) *
    (1 - pnorm(values.fs/sigma.nu + sigma.nu * alpha))
  P <- 1 - E.exp.mAlpha.max^nb.Q
  return(P)
}


library(tikzDevice)
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}",
                               "\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}",
                               "\\usepackage{amssymb}"))
tikz("formula.tex", width = 5, height = 3, standAlone = TRUE,
     packages = c("\\usepackage{tikz}",
                  "\\usepackage[active,tightpage,psfixbb]{preview}",
                  "\\PreviewEnvironment{pgfpicture}",
                  "\\setlength\\PreviewBorder{0pt}",
                  "\\usepackage{amssymb}"))



min.x <- -2
max.x <- 2
x <- seq(min.x,max.x,by=.01)
alpha <- c(.05,.1,2)

sigma.nu <- .2
nb.Q <- 4 # number of quarters

par(mfrow=c(1,1))
par(plt=c(.13,.95,.25,.97))

plot(x,proba.def.Q.quarters.fction.of.FS(x,alpha[1],sigma.nu,nb.Q),lwd=3,type="l",las=1,ylim=c(0,1),
     xlab="Fiscal space ($d^* - d_t$)",
     ylab="Probability of default")

polygon(c(x[1]-1,0,0,x[1]-1),c(-1,-1,2,2),
        col=rgb(.85,.85,.85),border = NA)
lines(x,proba.def.Q.quarters.fction.of.FS(x,alpha[1],sigma.nu,nb.Q),col="black",lwd=3)
lines(x,proba.def.Q.quarters.fction.of.FS(x,alpha[2],sigma.nu,nb.Q),col="dark grey",lwd=3)
lines(x,proba.def.Q.quarters.fction.of.FS(x,alpha[3],sigma.nu,nb.Q),col="black",lwd=3,lty=2)

text(x = min.x,y=.7,labels = "Negative fiscal",pos=4)
text(x = min.x,y=.6,labels = "space ($d^* < d_t$)",pos=4)
text(x = min.x+.6*(max.x-min.x),y=.4,labels = "Positive fiscal",pos=4)
text(x = min.x+.6*(max.x-min.x),y=.3,labels = "space ($d^* > d_t$)",pos=4)

legend("topright",
       c("$\\alpha=0.05$","$\\alpha=0.1$","$\\alpha=2$"),
       lty=c(1,1,2), # gives the legend appropriate symbols (lines)       
       lwd=3, # line width
       col=c("black","dark grey","black"),
       #title="For dates:",
       seg.len=3,
       cex=1,
       bg="white")


dev.off()


tools::texi2pdf("formula.tex")
system(paste(getOption("pdfviewer"), "formula.pdf"))

# file.copy("formula.pdf",
#           paste("../","formula.pdf", sep=""),
#           overwrite=TRUE)

source_file <- "formula.pdf"
destination_folder <- "figures/"

# Use file.copy to copy the file
file_copied <- file.copy(from = source_file, to = destination_folder,
                         overwrite=TRUE)

# List all files in the folder
all_files <- list.files(path = getwd(), full.names = TRUE)

# Filter files that start with "formula"
files_to_delete <- all_files[grepl("^formula", basename(all_files))]

# Remove the filtered files
file.remove(files_to_delete)