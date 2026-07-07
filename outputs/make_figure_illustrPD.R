# ==============================================================================
# Install illustrative default-probability chart.
# ==============================================================================

# This figure was produced with the original TikZ/LaTeX workflow used for the
# paper. We keep the publication-quality PDF as a package resource because a
# base-R fallback gives a visibly different result and the TikZ route depends on
# optional LaTeX packages that are not always installed.

message("Copying resources/formula.pdf to figures/formula.pdf.")

source_file <- file.path("resources", "formula.pdf")
destination_file <- file.path("figures", "formula.pdf")

if(!file.exists(source_file)){
  stop("Missing resources/formula.pdf.")
}

dir.create("figures", showWarnings = FALSE)
file.copy(from = source_file,
          to = destination_file,
          overwrite = TRUE)

message("Saved figures/formula.pdf.")
