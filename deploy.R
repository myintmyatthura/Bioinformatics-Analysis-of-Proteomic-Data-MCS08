
# Create a project-specific library directory
lib_dir <- "/srv/shiny-server/Bioinformatics-Analysis-of-Proteomic-Data-MCS08/r_libs"
dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)

# Define required packages
required_packages <- c("shiny", "ggplot2", "dplyr", "testthat", "shinythemes",
                       "limma", "pheatmap", "RColorBrewer")


# Install missing packages to the project library
new_packages <- required_packages[!required_packages %in% installed.packages(lib.loc = lib_dir)[,"Package"]]
if(length(new_packages) > 0) {
  install.packages(new_packages, lib = lib_dir, repos = "https://cran.rstudio.com/")
}

# Add this library to the search path
.libPaths(c(lib_dir, .libPaths()))
