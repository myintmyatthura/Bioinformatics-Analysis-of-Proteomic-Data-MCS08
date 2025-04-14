
# Create a project-specific library directory
lib_dir <- "/srv/shiny-server/Bioinformatics-Analysis-of-Proteomic-Data-MCS08/r_libs"
if (!dir.exists(lib_dir)) {
  dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
}

# Define required packages
if (!"BiocManager" %in% installed.packages()[,"Package"]) {
  install.packages("BiocManager", repos = "https://cran.rstudio.com/")
}
required_packages <- c("shiny", "ggplot2", "dplyr", "testthat", "shinythemes",
                       "pheatmap", "RColorBrewer", "httr", "jsonlite",
                       "pool", "DBI", "tinytex", "rmarkdown", "knitr","kableExtra")
bioc_packages <- c("limma")


# Install missing packages to the project library
for (pkg in required_packages) {
  if (!pkg %in% installed.packages(lib.loc = lib_dir)[,"Package"]) {
    tryCatch({
      install.packages(pkg, lib = lib_dir, repos = "https://cran.rstudio.com/", dependencies = TRUE)
      cat(paste0("Successfully installed ", pkg, "\n"))
    }, error = function(e) {
      cat(paste0("Failed to install ", pkg, ": ", e$message, "\n"))
    })
  }
}

if (!tinytex::is_tinytex()) {
  cat("Installing TinyTeX distribution...\n")
  tinytex::install_tinytex()
}

# Add this library to the search path
.libPaths(c(lib_dir, .libPaths()))
