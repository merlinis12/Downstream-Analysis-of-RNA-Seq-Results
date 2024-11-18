# Ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# List of CRAN packages
cran_packages <- c(
  "readr", 
  "formattable", 
  "ggplot2", 
  "data.table", 
  "dplyr", 
  "pander", 
  "ggrepel", 
  "ggnewscale"
)

# Install missing CRAN packages
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# List of Bioconductor packages
bioc_packages <- c(
  "fgsea", 
  "msigdbr", 
  "org.Hs.eg.db"
)

# Install missing Bioconductor packages
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Load all packages
all_packages <- c(cran_packages, bioc_packages)

for (pkg in all_packages) {
  library(pkg, character.only = TRUE)
}

message("All required packages have been successfully installed and loaded!")
