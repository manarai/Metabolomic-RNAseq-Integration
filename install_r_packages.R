#!/usr/bin/env Rscript
# One-shot installer for the R packages used by the v4 / v2 notebooks via %%R magic.
# Run once after installing R itself:
#     Rscript install_r_packages.R
#
# Packages installed:
#   CRAN:         BiocManager, r.jive
#   Bioconductor: limma, edgeR, mixOmics  (DIABLO is in mixOmics)

CRAN_REPO <- "https://cloud.r-project.org"

cat("=== Installing R dependencies for maize-multiomics ===\n\n")

# Step 1 — BiocManager (needed for Bioconductor packages)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager from CRAN...\n")
  install.packages("BiocManager", repos = CRAN_REPO)
}

# Step 2 — Bioconductor packages
bioc_pkgs <- c("limma", "edgeR", "mixOmics")
to_install_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install_bioc) > 0) {
  cat(sprintf("Installing Bioconductor packages: %s\n",
              paste(to_install_bioc, collapse = ", ")))
  BiocManager::install(to_install_bioc, update = FALSE, ask = FALSE)
} else {
  cat("All Bioconductor packages already installed.\n")
}

# Step 3 — CRAN-only packages
cran_pkgs <- c("r.jive")
to_install_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install_cran) > 0) {
  cat(sprintf("Installing CRAN packages: %s\n",
              paste(to_install_cran, collapse = ", ")))
  install.packages(to_install_cran, repos = CRAN_REPO)
} else {
  cat("All CRAN packages already installed.\n")
}

# Verify
cat("\n=== Verification ===\n")
all_pkgs <- c("BiocManager", "limma", "edgeR", "mixOmics", "r.jive")
for (p in all_pkgs) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING")
  cat(sprintf("  %-15s %s\n", p, v))
}

cat("\nDone. Activate the conda environment and launch the notebook:\n")
cat("  conda activate maize-multiomics\n")
cat("  jupyter lab multiomics_integration_v4.ipynb\n")
