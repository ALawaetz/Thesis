# Install and load the AnnotationHub package if it's not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Check if the AnnotationHub package is installed, and install it if not
if (!requireNamespace("AnnotationHub", quietly = TRUE)) {
  BiocManager::install("AnnotationHub")
}

if (!requireNamespace("AnnotationForge", quietly = TRUE)) {
  BiocManager::install("AnnotationForge")
}

# Load the AnnotationHub library
library(AnnotationForge)
library(AnnotationHub)

### this command takes a lot of time to run
### and required a lot of disk space (at least 50 GB).
### the majority can be deleted afterwards
makeOrgPackageFromNCBI(version = "0.1",
                       author = "Andreas Lawaetz",
                       maintainer = "Andreas Lawaetz <acl58@bath.ac.uk>",
                       outputDir = "/Users/andreas/Streptomyces/orgDB",
                       tax_id = "100226",  # Tax ID for Streptomyces coelicolor
                       genus = "Streptomyces",
                       species = "coelicolor")
