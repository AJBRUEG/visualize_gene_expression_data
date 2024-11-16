#script to capture, process, analyze, and vizualize transcriptomic data

#-----setup-----

# Install necessary packages
install.packages(c("tidyverse", "httr", "jsonlite", "pheatmap"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("GEOquery", "limma", "edgeR"))

#-----get data-----

if(!require("GEOquery"))install.packages("GEOquery")
library(GEOquery)

# Fetch dataset from GEO (example: GSEXXXXX)
#specify a GEO dataset ID (GSExxxxxx)
#look through the GEO database for a dataset you are interested in
gse_id <- "GSE223667"  # Replace with a plant-specific GEO dataset ID
gse_data <- getGEO(GEO = gse_id, GSEMatrix = TRUE, getGPL = FALSE)

# Extract expression data
expression_set <- gse_data[[1]]
expression_data <- exprs(expression_set)

# View metadata
metadata <- pData(expression_set)
