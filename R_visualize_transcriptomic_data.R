#script to capture, process, analyze, and vizualize transcriptomic data

#-----setup-----

# Install packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "Biobase", "limma"))
install.packages(c("ggplot2", "pheatmap", "reshape2", "dplyr"))

# Load libraries
library(GEOquery)
library(Biobase)
library(limma)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(dplyr)


#-----get data-----

# Download the GEO dataset
gse <- getGEO("GSE15521", GSEMatrix = TRUE)

# Access the expression set
exprSet <- gse[[1]]

# Extract expression data and metadata
exprData <- exprs(exprSet)         # Expression matrix
sampleMetadata <- pData(exprSet)   # Sample metadata
featureData <- fData(exprSet)      # Feature metadata

# View data structure
head(exprData)
head(sampleMetadata)

#-----PROCESS DATA-----

# Log-transform the data if needed
exprData <- log2(exprData + 1)

# Filter for highly variable genes (top 500 by variance)
topVarGenes <- head(order(apply(exprData, 1, var), decreasing = TRUE), 500)
filteredData <- exprData[topVarGenes, ]

#-----VISUALIZE WITH PCA-----

# Transpose for PCA: rows = samples, columns = genes
pca <- prcomp(t(filteredData), scale. = TRUE)

# Create a PCA data frame
pcaData <- as.data.frame(pca$x)
pcaData$Sample <- rownames(pcaData)
pcaData$Group <- sampleMetadata$group  # Adjust according to metadata

# Plot PCA
ggplot(pcaData, aes(PC1, PC2, color = Sample)) +
  geom_point(size = 3) +
  labs(title = "PCA of Transcriptomic Data", x = "PC1", y = "PC2") +
  theme_minimal()

#-----VISUALIZE WITH HEATMAP-----

# Scale data for heatmap
scaledData <- t(scale(t(filteredData)))

# Create heatmap
pheatmap(scaledData, annotation_col = sampleMetadata, 
         main = "Heatmap of Top Variable Genes",
         color = colorRampPalette(c("blue", "white", "red"))(50))


