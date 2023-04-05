#get correlation between seurat clusters
library(Seurat)
library(dplyr)
library(ggcorrplot)

get_average_correlation <- function(object, nfeatures=2000) {
  
  # Normalize data and find variable features
  object <- object %>% 
    NormalizeData() %>% 
    FindVariableFeatures(nfeatures = nfeatures)
  
  # Calculate average expression values
  average <- AverageExpression(object = object, slot = "data", group.by = "annotation", 
                               features = object@assays$RNA@var.features)$RNA
  
  # Calculate correlation
  correlation <- cor(average, use = "pairwise.complete.obs", method = "spearman")
  
  return(correlation)
}

#plotting corrplot
corr <- get_average_correlation(object = object, nfeatures = 2000)
ggcorrplot(corr, hc.order = TRUE, outline.color = 'white',tl.srt = 90)+
  scale_fill_gradient2(high="#C25539", mid="white", low = "#3F7F93",
                       breaks=c(0.5, 0.75, 1), limit=c(0.5, 1), midpoint = 0.75) 

