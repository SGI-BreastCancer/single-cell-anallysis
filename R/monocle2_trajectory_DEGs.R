##code for monocle2
##I used following code from public resources; http://cole-trapnell-lab.github.io/monocle-release/docs/
library(Seurat)
library(monocle)
cds <- as.CellDataSet(icaf) %>% ##convert seurat object to monocle object
  estimateSizeFactors() %>%
  estimateDispersions() %>%
  detectGenes(min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10))
cds <- reduceDimension(cds, max_components = 2, num_dim = 2, 
                       reduction_method = 'tSNE',
                       residualModelFormulaStr = "~orig.ident", ##batch correction against sample
                       verbose = T)
cds <- clusterCells(cds, num_clusters = 6)


diff_test_res <- differentialGeneTest(cds[expressed_genes,],
                                      fullModelFormulaStr = "~celltype") 
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, max_components = 2,  
                       reduction_method = 'DDRTree')
cds <- orderCells(cds)

##ploting
plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "celltype")

##diffrential gene test
diff_test_res <- differentialGeneTest(cds[expressed_genes,],
                                      fullModelFormulaStr = "~Pseudotime")
diff_test_res <- diff_test_res[!grepl(diff_test_res$gene_short_name, pattern = '^MT-'),]

sig_gene_names <- row.names(subset(diff_test_res, qval < 1e-20))
plot_pseudotime_heatmap(cds[sig_gene_names,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T)
