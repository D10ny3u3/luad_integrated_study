# code time: 20260603
# ref https://cole-trapnell-lab.github.io/monocle-release/docs/
# ref https://www.jianshu.com/p/5d6fd4561bc0

library(monocle)
library(igraph)

table(kim_epi_tumor$annotation)

kim_epi_tumor_selected <- subset(kim_epi_tumor, annotation %in% c(
  "CRABP2_pos", "AT2_like", "SCGB3A2_pos"
))

expr_matrix <- as(as.matrix(kim_epi_tumor_selected@assays$RNA@layers$counts), 'sparseMatrix')
p_data <- kim_epi_tumor_selected@meta.data 
p_data$celltype <- kim_epi_tumor_selected$annotation
f_data <- data.frame(gene_short_name = row.names(kim_epi_tumor_selected), row.names = row.names(kim_epi_tumor_selected))
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)

cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
rm(expr_matrix, p_data, f_data, pd, fd)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)

expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10))

diff_test_res <- differentialGeneTest(
  cds[expressed_genes,],
  fullModelFormulaStr = "~seurat_clusters")

saveRDS(diff_test_res, "diff_test_res_by_seurat_clusters.rds")

# ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

diff_test_res_selected <- subset(diff_test_res, qval < 0.01)
# diff_test_res_selected <- subset(diff_test_res_by_seurat_clusters, qval < 0.01)
diff_test_res_selected <- diff_test_res_selected[order(diff_test_res_selected$qval,decreasing=F),]
ordergene <- row.names(diff_test_res_selected)[order(diff_test_res_selected$qval)][1:3000]

# alternative ordergene
kim_epi_tumor_selected <- JoinLayers(kim_epi_tumor_selected)
kim_epi_tumor_selected[["RNA"]] <- split(kim_epi_tumor_selected[["RNA"]], f = kim_epi_tumor_selected$orig.ident)
kim_epi_tumor_selected <- NormalizeData(kim_epi_tumor_selected)
kim_epi_tumor_selected <- FindVariableFeatures(kim_epi_tumor_selected)
ordergene <- VariableFeatures(kim_epi_tumor_selected)

cds <- setOrderingFilter(cds, ordergene)
# plot_ordering_genes(cds)

cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
cds <- orderCells(cds)
plot_cell_trajectory(cds)

str(cds$celltype)
cds$celltype <- factor(cds$celltype)
levels(unique(cds$celltype))
new_order <- c("AT2_like", "SCGB3A2_pos", "CRABP2_pos")
cds$celltype <- factor(
  cds$celltype,
  levels = new_order
)

plot_cell_trajectory(cds, color_by = "celltype")
plot_cell_trajectory(cds, color_by = "celltype") +
  facet_wrap(~celltype, nrow = 1)

plot_cell_trajectory(cds, color_by = "State")

# occasionally
cds <- orderCells(cds, root_state = 3)

# modify pseudotime according to biological knowledge 

max(cds$Pseudotime)
# 20.05734
min(cds$Pseudotime)
cds$Pseudotime <- 20.05734 - cds$Pseudotime 
plot_cell_trajectory(cds, color_by = "Pseudotime")

saveRDS(cds, "kim_epi_tumor_selected_cds.rds")

# plot_genes_in_pseudotime ----

to_be_tested <- row.names(
  subset(fData(cds),
         gene_short_name %in% c("SFTPC", "SCGB3A2", "HPGD", "CRABP2", "MARCKSL1")))
cds_subset <- cds[to_be_tested,]
plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")

# BEAM ----

kim_epi_tumor_selected$State <- cds@phenoData@data$State
table(kim_epi_tumor_selected$State)
# DimPlot(kim_epi_tumor_selected, group.by = "State")
kim_epi_tumor_selected@active.ident <- kim_epi_tumor_selected$State
kim_epi_tumor_state_deg <- FindMarkers(
  kim_epi_tumor_selected, "1", "2")

top_genes <- kim_epi_tumor_state_deg %>%
  arrange(abs(p_val_adj)) %>%
  slice(1:80)
top_genes <- rownames(top_genes)
top_genes <- top_genes[!grepl("^RP", top_genes)]

library(colorRamps)
library(pheatmap)

plot_genes_branched_heatmap_new(
  cds[top_genes, ],
  branch_point = 1,
  num_clusters = 2,
  cores = 1,
  use_gene_short_name = T,
  show_rownames = T)
