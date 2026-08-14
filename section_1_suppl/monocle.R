# data ----

LUAD_tar <- merge(LUAD_ia, LUAD_mia)
LUAD_main <- subset(LUAD_tar, cluster_names == "at2_adca")

# process ----
LUAD_main <- JoinLayers(LUAD_main)
LUAD_main[["RNA"]] <- split(LUAD_main[["RNA"]], f = LUAD_main$orig.ident)
LUAD_main <- NormalizeData(LUAD_main)
LUAD_main <- FindVariableFeatures(LUAD_main)
LUAD_main <- ScaleData(LUAD_main)
LUAD_main <- ScaleData(LUAD_main, features = rownames(LUAD_main))
LUAD_main <- RunPCA(LUAD_main)
library(harmony)
LUAD_main <- RunHarmony(LUAD_main, "orig.ident")
LUAD_main <- FindNeighbors(LUAD_main, reduction = "harmony", dims = 1:30)
LUAD_main <- FindClusters(LUAD_main, resolution = 0.5, cluster.name = "harmony_clusters")
LUAD_main <- RunUMAP(LUAD_main, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

p<-DimPlot(LUAD_main, group.by = "seurat_clusters", label = T)
ggsave("umap_seurat.pdf", p, height = 4, width = 5)
dim(LUAD_main)

# explore data ----
FeaturePlot(LUAD_main, "SCGB3A2")
FeaturePlot(LUAD_main, "CRABP2")
FeaturePlot(LUAD_main, "SFTPC")
FeaturePlot(LUAD_main, "MARCKSL1")
FeaturePlot(LUAD_main, "nFeature_RNA")
FeaturePlot(LUAD_main, "SFTPC")
FeaturePlot(LUAD_main, "ROS1")
FeaturePlot(LUAD_main, "NEAT1")
DimPlot(LUAD_main, split.by = "inv_state")
p<-DimPlot(LUAD_main, group.by = "predicted.celltype" , split.by = "inv_state")
ggsave("umap_celltype_state.pdf", p, height = 4, width = 8)

VlnPlot(LUAD_main, "SFTPC", group.by = "predicted.celltype")
VlnPlot(LUAD_main, "SCGB3A2", group.by = "predicted.celltype")

# deg ----

LUAD_main <- JoinLayers(LUAD_main)
deg_LUAD_main <- FindAllMarkers(LUAD_main, only.pos = TRUE)
# saveRDS(deg_LUAD_main, "deg_LUAD_main.rds")
library(openxlsx)
wb <- createWorkbook()
clusters <- unique(deg_LUAD_main$cluster)
for (cluster in clusters) {
  cluster_data <- deg_LUAD_main %>% dplyr::filter(cluster == !!cluster)
  addWorksheet(wb, sheetName = paste0("C", cluster))
  writeData(wb, sheet = paste0("C", cluster), cluster_data)
}
saveWorkbook(wb, "deg_LUAD_main_clusters.xlsx", overwrite = TRUE)

# monocle 2 ----

library(monocle)
library(igraph)

expr_matrix <- as(as.matrix(LUAD_main@assays$RNA@layers$counts), 'sparseMatrix')
p_data <- LUAD_main@meta.data 
table(LUAD_main$predicted.celltype)
p_data$celltype <- LUAD_main$predicted.celltype
f_data <- data.frame(gene_short_name = row.names(LUAD_main), row.names = row.names(LUAD_main))
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

ordergene <- VariableFeatures(LUAD_main)
cds <- setOrderingFilter(cds, ordergene)
plot_ordering_genes(cds)
# reduce dim
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

## order cells ----
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 3)

plot_cell_trajectory(cds)

## vis ----

str(cds$celltype)
cds$celltype <- factor(cds$celltype)
levels(unique(cds$celltype))
new_order <- c("AT2_like", "SCGB3A2_pos", "CRABP2_pos")
cds$celltype <- factor(
  cds$celltype,
  levels = new_order
)
cds
plot_cell_trajectory(cds, color_by = "celltype")
plot_cell_trajectory(cds, color_by = "annotation")

cds$annotation <- LUAD_main$annotation
p <- plot_cell_trajectory(cds, color_by = "celltype") +
  facet_wrap(~seurat_clusters, nrow = 1)
ggsave("monocle_celltype_seurat.pdf", p, height = 4, width = 10)


p <- plot_cell_trajectory(cds, color_by = "annotation") +
  facet_wrap(~annotation, nrow = 1)
p

plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "Pseudotime")

max(cds$Pseudotime)
min(cds$Pseudotime)

## plot_genes_in_pseudotime ----

to_be_tested <- row.names(
  subset(fData(cds),
         gene_short_name %in% c("SFTPC", "SCGB3A2", "HPGD", "CRABP2", "MARCKSL1")))
cds_subset <- cds[to_be_tested,]
plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")


## plot single gene ----
colnames(pData(cds))
pData(cds)$SCGB3A2 =log2(exprs(cds)['SCGB3A2',]+1)
pData(cds)$SFTPC =log2(exprs(cds)['SFTPC',]+1)
pData(cds)$CRABP2 =log2(exprs(cds)['CRABP2',]+1)
library(ggsci)
p1=plot_cell_trajectory(cds, color_by="SCGB3A2") + scale_color_gsea()
p1
plot_cell_trajectory(cds, color_by="SFTPC") + scale_color_gsea()
plot_cell_trajectory(cds, color_by="CRABP2") + scale_color_gsea()


table(LUAD_main$annotation, LUAD_main$predicted.celltype)
