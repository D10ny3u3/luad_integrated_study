anchors <- FindTransferAnchors(
  reference = soil_epi_selected,
  query = xiong_mia_selected,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:50
)

soil_epi_selected <- RunUMAP(soil_epi_selected, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", return.model = TRUE)

DimPlot(soil_epi_selected, group.by = "seurat_clusters")

xiong_mia_selected <- MapQuery(
  anchorset = anchors,
  query = xiong_mia_selected,
  reference = soil_epi_selected,
  refdata = list(
    seurat_clusters = "seurat_clusters",
    celltype = "annotation",
    annotation_note = "annotation_note"),
  transferdata.args = list(
    prediction.assay = TRUE
  ),
  reference.reduction = "pca", 
  reduction.model = "umap.harmony"
)

p <- DimPlot(
  xiong_mia,
  reduction = "ref.umap", 
  group.by = "predicted.celltype", 
  combine     = TRUE
) +
  ggtitle("Query MIA") + 
  labs(
    x = "UMAP_1",
    y = "UMAP_2"
  ) +
  theme(
    axis.title = element_text(size = 16, face = "bold"), # axis title bold and large
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  # plot title bold and large
  ) 
p

ggsave(
  "UMAP_query_MIA.pdf",
  p,
  width=5,
  height=4
)
