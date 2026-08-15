# umap ----

DimPlot(bischoff_epi_tumor, label = T, group.by = "seurat_clusters", label.size = 5)
DimPlot(bischoff_epi_tumor, label = T, group.by = "orig.ident")

table(bischoff_epi_tumor$histotype)
bischoff_epi_tumor <- subset(bischoff_epi_tumor, histotype == "sarcomatoid", invert = T)

library(ggplot2)

p <- DimPlot(
  bischoff_epi_tumor,
  reduction   = "umap.harmony",
  group.by    = "histotype",
  # label = T,
  # label.size = 5,
  repel       = TRUE,
  combine     = TRUE
) +
  ggtitle("Bischoff ") + # Title descripting Cell Type
  labs(
    x = "UMAP_1",
    y = "UMAP_2"
  ) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5) 
  ) 
p

ggsave(
  "umap_bischoff.pdf",
  p,
  width=5.5,
  height=4
)

# markers ----

bischoff_epi_tumor$histotype <- factor(bischoff_epi_tumor$histotype)
levels(unique(bischoff_epi_tumor$histotype))
new_order <- c("lepedic", "acinar", "papillary", "micropapillary", "solid")
bischoff_epi_tumor$histotype <- factor(
  bischoff_epi_tumor$histotype,
  levels = new_order
)

Idents(bischoff_epi_tumor) <- bischoff_epi_tumor$histotype
markers <- FindAllMarkers(bischoff_epi_tumor, only.pos = T)

library(dplyr)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
bischoff_epi_tumor <- ScaleData(bischoff_epi_tumor, features = rownames(bischoff_epi_tumor))
DoHeatmap(bischoff_epi_tumor, features = top5$gene) + NoLegend()

# dotplot ----
DotPlot(bischoff_epi_tumor, features = unique(top5$gene)) + NoLegend() + RotatedAxis()
features = unique(top5$gene)

dot_plot <- DotPlot(
  bischoff_epi_tumor,
  features = features, 
  assay = "RNA"
) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 14),
    axis.text = element_text(size = 16, face = "bold")
  ) +
  labs(y = NULL, x = NULL) 

dot_plot

ggsave(
  "dot_plot.pdf",
  dot_plot,
  width=10,
  height=4
)

# vlnplot ----

p1 <- VlnPlot(bischoff_epi_tumor, features = "SFTPC", pt.size = 0) +
  NoLegend() +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold"),
    text = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_blank(),
  ) +
  labs(x = NULL, y = "SFTPC")

ggsave(
  "vlnplot_1.pdf",
  p1,
  width=4.5,
  height=1.75
)

p2 <- VlnPlot(bischoff_epi_tumor, features = "SCGB3A2", pt.size = 0) + 
  NoLegend() + 
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold"),
    text = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_blank(),
  ) +
  labs(x = NULL, y = "SCGB3A2")

ggsave(
  "vlnplot_2.pdf",
  p2,
  width=4.5,
  height=1.75
)

p3 <- VlnPlot(bischoff_epi_tumor, features = "CEACAM5", pt.size = 0) +
  NoLegend() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold"),
    text = element_text(size = 14),
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_blank(),
  ) +
  labs(x = NULL, y = "CEACAM5")

ggsave(
  "vlnplot_3.pdf",
  p3,
  width=4.5,
  height=3
)

