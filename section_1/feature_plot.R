markers <- c(
  "MUC5B",
  "MUC5AC",
  "SCGB1A1",
  "SCGB3A2",
  "SFTPC",
  "nFeature_RNA"
)

FeaturePlot(soil_club, c(
  "SCGB1A1", "SCGB3A2"
))

markers <- c(
  "MUC5B",
  "MUC5AC",
  "SCGB1A1",
  "SCGB3A2",
  "SFTPC"
)

p <- FeaturePlot(
  soil_club,
  features = markers,
  ncol = 3,
  order = TRUE,
  pt.size = 0.2
) & 
  labs(
    x = "UMAP_1",
    y = "UMAP_2"
  ) &
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
    # axis.line = element_blank()
  )

p

ggsave(
  "feat_club.pdf",
  p,
  width = 14,
  height = 8
)
