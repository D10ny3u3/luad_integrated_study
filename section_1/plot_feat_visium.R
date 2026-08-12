# Visium 14 ----

p <- SpatialFeaturePlot(
  sample_14, features = c("SCGB3A2", "SFTPC", "MARCKSL1", "CRABP2"), 
  images = "slice1",
  pt.size.factor = 800) &
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 22),     # legend标题
    plot.title = element_text(size = 26),       # feature名字
    axis.text = element_text(size = 20),        # 坐标轴文字
    axis.title = element_text(size = 22)
  ) 

ggsave(
  "sample_14_feat.pdf",
  p,
  width=15,
  height=15
)

# Visium 2c ----

p <- SpatialFeaturePlot(sample_2c, features = c("SCGB3A2", "SFTPC", "MARCKSL1", "MUC5B"), 
                   images = "slice1_FFPE_2C",
                   pt.size.factor = 1000) &
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 22),     # legend标题
    plot.title = element_text(size = 26),       # feature名字
    axis.text = element_text(size = 20),        # 坐标轴文字
    axis.title = element_text(size = 22)
  ) 

ggsave(
  "sample_2c_feat.pdf",
  p,
  width=15,
  height=15
)
