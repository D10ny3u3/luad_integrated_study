## spink1 sample_2c ----

p <- SpatialFeaturePlot(
  sample_2c, features = c("SPINK1"), 
  images = "slice1_FFPE_2C",
  pt.size.factor = 800) &
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 22)
  )

ggsave(
  "sample_2c_spink.pdf",
  p,
  width=8,
  height=8
)

## spink1 sample 16 ----

p <- SpatialFeaturePlot(
  sample_16, features = c("SPINK1"), 
  images = "slice1",
  pt.size.factor = 500) &
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 22)
  )

ggsave(
  "sample_16_spink.pdf",
  p,
  width=8,
  height=8
)

## ceacam5 sample 5 ----

p <- SpatialFeaturePlot(
  sample_5, features = c("CEACAM5"), 
  images = "slice1",
  pt.size.factor = 500) &
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 22)
  )

ggsave(
  "sample_5_cea.pdf",
  p,
  width=8,
  height=8
)

## ceacam5 sample 17 ----

p <- SpatialFeaturePlot(
  sample_17, features = c("CEACAM5"), 
  images = "slice1",
  pt.size.factor = 500) &
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 22)
  )

ggsave(
  "sample_17_cea.pdf",
  p,
  width=8,
  height=8
)
