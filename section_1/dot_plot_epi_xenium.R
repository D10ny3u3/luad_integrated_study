
# annotation demo 1 ----

table(xenium_14$annotation)

xenium_14_epi <- subset(xenium_14, annotation %in% c(
  "club_2_like", "epi", "club_1_like", "AT2_like" 
))

DotPlot(
  xenium_14_epi, 
  c("CDH1", "NKX2-1", "MET", "TNC", "SCGB3A2", "SCGB1A1", "SFTPC"))

features = c("CDH1", "NKX2-1", "MET", "TNC", "SCGB3A2", "SCGB1A1", "SFTPC")

xenium_14_epi$annotation <- factor(xenium_14_epi$annotation)
levels(xenium_14_epi$annotation)
new_order <- c("club_1_like", "club_2_like", "AT2_like", "epi")

xenium_14_epi$annotation <- factor(
  xenium_14_epi$annotation,
  levels = new_order
)

dot_plot <- DotPlot(
  xenium_14_epi,
  features = features, 
  # assay = "RNA",
  group.by = "annotation"
) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 14),
    axis.text = element_text(size = 16, face = "bold")
  ) +
  labs(y = NULL, x = NULL)

dot_plot

ggsave(
  "dot_plot_xenium_14_epi.pdf",
  dot_plot,
  width=6,
  height=4
)

# annotation demo 2 ----

table(xenium_tsu_21$annotation)
xenium_tsu_21_epi <- subset(xenium_tsu_21, annotation %in% c(
  "epi", "club_1_like", "club_2_like", "AT2_like" 
))

xenium_tsu_21_epi$annotation <- factor(xenium_tsu_21_epi$annotation)
levels(xenium_tsu_21_epi$annotation)
new_order <- c("club_1_like", "club_2_like", "AT2_like", "epi")

xenium_tsu_21_epi$annotation <- factor(
  xenium_tsu_21_epi$annotation,
  levels = new_order
)

DotPlot(
  xenium_tsu_21_epi, 
  c("CDH1", "NKX2-1", "SCGB3A2", "SCGB1A1", "SFTPC"))

features = c("CDH1", "NKX2-1", "SCGB3A2", "SCGB1A1", "SFTPC")

dot_plot <- DotPlot(
  xenium_tsu_21_epi,
  features = features, 
  # assay = "RNA",
  group.by = "annotation"
) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 14),
    axis.text = element_text(size = 16, face = "bold")
  ) +
  labs(y = NULL, x = NULL)

dot_plot

ggsave(
  "dot_plot_xenium_tsu_21_epi.pdf",
  dot_plot,
  width=6,
  height=4
)
