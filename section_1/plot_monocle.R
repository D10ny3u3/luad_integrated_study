# monocle type ----

p <- plot_cell_trajectory(cds, color_by = "celltype") +
  # ggtitle("Query MIA") + 
  # labs(
  #   x = "UMAP_1",
  #   y = "UMAP_2"
  # ) +
  theme(
    axis.title = element_text(size = 16, face = "bold"), # axis title bold and large
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size=12, face="bold"),
    legend.text = element_text(size=12),
    legend.position = "right", 
    # plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  # plot title bold and large
  )

p

ggsave(
  "monocle_2_celltype.pdf",
  p,
  width=6,
  height=4
)


# monocle state ----

p <- plot_cell_trajectory(cds, color_by = "State") +
  # ggtitle("Query MIA") + 
  # labs(
  #   x = "UMAP_1",
  #   y = "UMAP_2"
  # ) +
  theme(
    axis.title = element_text(size = 16, face = "bold"), # axis title bold and large
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size=12, face="bold"),
    legend.text = element_text(size=12),
    # plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  # plot title bold and large
  )

p

ggsave(
  "monocle_2_state.pdf",
  p,
  width=4,
  height=4
)


# monocle pseudotime ----

p <- plot_cell_trajectory(cds, color_by = "Pseudotime") +
  # ggtitle("Query MIA") + 
  # labs(
  #   x = "UMAP_1",
  #   y = "UMAP_2"
  # ) +
  theme(
    axis.title = element_text(size = 16, face = "bold"), # axis title bold and large
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size=12, face="bold"),
    legend.text = element_text(size=12),
    # plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  # plot title bold and large
  )

p

ggsave(
  "monocle_2_time.pdf",
  p,
  width=4,
  height=4
)

# pie ----

table(cds$State, cds$celltype)

library(dplyr)
library(tidyr)
library(ggplot2)

df <- as.data.frame(table(cds$State, cds$celltype))
colnames(df) <- c("State", "celltype", "Count")

df <- df %>%
  group_by(State) %>%
  mutate(
    Percent = Count / sum(Count) * 100
  )

library(ggplot2)
library(dplyr)

df <- as.data.frame(table(cds$State, cds$celltype))
colnames(df) <- c("State", "celltype", "Count")

for (s in unique(df$State)) {
  
  df_sub <- df %>%
    filter(State == s)
  
  p <- ggplot(df_sub, aes(x = "", y = Count, fill = celltype)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    theme_void() +
    theme(
      legend.position = "none"
    )
  
  ggsave(
    filename = paste0("State_", s, "_pie.pdf"),
    plot = p,
    width = 5,
    height = 5
  )
}
