library(Seurat)
library(dplyr)

create_xenium_annotation <- function(
    seurat_obj,
    annotation_col = "annotation",
    output_file = "xenium_cell_annotation.csv"
) {
  
  if (!annotation_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste0(
      "Column '", annotation_col, 
      "' not found in Seurat metadata."
    ))
  }
  
  xenium_annotation <- data.frame(
    cell_id = colnames(seurat_obj),
    group = seurat_obj[[annotation_col]][, 1],
    stringsAsFactors = FALSE
  )
  
  celltypes <- unique(xenium_annotation$group)
  
  n <- length(celltypes)
  polychrome_colors <- pals::polychrome(max(n))
  
  color_table <- data.frame(
    group = celltypes,
    color = polychrome_colors[seq_along(celltypes)],
    stringsAsFactors = FALSE
  )
  
  xenium_annotation <- xenium_annotation %>%
    dplyr::left_join(color_table, by = "group")
  
  write.csv(
    xenium_annotation,
    file = output_file,
    row.names = FALSE,
    quote = FALSE
  )
  
  message("Xenium annotation saved to: ", output_file)
  
  return(xenium_annotation)
}

create_xenium_annotation(
  seurat_obj = xenium_2,
  annotation_col = "annotation",
  output_file = "xenium_2_cell_annotation.csv"
)
