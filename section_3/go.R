# go alv ----

alv_genes <- subset(
  target_deg,
  cluster == "Alv_Fib" &
    p_val_adj < 0.05 &
    avg_log2FC > 1 &
    pct.1 > 0.1 &
    (pct.1 - pct.2) > 0.1
)

alv_genes <- alv_genes$gene

library(org.Hs.eg.db)
library(clusterProfiler)

alv_go <- enrichGO(
  gene = alv_genes,  
  OrgDb = org.Hs.eg.db,
  keyType = 'SYMBOL',  
  ont = 'BP',
  pAdjustMethod = 'fdr', 
  pvalueCutoff = 0.05,  
  qvalueCutoff = 0.2,  
  universe = rownames(target), 
  readable = FALSE
)

saveRDS(alv_go, "../alv_go.rds")

p <- dotplot(
  alv_go,
  # showCategory = 20,
  font.size = 8
) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

p

ggsave(
  "alv_go.pdf",
  p,
  width=8,
  height=7
)
