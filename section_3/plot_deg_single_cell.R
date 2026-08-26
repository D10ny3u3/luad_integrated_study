Idents(target) <- target$cell_type
target_deg <- FindAllMarkers(
  target, 
  only.pos = T,
  verbose = T)

target_deg$difference <- target_deg$pct.1 - target_deg$pct.2
target_deg$difference <- ifelse(target_deg$cluster == "Adv_Fib", -(target_deg$difference), target_deg$difference)
target_deg$avg_log2FC <- ifelse(target_deg$cluster == "Adv_Fib", -(target_deg$avg_log2FC), target_deg$avg_log2FC)

deg_sig <- target_deg[which(target_deg$p_val_adj<0.05 & abs(target_deg$avg_log2FC) >0.25),]
deg_sig$label <- rownames(deg_sig)

library(ggplot2)
library(ggrepel)

top_genes <- target_deg %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 25) %>%
  ungroup()

str(top_genes)

deg_fib <- ggplot(deg_sig, aes(x=difference, y=avg_log2FC)) + 
  geom_point(size=2, color="grey60") + 
  geom_point(data=top_genes[which(top_genes$p_val_adj<0.05 & top_genes$avg_log2FC>0.1),],
             aes(x=difference, y=avg_log2FC),
             size=2, color="#f0988c")+
  geom_point(data=top_genes[which(top_genes$p_val_adj<0.05 & top_genes$avg_log2FC< -0.1),],
             aes(x=difference, y=avg_log2FC),
             size=2, color="#b883d3")+
  geom_text_repel(data=top_genes, aes(label=gene),
                  color="black",fontface="italic", size=3)+
  theme_classic()+
  theme(axis.text.x = element_text(colour = 'black',size = 12),
        axis.text.y = element_text(colour = 'black',size = 12),
        axis.title = element_text(colour = 'black',size = 15),
        axis.line = element_line(color = 'black', size = 1))+
  geom_hline(yintercept = 0,lty=2,lwd = 1)+
  geom_vline(xintercept = 0,lty=2,lwd = 1)+
  ylab("Log-fold Change")+
  xlab("Delta Percent")

deg_fib

ggsave(
  "deg_fib.pdf",
  deg_fib,
  width=5,
  height=4
)
