library(ggrepel)

# Separate significant genes based on logFC
open_chromatin_genes <- significant_hits %>%
  filter(logFC > 0) %>%
  arrange(adj.P.Val) %>%
  head(10)

closed_chromatin_genes <- significant_hits %>%
  filter(logFC < 0) %>%
  arrange(adj.P.Val) %>%
  head(10)

# Combine the top 10 from each group
top_genes <- bind_rows(open_chromatin_genes, closed_chromatin_genes)

# Volcano Plot with Annotations
pdf("OutATAC-seq_eosinophils/volcano_plot_annotated.pdf")
ggplot(limmaRes_with_gene, aes(x = logFC, y = -log10(P.Value), color = adj.P.Val < 0.05)) +
  geom_point() +
  facet_wrap(~coef) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
  labs(x = "Log Fold Change", y = "-Log10 P-value", title = "Volcano Plot of Differential Expression") +
  theme_minimal() +
  # Add labels for the top 10 genes in each group
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.5,
    segment.color = "black"
  )
dev.off()


# Filter out rows where the gene column is NA
mapped_genes_no_na <- limmaRes_with_gene %>%
  filter(!is.na(gene))

# Save the filtered results to an Excel file
library(openxlsx)
write.xlsx(mapped_genes_no_na, file = "OutATAC-seq_eosinophils/mapped_genes_no_na.xlsx", rowNames = FALSE)
