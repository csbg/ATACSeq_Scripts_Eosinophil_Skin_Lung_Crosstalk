library(ComplexHeatmap)
library(circlize)
library(dplyr)

# Assuming 'v$E' contains normalized expression data
scaled_vE <- t(scale(t(v$E)))

# Use the significant_hits object that has the correct gene names mapped
# Select the top 50 significant probes based on adj.P.Val
significant_hits_subset <- significant_hits %>%
  arrange(adj.P.Val) %>%
  head(50)  # Select the top 50 significant probes

# Subset the scaled_vE matrix to include only the significant probes
scaled_vE_subset <- scaled_vE[significant_hits_subset$probe, ]

# Create a mapping of probes to genes using the correct gene names from significant_hits
probe_to_gene <- setNames(significant_hits_subset$gene, significant_hits_subset$probe)

# Update row names in scaled_vE_subset to gene names from the significant_hits object
rownames(scaled_vE_subset) <- probe_to_gene[rownames(scaled_vE_subset)]

# Ensure the order of columns in scaled_vE_subset matches the sample order in meta.table
meta.table$sample_name <- make.names(meta.table$sample_name)
scaled_vE_subset <- scaled_vE_subset[, meta.table$sample_name]

# Define group colors for disease
group_colors <- c(
  "PBS.Control." = "blue",
  "S.aureus.skin.infection." = "red"
)

# Update meta.table to use Disease as the grouping factor
meta.table$Group <- meta.table$Disease

# Create the row annotation for the samples based on the disease group
row_annotation_samples <- rowAnnotation(
  Samples = factor(meta.table$Group, levels = unique(meta.table$Group)),
  col = list(Samples = group_colors),
  show_annotation_name = TRUE
)

# Create and draw the heatmap with gene names as row labels
ht_list <- Heatmap(
  t(scaled_vE_subset),  # Transpose to have genes on rows, samples on columns
  name = "Expression",
  left_annotation = row_annotation_samples,
  show_row_names = TRUE,  # Ensure row names (gene names) are shown
  show_column_names = TRUE,
  col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  show_row_dend = FALSE,
  show_column_dend = FALSE
)

# Draw the heatmap
draw(ht_list, heatmap_legend_side = "bot", annotation_legend_side = "bot")
