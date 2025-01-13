# Load necessary libraries
library(limma)
library(edgeR)
library(GenomicRanges)
library(tidyverse)
library(enrichR)
library(pheatmap)
library(ggplot2)

# Define the datasets for enrichment analysis
enrich_dbs <- c("WikiPathways_2024_Mouse",
                "TRRUST_Transcription_Factors_2019",
                "GO_Biological_Process_2023",
                "GO_Cellular_Component_2023",
                "Reactome_2022",
                "KEGG_2019_Mouse")

# Load data
counts <- read.csv("OutATAC-seq_eosinophils/Counts.csv", row.names = 1)
colnames(counts) <- sub(".filtered.bam$", "", colnames(counts))

# Load metadata
meta.table <- read.table("meta.table.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t")
meta.table$Disease <- gsub(" ", ".", meta.table$Disease)
meta.table$Cell_Type <- gsub(" ", ".", meta.table$Cell_Type)
rownames(meta.table) <- meta.table$sample_name
counts <- counts[, rownames(meta.table)] 

# Voom normalization and design matrix
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)
threshold <- 10
drop <- which(apply(cpm(dge), 1, max) < threshold)
dge <- dge[-drop, ]
v <- voom(dge, design = NULL, plot = TRUE)
design <- model.matrix(~ Disease, data = meta.table)
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Extract the significant results for the coefficient of interest (PBS vs SA)
limmaRes <- topTable(fit, coef = 2, number = Inf) # Adjust coef based on experiment
limmaRes <- filter(limmaRes, adj.P.Val < 0.05) # Include all significant genes

# Load gene-probe mappings
probes <- read.csv("OutATAC-seq_eosinophils/Probes.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
limmaRes <- limmaRes %>%
  rownames_to_column(var = "probe")

limmaRes <- merge(limmaRes, probes[, c("probe", "gene")], by = "probe", all.x = TRUE)

# Divide genes by chromatin state
open_chromatin_genes <- unique(limmaRes$gene[limmaRes$logFC > 0 & !is.na(limmaRes$gene)])  # Open chromatin (upregulated in SA)
closed_chromatin_genes <- unique(limmaRes$gene[limmaRes$logFC < 0 & !is.na(limmaRes$gene)]) # Closed chromatin (downregulated in SA)

# Perform enrichment analysis for each group (open and closed chromatin)
all_enrich_results_open <- list()
all_enrich_results_closed <- list()

# Enrichment for open chromatin genes
if (length(open_chromatin_genes) > 0) {
  for (db in enrich_dbs) {
    enrich_results <- enrichr(open_chromatin_genes, database = db)
    
    if (length(enrich_results) > 0) {
      if ("Adjusted.P.value" %in% colnames(enrich_results[[db]])) {
        all_enrich_results_open[[db]] <- enrich_results[[db]]
        write.csv(enrich_results[[db]], paste0("enrichment_open_", db, ".csv"), row.names = FALSE)
      }
    }
  }
}

# Enrichment for closed chromatin genes
if (length(closed_chromatin_genes) > 0) {
  for (db in enrich_dbs) {
    enrich_results <- enrichr(closed_chromatin_genes, database = db)
    
    if (length(enrich_results) > 0) {
      if ("Adjusted.P.value" %in% colnames(enrich_results[[db]])) {
        all_enrich_results_closed[[db]] <- enrich_results[[db]]
        write.csv(enrich_results[[db]], paste0("enrichment_closed_", db, ".csv"), row.names = FALSE)
      }
    }
  }
}

# Combine and filter results for top pathways (open and closed chromatin)
combined_enrichment <- bind_rows(
  imap_dfr(all_enrich_results_open, ~ .x %>%
             filter(Adjusted.P.value < 0.05) %>%
             mutate(chromatin_state = "Open", Database = .y)),
  imap_dfr(all_enrich_results_closed, ~ .x %>%
             filter(Adjusted.P.value < 0.05) %>%
             mutate(chromatin_state = "Closed", Database = .y))
)

# Prepare data for heatmap plot
plot_data <- combined_enrichment %>%
  mutate(
    p_size = -log10(Adjusted.P.value),
    logFC = log2(Odds.Ratio),  # Use log2 of the Odds Ratio as the log fold change value
    Pathway = Term
  )

# Save plot_data as CSV
write.csv(plot_data, "OutATAC-seq_eosinophils/plot_data_combined_enrichment_results.csv", row.names = FALSE)

# Create heatmap with ggplot
p <- ggplot(plot_data, aes(x = chromatin_state, y = Pathway, size = p_size, fill = logFC)) +
  geom_point(shape = 21, color = "black") +
  scale_size_continuous(range = c(3, 12), name = "-log10(Adjusted P-value)") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Log2(Odds Ratio)") +
  facet_grid(Database ~ ., scales = "free_y", space = "free_y") +  # Facet by database
  theme_bw(base_size = 12) +
  theme(strip.text.y = element_text(angle = 0), # Make facet labels horizontal
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        panel.spacing.y = unit(0.8, "lines"),
        legend.position = "right") +
  labs(x = 'Chromatin State', y = 'Pathway', fill = "Log2(Odds Ratio)") +
  ggtitle("Top Enriched Pathways across Open and Closed Chromatin States")

# Save the plot as PDF
ggsave("OutATAC-seq_eosinophils/combined_enrichment_heatmap.pdf", plot = p, width = 12, height = 18)

# Display the plot
print(p)
