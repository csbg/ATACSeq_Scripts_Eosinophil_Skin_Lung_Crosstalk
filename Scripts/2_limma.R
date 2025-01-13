library(limma)
library(edgeR)
library(GenomicRanges)
library(tidyverse)
library(enrichR)
library(pheatmap)

# Load Counts
counts <- read.csv("OutATAC-seq_eosinophils/Counts.csv", row.names = 1)
colnames(counts) <- sub(".filtered.bam$", "", colnames(counts))

# meta.table
meta.table <- read.table("meta.table.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t")

# Replace spaces with dots in Disease and Cell_Type columns
meta.table$Disease <- gsub(" ", ".", meta.table$Disease)
meta.table$Cell_Type <- gsub(" ", ".", meta.table$Cell_Type)

## Making sure that the rownames of the meta.table are the same order as the colnames of counts
rownames(meta.table) <- meta.table$sample_name
counts <- counts[, rownames(meta.table)]

# Voom Normalization
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge) # add filter after here to drop reads high logFC and lower counts
threshold <- 10
drop <- which(apply(cpm(dge), 1, max) < threshold)
dge <- dge[-drop, ]
v <- voom(dge, design = NULL, plot = TRUE)

# Save Voom plot as PDF
pdf("./OutATAC-seq_eosinophils/voom_normalization.pdf")
v <- voom(dge, design = NULL, plot = TRUE)
dev.off()

# Design Matrix
# Creating the design matrix with only Disease as a factor
design <- model.matrix(~ Disease, data = meta.table)

# Visualize the design matrix and save to PDF
pdf("./OutATAC-seq_eosinophils/design_matrix_heatmap.pdf")
pheatmap(design)
dev.off()

# Confirm that the column names of the design matrix are as expected
colnames(design)
row.names(design)

# Fit the model with the interaction term
fit <- lmFit(v, design)
fit <- eBayes(fit)

head(coef(fit))

# Read probes data
probes <- read.csv("OutATAC-seq_eosinophils/Probes.csv", stringsAsFactors = FALSE, header = TRUE, sep = "\t")

# Generate topTable for each coefficient
limmaRes <- list()
for (coefx in colnames(coef(fit))) {
  limmaRes[[coefx]] <- topTable(fit, coef = coefx, number = Inf)
}
limmaRes <- bind_rows(limmaRes, .id = "coef")
limmaRes <- filter(limmaRes, coef != "(Intercept)")

#### Plots and Visualization ########################################################

### Volcano Plot 
pdf("OutATAC-seq_eosinophils/volcano_plot.pdf")
ggplot(limmaRes, aes(x = logFC, y = -log10(P.Value), color = adj.P.Val < 0.05)) +
  geom_point() +
  facet_wrap(~coef) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
  labs(x = "Log Fold Change", y = "-Log10 P-value", title = "Volcano Plot of Differential Expression") +
  theme_minimal()
dev.off()

### Histogram of P-Values
pdf("OutATAC-seq_eosinophils/p_value_histogram.pdf")
ggplot(limmaRes, aes(x = P.Value)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  facet_wrap(~coef) +
  labs(x = "P-value", y = "Frequency", title = "Histogram of P-values") +
  theme_minimal()
dev.off()

#### Significant Hits with Gene Names ####
# Convert the row names of limmaRes to a column called 'probe'
limmaRes <- limmaRes %>%
  rownames_to_column(var = "probe")

# Merge limmaRes with probes based on the 'probe' column
limmaRes_with_gene <- merge(limmaRes, probes[, c("probe", "gene")], by = "probe", all.x = TRUE)

# Check if the merge worked and gene names have been added
head(limmaRes_with_gene)

# Filter for significant hits with adj.P.Val < 0.05
significant_hits <- limmaRes_with_gene %>%
  filter(adj.P.Val < 0.05)

# Keep only the 'probe', 'gene', and other relevant columns (logFC, adj.P.Val, etc.)
significant_hits <- significant_hits %>%
  select(probe, gene, logFC, adj.P.Val, P.Value, AveExpr, t, B)

# Check the structure of the significant hits
# Get a logical vector that indicates which rows are complete cases
complete_sig_genes_indices <- complete.cases(significant_hits)

# Use this vector to subset the original dataframe
complete_sig_genes <- significant_hits[complete_sig_genes_indices, ]

# Check the structure of the new dataframe
str(complete_sig_genes)

# Optionally, view the first few rows to confirm
head(complete_sig_genes)

# Save significant hits to CSV
write.csv(significant_hits, "OutATAC-seq_eosinophils/significant_hits_with_gene_names.csv", row.names = FALSE)

#Excel
library(openxlsx)
write.xlsx(complete_sig_genes, file = "OutATAC-seq_eosinophils/complete_sig_genes.xlsx", rowNames = FALSE)
