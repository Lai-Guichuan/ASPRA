## Installation Process for the ASPRA Package

After the article is accepted, the ASPRA package will be uploaded. To install the ASPRA package, use the following steps:

```R
library(devtools)
devtools::install_github("Lai-Guichuan/ASPRA")

Workflow for Calculating PBIS-TP53 Scores in Single-Cell Data Using the ASPRA Package:
counts <- GetAssayData(object = scobj, slot = "counts")
gene_lengths <- read.table("gencode.v22.annotation.txt", header = TRUE, sep = "\t")
genes_in_counts <- rownames(counts)
matched_gene_lengths <- gene_lengths[gene_lengths$gene %in% genes_in_counts, ]
counts <- counts[matched_gene_lengths$gene, ]
gene_lengths_vector <- matched_gene_lengths$length
RPK <- counts / (gene_lengths_vector / 1000)
TPM <- sweep(RPK, 2, colSums(RPK), FUN = "/") * 1e6

library(ASPRA)
filtered_gene_sets <- readRDS("filtered_gene_sets.rds")
coef_file <- readRDS("coef_file.rds")
pos_file <- readRDS("pos_file.rds")
neg_file <- readRDS("neg_file.rds")
expr_file <- "TPM.txt"
result <- calculate_total_score(expr_file, filtered_gene_sets, coef_file, pos_file, neg_file)
write.table(result, file = "PBIS-TP53.txt", sep = "\t", quote = FALSE)
