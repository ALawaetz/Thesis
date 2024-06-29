library('DESeq2')
library('ggplot2')
library( "gplots" )
library( "RColorBrewer" )
library('genefilter')
library('pheatmap')
library('FactoMinheR')
library(ggvenn)
library(venneuler)
library(ggrepel)
library(tibble)

directory <- '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpA_kre/X204SC22071444-Z01-F001/Analysis/Coverage_bed'
sampleFiles <- list.files(path=directory, pattern='WT|kre')
# Filter out files that match the exclude pattern
sampleFiles <- sampleFiles[!grepl('ALL_merged', sampleFiles)]
sampleFiles
sampleNames <- c('kre_1', 'kre_2', 'kre_3', 'WT_1', 'WT_2', 'WT_3')
OD = c(1, 1, 1, 2, 2, 1)
Names <- c('kre_1', 'kre_2', 'kre_3', 'WT_1', 'WT_2', 'WT_3')
Strain <- c('kre', 'kre', 'kre', 'WT', 'WT', 'WT')
replicate <- c(1, 2, 3, 1, 2, 3)

sampleTable<-data.frame(sampleName=sampleNames, fileName=sampleFiles, names = Names, strain=Strain, OD = OD, replicate = replicate)

sampleTable$OD <- as.factor(sampleTable$OD)
sampleTable$strain <- as.factor(sampleTable$strain)
sampleTable

dds<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~ OD + strain)

colData(dds)$strain<-factor(colData(dds)$strain, levels=c('WT','kre'))

dds <- DESeq(dds)

### PCA plot
rld <- rlogTransformation(dds, blind=TRUE) 
pca_plot <- plotPCA(rld, intgroup = c("strain", 'replicate'), ntop=5000)
pca_data <- pca_plot$data

pca_plot <- ggplot(data = pca_data, aes(x = PC1, y = PC2, color = strain, shape = factor(replicate))) +
  geom_point(size = 3) +
  labs(title = "kre") +
  xlab(pca_plot$labels$x) +
  ylab(pca_plot$labels$y) +
  labs(colour = 'strain', shape = 'replicate') +
  theme_minimal() + 
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  ) + 
  guides(colour = guide_legend(order = 1))

pca_plot

# Specify the output directory and file name
output_dir <- "/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpA_kre/X204SC22071444-Z01-F001/Analysis/figures/"
output_file <- file.path(output_dir, "PCA_plot_kre.pdf")
# Save the plot as a PDF file
ggsave(output_file, plot = pca_plot, width = 8, height = 6)


### kre
### compute fold changes and statistics
res <- results(dds, contrast = c("strain", "kre", "WT")) 
### remove rows with padj NA values
res <-subset(res, (!is.na(res[,6])))
# Extract the base gene names
res$Gene_base <- sub("\\..*$", "", rownames(res))
# Remove duplicate rows based on the Gene_base column and keep the first occurrence
res <- res[!duplicated(res$Gene_base), ]
### omit rows
rows_to_remove <- c("kre")
# Use negative indexing to exclude rows by their names
res <- subset(res, Gene_base != rows_to_remove)
res <- as.data.frame(res)
volcano_plot <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = (padj <= 0.05) & (abs(log2FoldChange) >= 0.5)), size = 2) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(
    title = "kre",
    x = "Log2(Fold Change)",
    y = "-log10(P-Value)"
  ) +
  guides(color = 'none') +
  geom_text_repel(
    aes(label = ifelse(padj <= 0.05 & abs(log2FoldChange) >= 0.5, Gene_base, "")),
    size = 3, nudge_y = 0.1  # Adjust label size and position as needed
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(c(-max(abs(res$log2FoldChange)), max(abs(res$log2FoldChange)))) 
volcano_plot
# Specify the output directory and file name
output_dir <- "/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpA_kre/X204SC22071444-Z01-F001/Analysis/figures/"
output_file <- file.path(output_dir, "Volcano_plot_kre.pdf")
# Save the plot as a PDF file
ggsave(output_file, plot = volcano_plot, width = 8, height = 6)


### Save DESeq2 table
### compute fold changes and statistics
res <- results(dds, contrast = c("strain", "kre", "WT")) 
res$sig <- ifelse(res$padj <= 0.05, 1, 0)
#Calculate the order based on 'sig' and absolute 'log2FoldChange'
order_indices <- with(res, order(-sig, -abs(log2FoldChange)))
# Reorder the DESeqResults object
res <- res[order_indices, ]

# Extract the base gene names
res$Gene_base <- sub("\\..*$", "", rownames(res))
### omit rows
rows_to_remove <- c("kre")
# Use negative indexing to exclude rows by their names
res <- subset(res, Gene_base != rows_to_remove)
res <- as.data.frame(res)
res <- rownames_to_column(res, var = "Gene")
# save table
# Specify the output directory and file name
output_dir <- "/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpA_kre/X204SC22071444-Z01-F001/Analysis/DEG/"
output_file <- file.path(output_dir, "kre_DEG.txt")
write.table(res, file = output_file, sep = '\t', row.names = FALSE, quote = FALSE)


### At this point, run python script "Prepare_DEG_for_GO_analysis.ipynb"

### Then do GO term analysis and KEGG
### Make custom GO terms from SubtiWiki to input ClusterProfiler
# Load necessary libraries
library(clusterProfiler)
library(dplyr)

# Load your custom gene-to-category data
gene2category <- read.csv("/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpA_kre/X204SC22071444-Z01-F001/Analysis/Gene_categories/GeneID_Category_All.csv")
colnames(gene2category) <- c("Category", "GeneID")

# Create a TERM2GENE data frame
TERM2GENE <- gene2category %>%
  rename(term = Category, gene = GeneID)

# Create a TERM2NAME data frame (optional, but useful for readability)
TERM2NAME <- gene2category %>%
  select(Category) %>%
  distinct() %>%
  mutate(name = Category) %>%
  rename(term = Category)

# Specify the path to your gene list file
### RNAseq
gene_list_file <- "/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpA_kre/X204SC22071444-Z01-F001/Analysis/DEG_significant/kre_DEG.txt_sig_genes_BSU_Numbers.txt"
# Read the gene list into R
genes_of_interest <- read.table(gene_list_file, header = FALSE, stringsAsFactors = FALSE)
genes_of_interest <- as.character(genes_of_interest$V1)

### CRAC
# Specify the path to your gene list file
gene_list_file <- "/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/ratios/prepared/kre_ratios.txt_sig_genes_BSU_Numbers.txt"
# Read the gene list into R
genes_of_interest <- read.table(gene_list_file, header = FALSE, stringsAsFactors = FALSE)
genes_of_interest <- as.character(genes_of_interest$V1)
genes_of_interest <- unique(genes_of_interest)
genes_of_interest <- genes_of_interest[1:20]


# Perform GO-like enrichment analysis
go_enrichment <- enricher(genes_of_interest, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME)

# View results
print(go_enrichment)
# go_enrichment[go_enrichment$ID == "Ribosomal proteins", ]
# View(go_enrichment)
# 
# # Access individual gene sets
# gene_set_name <- "Biosynthesis of antibacterial compounds"
# if (gene_set_name %in% names(go_enrichment@geneSets)) {
#   individual_gene_set <- go_enrichment@geneSets[[gene_set_name]]
#   print(individual_gene_set)
# } else {
#   cat("Gene set not found:", gene_set_name, "\n")
# }

# Dot plot of the enrichment results
my_dotplot <- dotplot(go_enrichment, showCategory = 20) + 
  labs(title = "Kre") +
  theme(plot.title = element_text(hjust = 0.5))
my_dotplot

# Specify the output directory and file name
#output_dir <- "/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/figures/"
output_dir <- "/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/figures"
#output_file <- file.path(output_dir, "Dotplot_GSEA_jag_RNAseq.pdf")
output_file <- file.path(output_dir, "Dotplot_GSEA_kre_CRAC.pdf")
# Save the plot as a PDF file
ggsave(output_file, plot = my_dotplot, width = 8, height = 6)




# Load required libraries
library(org.Bsubtilis.eg.db)  # Replace 'Xx' with the appropriate organism code
# Specify the path to your gene list file
gene_list_file <- "/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpA_kre/X204SC22071444-Z01-F001/Analysis/DEG_significant/kre_DEG.txt_sig_genes_BSUNumbers_no_underscore.txt"
# Read the gene list into R
genes_of_interest <- read.table(gene_list_file, header = FALSE, stringsAsFactors = FALSE)
genes_of_interest <- as.character(genes_of_interest$V1)

# Perform KEGG pathway enrichment analysis
kegg_enrichment <- enrichKEGG(gene = gene_listBSU,
                              organism = "bsu",  # Replace 'Xx' with the appropriate organism code
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05)

# View the results
head(kegg_enrichment)
# Create a dot plot
keggplot <- dotplot(kegg_enrichment) + ggtitle("KEGG Enrichment Analysis")

# Specify the output directory and file name
output_dir <- "/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpA_kre/X204SC22071444-Z01-F001/Analysis/figures/"
output_file <- file.path(output_dir, "Dotplot_KEGG_kre.pdf")
# Save the plot as a PDF file
ggsave(output_file, plot = keggplot, width = 8, height = 6)

