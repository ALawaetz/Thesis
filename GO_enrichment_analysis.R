#install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
library(clusterProfiler)

# Load BiocManager and install AnnotationHub
#BiocManager::install("AnnotationHub")
# Load AnnotationHub
library(AnnotationHub)

# Install biomaRt using BiocManager
#BiocManager::install("biomaRt")

#BiocManager::install("AnnotationForge")
library(AnnotationForge)
### Make custom bacillus subtilis orgDb project
makeOrgPackageFromNCBI(version = "0.1",
                       author = "AC Lawaetz",
                       maintainer = "AC Lawaetz <acl58@bath.ac.uk",
                       outputDir = "/Users/andreas/Bacillus/Bioinformatics/ClusterProfiler_orgDB",
                       tax_id = "224308",
                       genus = "Bacillus",
                       species = "subtilis")


### Make custom GO terms from SubtiWiki to input ClusterProfiler
# Load necessary libraries
library(clusterProfiler)
library(dplyr)

# Load your custom gene-to-category data
gene2category <- read.csv("/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/Gene_categories/GeneID_Category_All.csv")
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
gene_list_file <- "/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpA/sig_DEG/pyRAP_Edinburgh_jag_all.csv_sig_genes_BSU_Numbers.txt"
# Read the gene list into R
genes_of_interest <- read.table(gene_list_file, header = FALSE, stringsAsFactors = FALSE)
genes_of_interest <- as.character(genes_of_interest$V1)

# Perform GO-like enrichment analysis
go_enrichment <- enricher(genes_of_interest, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME)

# View results
print(go_enrichment)

# Dot plot of the enrichment results
dotplot(go_enrichment, showCategory = 20) + 
  labs(title = "GO Enrichment Analysis Jag") +
  theme(plot.title = element_text(hjust = 0.5))



# Print the current cache directory
current_cache <- getAnnotationHubOption("CACHE")
print(current_cache)

# Load your custom OrgDb package
install.packages("/Users/andreas/Bacillus/Bioinformatics/ClusterProfiler_orgDB/org.Bsubtilis.eg.db", repos = NULL, type = "source")

library(clusterProfiler)
library(org.Bsubtilis.eg.db)
library(ggplot2)

# Specify the path to your gene list file
gene_list_file <- "/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpA/sig_DEG/pyRAP_Edinburgh_jag_all.csv_sig_genes.txt"
gene_list_BSUfile <- "/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpA/sig_DEG/pyRAP_Edinburgh_jag_all.csv_sig_genes_BSUNumbers_no_underscore.txt"

# Read the gene list into R
gene_list <- read.table(gene_list_file, header = FALSE, stringsAsFactors = FALSE)
gene_listBSU <- read.table(gene_list_BSUfile, header = FALSE, stringsAsFactors = FALSE)


# Convert to a character vector
gene_list <- gene_list$V1
gene_listBSU <- gene_listBSU$V1

# Check supported key types
#keytypes(org.Bsubtilis.eg.db)


# Perform GO enrichment analysis
ego <- enrichGO(gene          = gene_list,
                OrgDb         = org.Bsubtilis.eg.db,
                keyType       = "SYMBOL",  # Adjust the keyType based on your gene identifiers
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)

# View the results
head(ego)

# Create a dot plot
dotplot(ego) + ggtitle("GO Enrichment Analysis")

# Load required libraries
library(clusterProfiler)
library(org.Bsubtilis.eg.db)  # Replace 'Xx' with the appropriate organism code

# Perform KEGG pathway enrichment analysis
kegg_enrichment <- enrichKEGG(gene = gene_listBSU,
                              organism = "bsu",  # Replace 'Xx' with the appropriate organism code
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05)

# View the results
head(kegg_enrichment)
# Create a dot plot
dotplot(kegg_enrichment) + ggtitle("KEGG Enrichment Analysis")



# Extract genes used in the enrichment analysis
included_genes <- kegg_enrichment@result$geneID

# Convert the included genes from a single string to a vector
included_genes_vector <- unlist(strsplit(included_genes, "/"))

# Identify the excluded genes
excluded_genes <- setdiff(gene_list, included_genes_vector)

# Display excluded genes
cat("Excluded genes:\n")
print(excluded_genes)

# Test if "dps" is in the excluded genes list
is_dps_excluded <- "dps" %in% excluded_genes
# Print the result
if (is_dps_excluded) {
  cat("The gene 'dps' is in the excluded genes list.\n")
} else {
  cat("The gene 'dps' is not in the excluded genes list.\n")
}



### Examine orgDb object
library(AnnotationDbi)
# 1. List available keys
available_keys <- keytypes(org.Bsubtilis.eg.db)
print(available_keys)


# 2. Extract gene annotations
# For example, you can extract gene symbols and names
genes <- AnnotationDbi::select(org.Bsubtilis.eg.db, keys = keys(org.Bsubtilis.eg.db), columns = c("SYMBOL", "GENENAME", 'GO', 'ONTOLOGY'), keytype = "ENTREZID")

# Print the first few rows of the gene annotations
print(head(genes))

# Save the gene annotations to a CSV file
write.csv(genes, file = "/Users/andreas/Bacillus/Bioinformatics/RNAseq_spoVG_LB_OD3/Analysis/gene_annotations_orgDb.csv", row.names = FALSE)


### Investigate why the ont='ALL' doesn't include 'dps' when it is included in ont='MF'
#BiocManager::install("GO.db")
# Assuming "dps" is a symbol, get the ENTREZID for it
dps_entrez <- AnnotationDbi::select(org.Bsubtilis.eg.db, keys = "dps", keytype = "SYMBOL", columns = "ENTREZID")

# Get GO annotations for "dps"
dps_go_annotations <- AnnotationDbi::select(org.Bsubtilis.eg.db, keys = dps_entrez$ENTREZID, columns = c("GO", "ONTOLOGY"), keytype = "ENTREZID")

print(dps_go_annotations)

# Perform enrichment analysis for 'MF'
ego_mf <- enrichGO(gene = gene_list,
                   OrgDb = org.Bsubtilis.eg.db,
                   keyType = "SYMBOL",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

# Perform enrichment analysis for 'BP'
ego_bp <- enrichGO(gene = gene_list,
                   OrgDb = org.Bsubtilis.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

# Perform enrichment analysis for 'ALL'
ego_all <- enrichGO(gene = gene_list,
                    OrgDb = org.Bsubtilis.eg.db,
                    keyType = "SYMBOL",
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)

# Compare results
print(head(ego_mf))
print(head(ego_bp))
print(head(ego_all))

# Check if "dps" is included in the results
included_genes_mf <- unique(unlist(strsplit(ego_mf@result$geneID, "/")))
included_genes_bp <- unique(unlist(strsplit(ego_bp@result$geneID, "/")))
included_genes_all <- unique(unlist(strsplit(ego_all@result$geneID, "/")))

"dps" %in% included_genes_mf
"dps" %in% included_genes_bp
"dps" %in% included_genes_all

