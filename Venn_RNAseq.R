if (!require(devtools)) install.packages("devtools")
if (!require(devtools)) devtools::install_github("yanlinlin82/ggvenn")

library(ggvenn)

### set working directory
setwd('/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/figures')

##### RNAseq
jag_RNA = read.csv('/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/DEG_significant/jag_DEG.txt_sig_genes_BSU_Numbers.txt', sep = '\t', header = FALSE)
khpA_RNA = read.csv('/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/DEG_significant/khpA_DEG.txt_sig_genes_BSU_Numbers.txt', sep = '\t', header = FALSE)
kre_RNA = read.csv('/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/DEG_significant/kre_DEG.txt_sig_genes_BSU_Numbers.txt', sep = '\t', header = FALSE)

x <- list(jag=jag_RNA$V1, khpA=khpA_RNA$V1, kre=kre_RNA$V1)

ggp <- ggvenn(
  x,
  fill_color = c('#E31A1C', '#1F78B4', "#33A02C"),
  stroke_size = 0.5, set_name_size = 5, show_percentage = FALSE
)
ggp

# Specify the output directory and file name
output_dir <- "/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpA_kre/X204SC22071444-Z01-F001/Analysis/figures/"
output_file <- file.path(output_dir, "Venn_RNAseq.pdf")
# Save the plot as a PDF file
ggsave(output_file, plot = ggp, width = 8, height = 6)


### identify overlapping genes
overlap_All <- Reduce(intersect, list(jag_df$V1, khpA_df$V1, kre_df$V1))
overlap_All

overlap_jag_khpA <- setdiff(Reduce(intersect, list(jag_df$V1, khpA_df$V1)), overlap_All)
overlap_jag_khpA

overlap_jag_kre <- setdiff(Reduce(intersect, list(jag_df$V1, kre_df$V1)), overlap_All)
overlap_jag_kre

overlap_khpA_kre <- setdiff(Reduce(intersect, list(khpA_df$V1, kre_df$V1)), overlap_All)
overlap_khpA_kre



#### CRAC
jag_CRAC = read.csv('/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/ratios/prepared/jag_ratios.txt_sig_genes_BSU_Numbers.txt', sep = '\t', header = FALSE)
khpA_CRAC = read.csv('/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/ratios/prepared/khpA_ratios.txt_sig_genes_BSU_Numbers.txt', sep = '\t', header = FALSE)
kre_CRAC = read.csv('/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/ratios/prepared/kre_ratios.txt_sig_genes_BSU_Numbers.txt', sep = '\t', header = FALSE)

x <- list(jag=jag_CRAC$V1, khpA=khpA_CRAC$V1[1:200], kre=kre_CRAC$V1)

ggp <- ggvenn(
  x,
  fill_color = c('#E31A1C', '#1F78B4', "#33A02C"),
  stroke_size = 0.5, set_name_size = 5, show_percentage = FALSE
)
ggp

# Specify the output directory and file name
output_dir <- "/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/figures/"
output_file <- file.path(output_dir, "Venn_CRAC.pdf")
# Save the plot as a PDF file
ggsave(output_file, plot = ggp, width = 8, height = 6)


### identify overlapping genes
overlap_All <- Reduce(intersect, list(jag_df$V1, khpA_df$V1, kre_df$V1))
overlap_All

overlap_jag_khpA <- Reduce(intersect, list(jag_CRAC$V1, khpA_CRAC$V1))
overlap_jag_khpA

overlap_jag_kre <- setdiff(Reduce(intersect, list(jag_df$V1, kre_df$V1)), overlap_All)
overlap_jag_kre

overlap_khpA_kre <- setdiff(Reduce(intersect, list(khpA_df$V1, kre_df$V1)), overlap_All)
overlap_khpA_kre


### overlap between RNAseq and CRAC
## jag
x <- list(jag_RNA=jag_RNA$V1, jag_CRAC=jag_CRAC$V1)

ggp <- ggvenn(
  x,
  fill_color = c('#E31A1C', '#1F78B4', "#33A02C"),
  stroke_size = 0.5, set_name_size = 5, show_percentage = FALSE
)
ggp

overlap_jag <- Reduce(intersect, list(jag_RNA$V1, jag_CRAC$V1))
overlap_jag

# Specify the output directory and file name
output_dir <- "/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_30thMay2024/figures/"
output_file <- file.path(output_dir, "Venn_jag.pdf")
# Save the plot as a PDF file
ggsave(output_file, plot = ggp, width = 8, height = 6)

## khpA
x <- list(khpA_RNA=khpA_RNA$V1, khpA_CRAC=khpA_CRAC$V1)

ggp <- ggvenn(
  x,
  fill_color = c('#E31A1C', '#1F78B4', "#33A02C"),
  stroke_size = 0.5, set_name_size = 5, show_percentage = FALSE
)
ggp

overlap_khpA <- Reduce(intersect, list(khpA_RNA$V1, khpA_CRAC$V1))
overlap_khpA

# Specify the output directory and file name
output_dir <- "/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_30thMay2024/figures/"
output_file <- file.path(output_dir, "Venn_khpA.pdf")
# Save the plot as a PDF file
ggsave(output_file, plot = ggp, width = 8, height = 6)

## kre
x <- list(kre_RNA=kre_RNA$V1, kre_CRAC=kre_CRAC$V1)

ggp <- ggvenn(
  x,
  fill_color = c('#E31A1C', '#1F78B4', "#33A02C"),
  stroke_size = 0.5, set_name_size = 5, show_percentage = FALSE
)
ggp

overlap_kre <- Reduce(intersect, list(kre_RNA$V1, kre_CRAC$V1))
overlap_kre

# Specify the output directory and file name
output_dir <- "/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_30thMay2024/figures/"
output_file <- file.path(output_dir, "Venn_kre.pdf")
# Save the plot as a PDF file
ggsave(output_file, plot = ggp, width = 8, height = 6)
