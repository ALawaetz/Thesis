if (!require(devtools)) install.packages("devtools")
if (!require(devtools)) devtools::install_github("yanlinlin82/ggvenn")

library(ggvenn)

### set working directory
setwd('/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/figures')

##### RNAseq
jag_RNA = read.csv('/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/DEG_significant/jag_DEG.txt_sig_genes_BSU_Numbers.txt', sep = '\t', header = FALSE)
khpA_RNA = read.csv('/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/DEG_significant/khpA_DEG.txt_sig_genes_BSU_Numbers.txt', sep = '\t', header = FALSE)
kre_RNA = read.csv('/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/DEG_significant/kre_DEG.txt_sig_genes_BSU_Numbers.txt', sep = '\t', header = FALSE)
#kre_feedback = read.csv('/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/DEG_significant/kre_feedback_loop_sig2_genes.tsv', sep = '\t', header = FALSE)
kre_feedback = read.csv('/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/DEG_significant/kre_feedback_loop_sig1_genes.tsv', sep = '\t', header = FALSE)
#kre_feedback = read.csv('/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/DEG_significant/kre_feedback_loop_sig_less_than0_genes.tsv', sep = '\t', header = FALSE)
#kre_feedback = read.csv('/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/DEG_significant/kre_feedback_loop_sig_more_than0_genes.tsv', sep = '\t', header = FALSE)


x <- list(jag=jag_RNA$V1, khpA=khpA_RNA$V1, kre=kre_RNA$V1)
x <- list(kre=kre_RNA$V1, kre_feedback=kre_feedback$V1)

ggp <- ggvenn(
  x,
  fill_color = c('#E31A1C', '#1F78B4', "#33A02C"),
  stroke_size = 0.5, set_name_size = 8, text_size = 8, show_percentage = FALSE
)
ggp

# Specify the output directory and file name
output_dir <- "/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpA_kre/X204SC22071444-Z01-F001/Analysis/figures/"
output_file <- file.path(output_dir, "Venn_RNAseq_foldchange_1.pdf")
# Save the plot as a PDF file
ggsave(output_file, plot = ggp, width = 8, height = 6)


### identify overlapping genes
overlap_kre <- Reduce(intersect, list(kre_RNA$V1, kre_feedback$V1))
overlap_kre

overlap_All <- Reduce(intersect, list(jag_RNA$V1, khpA_RNA$V1, kre_RNA$V1))
overlap_All

overlap_jag_khpA <- setdiff(Reduce(intersect, list(jag_RNA$V1, khpA_RNA$V1)), overlap_All)
overlap_jag_khpA

overlap_jag_kre <- setdiff(Reduce(intersect, list(jag_RNA$V1, kre_RNA$V1)), overlap_All)
overlap_jag_kre

overlap_khpA_kre <- setdiff(Reduce(intersect, list(khpA_RNA$V1, kre_RNA$V1)), overlap_All)
overlap_khpA_kre


# Save the vector to a .txt file
### set working directory
setwd('/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/figures')
writeLines(overlap_jag_khpA, "overlap_jag_khpA_RNAseq.txt")
writeLines(overlap_All, "overlap_All_RNAseq.txt")
writeLines(overlap_jag_kre, "overlap_jag_kre_RNAseq.txt")
writeLines(overlap_khpA_kre, "overlap_khpA_kre_RNAseq.txt")

#### CRAC
jag_CRAC = read.csv('/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/ratios/prepared/jag_ratios.txt_sig_genes_BSU_Numbers.txt', sep = '\t', header = FALSE)
khpA_CRAC = read.csv('/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/ratios/prepared/khpA_ratios.txt_sig_genes_BSU_Numbers.txt', sep = '\t', header = FALSE)
kre_CRAC = read.csv('/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/ratios/prepared/kre_ratios.txt_sig_genes_BSU_Numbers.txt', sep = '\t', header = FALSE)

x <- list(Jag=jag_CRAC$V1, KhpA=khpA_CRAC$V1, Kre=kre_CRAC$V1)
#x <- list(Jag=jag_CRAC$V1, KhpA=khpA_CRAC$V1, Kre=kre_CRAC$V1)


ggp <- ggvenn(
  x,
  fill_color = c('#E31A1C', '#1F78B4', "#33A02C"),
  stroke_size = 0.5, set_name_size = 8, text_size = 8, show_percentage = FALSE
)

# Adjust the margins using theme
ggp <- ggp + theme(
  plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")  # Adjust these values as needed
)

ggp

# Specify the output directory and file name
output_dir <- "/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/figures/"
output_file <- file.path(output_dir, "Venn_CRAC.pdf")
# Save the plot as a PDF file
ggsave(output_file, plot = ggp, width = 8, height = 6)


### identify overlapping genes
overlap_All <- Reduce(intersect, list(jag_CRAC$V1[1:100], khpA_CRAC$V1[1:100], kre_CRAC$V1[1:100]))
overlap_All

overlap_jag_khpA <- Reduce(intersect, list(jag_CRAC$V1, khpA_CRAC$V1))
overlap_jag_khpA

overlap_jag_kre <- Reduce(intersect, list(jag_CRAC$V1, kre_CRAC$V1))
overlap_jag_kre

overlap_khpA_kre <- Reduce(intersect, list(khpA_CRAC$V1, kre_CRAC$V1))
overlap_khpA_kre

jag_unique <- setdiff(setdiff(jag_CRAC$V1, overlap_jag_khpA), overlap_jag_kre)
jag_unique

khpA_unique <- setdiff(setdiff(khpA_CRAC$V1, overlap_jag_khpA), overlap_khpA_kre)
khpA_unique

kre_unique <- setdiff(setdiff(kre_CRAC$V1, overlap_khpA_kre), overlap_jag_kre)
kre_unique


### overlap between RNAseq and CRAC
## jag
x <- list("jag (RNA-seq)"=jag_RNA$V1, "Jag (CRAC)"=jag_CRAC$V1)

ggp <- ggvenn(
  x,
  fill_color = c('#E31A1C', '#1F78B4', "#33A02C"),
  stroke_size = 0.5, set_name_size = 8, text_size = 8, show_percentage = FALSE
)

ggp

overlap_jag <- Reduce(intersect, list(jag_RNA$V1, jag_CRAC$V1))
overlap_jag

# Specify the output directory and file name
output_dir <- "/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/figures/"
output_file <- file.path(output_dir, "Venn_jag.pdf")
# Save the plot as a PDF file
ggsave(output_file, plot = ggp, width = 8, height = 6)


## khpA
khpA_CRAC_list <- setdiff(khpA_CRAC$V1, "BSU_29340")
x <- list("khpA (RNA-seq)"=khpA_RNA$V1, "KhpA (CRAC)"=khpA_CRAC_list)

ggp <- ggvenn(
  x,
  fill_color = c('#E31A1C', '#1F78B4', "#33A02C"),
  stroke_size = 0.5, set_name_size = 8, text_size = 8, show_percentage = FALSE
)
ggp

overlap_khpA <- Reduce(intersect, list(khpA_RNA$V1, khpA_CRAC$V1))
overlap_khpA

# Specify the output directory and file name
output_dir <- "/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/figures/"
output_file <- file.path(output_dir, "Venn_khpA.pdf")
# Save the plot as a PDF file
ggsave(output_file, plot = ggp, width = 8, height = 6)

## kre
kre_RNA_list <- setdiff(kre_RNA$V1, "BSU_14020")
x <- list("kre (RNA-seq)"=kre_RNA_list, "Kre (CRAC)"=kre_CRAC$V1)
x <- list("kre (feedback_loop)"=kre_feedback$V1, "Kre (CRAC)"=kre_CRAC$V1)


ggp <- ggvenn(
  x,
  fill_color = c('#E31A1C', '#1F78B4', "#33A02C"),
  stroke_size = 0.5, set_name_size = 8, text_size = 8, show_percentage = FALSE
)
ggp

overlap_kre <- Reduce(intersect, list(kre_RNA$V1, kre_CRAC$V1))
overlap_kre <- Reduce(intersect, list(kre_feedback$V1, kre_CRAC$V1))
overlap_kre

# Specify the output directory and file name
output_dir <- "/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/figures/"
output_file <- file.path(output_dir, "Venn_kre.pdf")
# Save the plot as a PDF file
ggsave(output_file, plot = ggp, width = 8, height = 6)


### overlap in Top cross-linked genes NO RATIOS
##### RNAseq
jag_RNA = read.csv('/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/ratios_NoRatios/RPKMPW_Jag_no_ratios.txt', sep = '\t', header = FALSE)
khpA_RNA = read.csv('/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/ratios_NoRatios/RPKMPW_KhpA_no_ratios.txt', sep = '\t', header = FALSE)
kre_RNA = read.csv('/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/ratios_NoRatios/RPKMPW_Kre_no_ratios.txt', sep = '\t', header = FALSE)

x <- list(jag=jag_RNA$V1, khpA=khpA_RNA$V1, kre=kre_RNA$V1)

ggp <- ggvenn(
  x,
  fill_color = c('#E31A1C', '#1F78B4', "#33A02C"),
  stroke_size = 0.5, set_name_size = 8, text_size = 8, show_percentage = FALSE
)
ggp

### identify overlapping genes
overlap_All <- Reduce(intersect, list(jag_RNA$V1, khpA_RNA$V1, kre_RNA$V1))
overlap_All

setwd('/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/ratios_NoRatios')
writeLines(overlap_All, "overlap_All.txt")
