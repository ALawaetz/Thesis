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

###############################################################################
### Estimate sice factors for merged files
directory <- '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/Coverage_bed'
sampleFiles <- list.files(path=directory, pattern='ALL_merged_cut.txt')
sampleFiles
sampleNames <- c('jag', 'khpA', 'kre', 'WT')
sampleTable<-data.frame(sampleName=sampleNames, fileName=sampleFiles)
sampleTable$strain <- as.factor(sampleTable$sampleName)
sampleTable

dds<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~ strain)

### View orginal counts table
#View(counts(dds))

### estimate size factors for normalising bedgraph files
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
### beware that KhpA is much smaller because it has been merged from only 2 and not 3 samples
#normalized_counts <- counts(dds, normalized=TRUE)
#View(normalized_counts)

#colData(dds)$strain<-factor(colData(dds)$strain, levels=c('WT','kre'))

###############################################################################
### Estimate size factors for all strains/replicates in the same comparison
directory <- '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/Coverage_bed'
sampleFiles <- list.files(path=directory, pattern='txt')
# Filter out files that match the exclude pattern
sampleFiles <- sampleFiles[!grepl('ALL_merged', sampleFiles)]
sampleFiles <- sampleFiles[!grepl('khpA_1', sampleFiles)]
sampleFiles
sampleNames <- c('jag_1', 'jag_2', 'jag_3', 'khpA_2', 'khpA_3', 'kre_1', 'kre_2', 'kre_3', 'WT_1', 'WT_2', 'WT_3')
OD = c(2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 1)
Names <- c('jag_1', 'jag_2', 'jag_3', 'khpA_2', 'khpA_3', 'kre_1', 'kre_2', 'kre_3', 'WT_1', 'WT_2', 'WT_3')
Strain <- c('jag', 'jag', 'jag', 'khpA', 'khpA', 'kre', 'kre', 'kre', 'WT', 'WT', 'WT')
replicate <- c(1, 2, 3, 2, 3, 1, 2, 3, 1, 2, 3)

sampleTable<-data.frame(sampleName=sampleNames, fileName=sampleFiles, names = Names, strain=Strain, OD = OD, replicate = replicate)
sampleTable$OD <- as.factor(sampleTable$OD)
sampleTable$strain <- as.factor(sampleTable$strain)
sampleTable

dds<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~ OD + strain)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

###############################################################################
### Make PCA plot with all strains in the same plot #niiice
dds<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~ OD + strain)
dds <- DESeq(dds)

### PCA plot
rld <- rlogTransformation(dds, blind=TRUE) 
pca_plot <- plotPCA(rld, intgroup = c("strain"), ntop=5000)
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
output_file <- file.path(output_dir, "PCA_plot_ALL.pdf")
# Save the plot as a PDF file
ggsave(output_file, plot = pca_plot, width = 8, height = 6)


###############################################################################
### Estimate size factors for WT/jag_khpA merged files
directory <- '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/Coverage_bed'
sampleFiles <- list.files(path=directory, pattern='ALL_merged_cut.txt')
sampleFiles <- sampleFiles[!grepl('kre', sampleFiles)]
sampleFiles
sampleNames <- c('jag', 'khpA', 'WT')
sampleTable<-data.frame(sampleName=sampleNames, fileName=sampleFiles)
sampleTable$strain <- as.factor(sampleTable$sampleName)
sampleTable

dds<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~ strain)

### View orginal counts table
#View(counts(dds))

### estimate size factors for normalising bedgraph files
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

###############################################################################
### Estimate size factors for WT/jag_khpA NOT merged files
directory <- '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/Coverage_bed'
sampleFiles <- list.files(path=directory, pattern='txt')
exclude_patterns <- 'kre|khpA_1|ALL_merged'
sampleFiles <- sampleFiles[!grepl(exclude_patterns, sampleFiles)]
sampleNames <- c('jag_1', 'jag_2', 'jag_3', 'khpA_2', 'khpA_3', 'WT_1', 'WT_2', 'WT_3')
Strain <- c('jag', 'jag', 'jag', 'khpA', 'khpA', 'WT', 'WT', 'WT')
sampleTable<-data.frame(sampleName=sampleNames, fileName=sampleFiles, strain=Strain)
sampleTable$strain <- as.factor(sampleTable$strain)
sampleTable

dds<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~ strain)

### View orginal counts table
#View(counts(dds))

### estimate size factors for normalising bedgraph files
dds <- estimateSizeFactors(dds)
sizeFactors(dds)


###############################################################################
### Estimate size factors for WT/kre merged files
directory <- '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/Coverage_bed'
sampleFiles <- list.files(path=directory, pattern='ALL_merged_cut.txt')
exclude_patterns <- 'jag|khpA'
sampleFiles <- sampleFiles[!grepl(exclude_patterns, sampleFiles)]
sampleFiles
sampleNames <- c('kre', 'WT')
sampleTable<-data.frame(sampleName=sampleNames, fileName=sampleFiles)
sampleTable$strain <- as.factor(sampleTable$sampleName)
sampleTable

dds<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~ strain)

### View orginal counts table
#View(counts(dds))

### estimate size factors for normalising bedgraph files
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

###############################################################################
### Estimate size factors for WT/jag_khpA NOT merged files
directory <- '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/Coverage_bed'
sampleFiles <- list.files(path=directory, pattern='txt')
exclude_patterns <- 'jag|khpA|ALL_merged'
sampleFiles <- sampleFiles[!grepl(exclude_patterns, sampleFiles)]
sampleNames <- c('kre_1', 'kre_2', 'kre_3', 'WT_1', 'WT_2', 'WT_3')
Strain <- c('kre', 'kre', 'kre', 'WT', 'WT', 'WT')
sampleTable<-data.frame(sampleName=sampleNames, fileName=sampleFiles, strain=Strain)
sampleTable$strain <- as.factor(sampleTable$strain)
sampleTable

dds<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~ strain)

### View orginal counts table
#View(counts(dds))

### estimate size factors for normalising bedgraph files
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
