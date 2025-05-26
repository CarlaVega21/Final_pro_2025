

BiocManager::install("DESeq2")
library("DESeq2")
library(phyloseq)
library(microbiome)
library(tidyverse)
library(dplyr)
library(dplyr)
library(tibble)


comme_counts <- read.table("GSE286382_rawreads_aligned_to_apis.txt/GSE286382_rawreads_aligned_to_apis1.txt", 
                           header = TRUE, row.names = 1, sep = "\t")

readLines("GSE286382_rawreads_aligned_to_apis.txt/GSE286382_rawreads_aligned_to_apis1.txt", n = 10)

header_line <- readLines("GSE286382_rawreads_aligned_to_apis.txt/GSE286382_rawreads_aligned_to_apis1.txt", n = 1)
header <- strsplit(header_line, "\t")[[1]]



meta_cols <- grep("^#CLASS:", header)
raw_data <- read.table("GSE286382_rawreads_aligned_to_apis.txt/GSE286382_rawreads_aligned_to_apis1.txt",
                       sep = "\t", header = FALSE, skip = 1, fill = TRUE, quote = "", comment.char = "")



colnames(raw_data) <- header

metadata <- raw_data[ , 1:(length(meta_cols) + 1)]  # +1 por "Sample"
expr_data <- raw_data[ , -(1:(length(meta_cols) + 1))]

dim(expr_data)
head(expr_data[ , 1:5])


expr_matrix <- as.matrix(sapply(expr_data, as.numeric))
rownames(expr_matrix) <- metadata$Sample

colnames(metadata) <- gsub("#CLASS:", "", colnames(metadata)) 
colnames(metadata) <- make.names(colnames(metadata)) 
colnames(metadata)
class(metadata)
metadata <- as.data.frame(metadata)




meta_selected <- metadata_df %>%
  dplyr::select(Sample, Age, Carb)

sample_metadata <- column_to_rownames(meta_selected, var = "Sample")

head(sample_metadata)

dim(expr_matrix)


head(rownames(expr_matrix))
head(colnames(expr_matrix))


expr_matrix_t <- t(expr_matrix)

library(DESeq2)
colnames(sample_metadata)


all(rownames(sample_metadata) %in% colnames(expr_matrix_t))
all(colnames(expr_matrix_t) %in% rownames(sample_metadata))
head(rownames(sample_metadata))
head(colnames(expr_matrix_t))

print(colnames(sample_metadata))


print(rownames(sample_metadata))

print(colnames(expr_matrix_t))

print(dim(sample_metadata))
print(dim(expr_matrix_t))

metadata_df <- as.data.frame(metadata)

sample_metadata <- metadata_df %>%
  dplyr::select(Sample, Age, Carb) %>%
  column_to_rownames(var = "Sample")

print(colnames(sample_metadata))

dds <- DESeqDataSetFromMatrix(countData = expr_matrix_t,
                              colData = sample_metadata,
                              design = ~ Age + Carb)
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj), ]
head(resOrdered)


resSig <- subset(resOrdered, padj < 0.05)
head(resSig)

library(ggplot2)
res_df <- as.data.frame(res)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  labs(title = "Volcano plot", x = "log2 Fold Change (Carb)", y = "-log10(padj)")



resAge <- results(dds, contrast = c("YOUNG", "OLD"))

resAgeSig <- subset(resAge, padj < 0.05)


#################################3
dada
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dada2")
library(dada2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ShortRead)
library(Biostrings)
library(stats)
path <- "C:/Users/crism/Documents/proyecto final robert/proyecto final"
setwd(path)
list.files()
fnFs <- sort(list.files(path, pattern = "\\.fastq\\.gz$", full.names = TRUE))


length(fnFs)  
head(fnFs)  


filt_path <- file.path(path, "filtered")
dir.create(filt_path, showWarnings = FALSE)
filtFs <- file.path(filt_path, basename(fnFs))

plotQualityProfile(fnFs[1:2])

out <- filterAndTrim(fnFs, filtFs,
                     truncLen=250,  
                     maxN=0, maxEE=2, truncQ=2, 
                     compress=TRUE, multithread=TRUE)



errF <- learnErrors(filtFs, multithread=TRUE)


dds <- dada(filtFs, err=errF, multithread=TRUE)


seqtab <- makeSequenceTable(dds)


seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)


summary(seqtab.nochim)
raw
################otra vez
readLines("GSE286382_rawreads_aligned_to_apis.txt/GSE286382_rawreads_aligned_to_apis1.txt", n = 10)

header_line <- readLines("GSE286382_rawreads_aligned_to_apis.txt/GSE286382_rawreads_aligned_to_apis1.txt", n = 1)
header <- strsplit(header_line, "\t")[[1]]



meta_cols <- grep("^#CLASS:", header)
raw_data <- read.table("GSE286382_rawreads_aligned_to_apis.txt/GSE286382_rawreads_aligned_to_apis1.txt",
                       sep = "\t", header = FALSE, skip = 1, fill = TRUE, quote = "", comment.char = "")



colnames(raw_data) <- header

metadata <- raw_data[ , 1:(length(meta_cols) + 1)]  # +1 por "Sample"
expr_data <- raw_data[ , -(1:(length(meta_cols) + 1))]

dim(expr_data)
head(expr_data[ , 1:5])


expr_matrix <- as.matrix(sapply(expr_data, as.numeric))
rownames(expr_matrix) <- metadata$Sample

colnames(metadata) <- gsub("#CLASS:", "", colnames(metadata))  # elimina "#CLASS:"
colnames(metadata) <- make.names(colnames(metadata))  # hace nombres válidos en R
colnames(metadata)
class(metadata)
metadata <- as.data.frame(metadata)

dim(expr_matrix)


head(rownames(expr_matrix))
head(colnames(expr_matrix))


expr_matrix_t <- t(expr_matrix)

library(DESeq2)


metadata_df <- as.data.frame(metadata)

sample_metadata <- metadata_df %>%
  dplyr::select(Sample, Age, Carb) %>%
  column_to_rownames(var = "Sample")

print(colnames(sample_metadata))

dds <- DESeqDataSetFromMatrix(countData = expr_matrix_t,
                              colData = sample_metadata,
                              design = ~ Age + Carb)
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj), ]
head(resOrdered)


resSig <- subset(resOrdered, padj < 0.05)
head(resSig)

library(ggplot2)
res_df <- as.data.frame(res)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  labs(title = "Volcano plot", x = "log2 Fold Change (Carb)", y = "-log10(padj)")

library(pheatmap)

