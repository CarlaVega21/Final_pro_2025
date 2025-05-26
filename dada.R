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
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")
BiocManager::install("dada2")

plotQualityProfile(fnFs[1:2])

filt_path <- file.path(path, "filtered")
dir.create(filt_path, showWarnings = FALSE)
filtFs <- file.path(filt_path, basename(fnFs))


out <- filterAndTrim(fnFs, filtFs,
                     truncLen=250,  
                     maxN=0, maxEE=2, truncQ=2, 
                     compress=TRUE, multithread=TRUE)



errF <- learnErrors(filtFs, multithread=TRUE)


dds <- dada(filtFs, err=errF, multithread=TRUE)


seqtab <- makeSequenceTable(dds)


seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)


summary(seqtab.nochim)

