# ssGSEA

## import library
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(GSVA)

## import geneset
### https://www.nature.com/articles/s41591-020-1082-2 Sup.Table8 and Treg geneset from https://www.sciencedirect.com/science/article/pii/S0092867414016390
geneset <- "path/to/geneset/geneset.txt"
geneset <- fread(geneset,header=FALSE)
geneset <- rbind(geneset,data.frame("CD4_Regulatory_T_cells","FOXP3,LINC02694,IL5,CTLA4,IL32,GPR15,IL4"),use.names=FALSE)

colnames(geneset) <- c("geneset","genes")
gs <- as.list(geneset$genes) %>% str_split(pattern=",")
names(gs) <- geneset$geneset

## import data
TPM <- "path/to/tpm/TPM.txt"
TPM <- fread(TPM)

## data formatting
x <- (TPM[,-1] != 0) %>% rowSums()
x <- x != 0
TPM <- TPM[x,] 
logTPM <- log10(TPM[,-1] + 1) %>% as.matrix()
rownames(logTPM) <- TPM$Gene.Name

## run ssGSEA
results <- gsva(logTPM, gs, method=c("ssgsea"),verbose=FALSE)

## save results
OUTPUT_PATH <- "path/to/output"
out <- results %>% data.frame()

fwrite(out,paste0(OUTPUT_PATH,"/ssGSEA_results.txt"), row.names = TRUE,sep="\t")
