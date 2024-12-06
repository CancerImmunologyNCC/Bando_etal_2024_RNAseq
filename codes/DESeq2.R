# Differentially Expressed Genes analysis

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(magrittr)
library(edgeR)
library(ggplot2)
library(DESeq2)

GeneCount <- "path/to/GeneCount.csv"
metadata <- "path/to/metadata"
OUTPUT_PATH <- "path/to/output"

################################################################################
# data import
dat <- fread(GeneCount)
metadata <- fread(metadata)

## Metadata Contents
##SampleID  PatientID Treatment CR  SequenceBatch
##TENERGY1-1  TENERGY1  1 1 HN00156398
##TENERGY1-2  TENERGY1  2 1 HN00156398
##TENERGY1-3  TENERGY1  3 1 HN00156398

## Treatment : 1 -> pre, 2 -> afterCRT, 3 -> afterAtezolizumab
## CR : 1 -> CR, 0 -> nonCR

################################################################################
# data formatting
dat <- dat %>% separate(gene_id,c("ENSG","SYMBOL"),sep="\\|") %>% select(- ENSG)
dat <- dat %>% group_by(SYMBOL) %>% summarise(across(everything(), sum) )


################################################################################
# Calculation differential expressed genes by DESeq2

## the example about CR vs nonCR in Pre patients

### Select samples
group <- metadata %>% 
         filter(Treatment == 1) %>% select(SampleID,CR)  %>% arrange(- CR) %>% mutate(CR=factor(CR))
testdat <- dat[,group$SampleID]

### filtered out low expressed genes
countfilter <- testdat %>% apply(1, mean) 
countfilter <- countfilter > 1
testdat <- set_rownames(x = testdat[countfilter,], value = rownames(dat)[countfilter])

### runDEseq2
group <- group %>% select(CR)
dds <- DESeqDataSetFromMatrix(countData = testdat, colData = group, design = ~ CR)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomLRT(dds, full = ~ CR, reduced = ~ 1)
res <- results(dds,tidy=TRUE)
res <- cbind(res,testdat) 
colnames(res)[1] <- "SYMBOL"

### save file
fwrite(res,paste0(OUTPUT_PATH,"CRvsNonCRinPre_DESeq2_results.txt"),sep="\t")

################################################################################
## the example about afterCRT vs Pre in all patients

### Select samples
group <- metadata %>% 
         filter(Treatment!=3) %>% group_by(PatientID) %>% mutate(n=n()) %>%  filter(n == 2) %>%
         select(SampleID,PatientID, Treatment) %>% arrange( - Treatment) %>% mutate(Treatment=factor(Treatment))
testdat <- dat[,group$SampleID]

### filtered out low expressed genes
countfilter <- testdat %>% apply(1, mean) 
countfilter <- countfilter > 1
testdat <- set_rownames(x = testdat[countfilter,], value = rownames(dat)[countfilter])
group <- group %>% select(PatientID, Treatment) %>% mutate(PatientID=factor(PatientID))

### runDEseq2
dds <- DESeqDataSetFromMatrix(countData = testdat, colData = group, design = ~ PatientID + Treatment)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomLRT(dds, full = ~ PatientID +  Treatment, reduced = ~ PatientID)
res <- results(dds,tidy=TRUE)
res <- cbind(res,testdat) 
res %>% arrange(padj,pvalue) %>% head(10)
colnames(res)[1] <- "SYMBOL"

fwrite(res,paste0(OUTPUT_PATH,"afterCRTvsPreinall_DESeq2_results.txt"),sep="\t")







