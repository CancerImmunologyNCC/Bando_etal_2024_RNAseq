# Volcano plot codes

library(data.table)
library(dplyr)
library(ggplot2)

DESeq2result <- "/path/to/DESeq2"
dat <- paste0(DESeq2result,"/DESeq2_result.txt") %>% fread()

# the code example about afterCRT vs pre in all patients 

dat %>% filter(!is.na(padj)) %>%
        mutate(minuslog10FDR = - log10(padj)) %>%
        mutate(category = case_when(
                            log2FoldChange > 1  & padj < 0.05 ~ "postCRT",
                            log2FoldChange < -1 & padj < 0.05 ~ "Pre",
                            TRUE ~ "NotSig")) %>% 
        mutate(category=factor(category,levels=c("postCRT","Pre","NotSig"))) %>%
    ggplot(aes(x=log2FoldChange,y=minuslog10FDR,color=category,size=category,label=label)) +
        geom_point() +
        scale_size_manual(values=c("postCRT"=1,"Pre"=1,"NotSig"=0.2)) + 
        scale_colour_manual(values=c("postCRT"="darkred", "Pre"="darkblue", "NotSig"="gray")) +
        theme_classic() +
        geom_hline(yintercept = -log10(0.05),linetype = "dotted") +
        geom_vline(xintercept = c(-1,1),linetype = "dotted") +
        labs(y = "-log10FDR",
             x = "log2FC\n\nPre <---> postCRT") + 
        ggtitle("postCRT vs Pre in all")


