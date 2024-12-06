# calc. GSEA

## import library
library(data.table)
library(dplyr)
library(stringr)
library(msigdbr)
library(org.Hs.eg.db)
library(clusterProfiler)

## formatteing gene sets
all_gene_sets = msigdbr(species = "Homo sapiens")
H <- all_gene_sets %>% filter(gs_cat == "H")
others <- all_gene_sets %>% filter(gs_subcat %in% c("CP:BIOCARTA","CP:KEGG","CP:PID","CP:REACTOME","CP:WIKIPATHWAYS","GO:BP","GO:CC","GO:MF"))
test_gene_set <- rbind(H,others)
test_gene_set <- test_gene_set %>% dplyr::select(gs_name,entrez_gene)

## run ClusterProfiler
input <- path/to/DESeq2Results/DESeq2Results.txt
table <- fread(input)

geneList <- table$SYMBOL
geneList.df <- bitr(geneList, fromType = "SYMBOL",
                    toType = c("ENSEMBL", "ENTREZID"),
                    OrgDb = org.Hs.eg.db)
table <- left_join(table,geneList.df,by="SYMBOL")
table <- table %>%  distinct(ENTREZID,.keep_all=TRUE) %>% na.omit() 
gsea_input <- table$log2FoldChange
names(gsea_input) <- table$ENTREZID
gsea_input <- sort(gsea_input,decreasing = T)
gsea_out <- GSEA(gsea_input, TERM2GENE = test_gene_set,eps=0,seed=1,pvalueCutoff = 1,pAdjustMethod = "BH")
gsea_out <- setReadable(gsea_out, 'org.Hs.eg.db', 'ENTREZID')

## save results
output <- gsea_out %>% data.frame()
x <- output$ID %>% str_split("_",simplify=T)
output <- output %>% mutate(geneset_name = factor(x[,1]))

OUTPUT_PATH <- path/to/output
fwrite(output,paste0(OUTPUT_PATH,"/GSEA_results.txt"))
saveRDS(gsea_out,paste0(OUTPUT_PATH,"/GSEA_results.rds"))
