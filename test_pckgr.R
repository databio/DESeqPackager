source("~/Documents/Databio/DESeq-Packager/DESeq-Packager.R")

countDataSet <- DESeq_Packager("~/Documents/Databio/DESeq-Packager/project_config.yaml", "data_source", "ensembl_gene_id", "FPKM")


#the user will conduct the DESeq analysis themselves
#this code is for my testing purposes
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq")
library("DESeq")


condition <- factor(c("controlDMSO","controlDMSO","knockoutDMSO","knockoutDMSO"))
cds <- newCountDataSet(countDataSet, condition)
cds <- estimateSizeFactors(cds)

cds <- estimateDispersions(cds, fitType = "local")
plotDispEsts(cds)

res <- nbinomTest(cds, "controlDMSO", "knockoutDMSO")
plotMA(res)

hist(res$pval, breaks =100, col='skyblue', border = 'slateblue', main = '')

#filter for significant genes
resSig = res[res$padj<0.1, ]
head(resSig[order(resSig$pval), ])