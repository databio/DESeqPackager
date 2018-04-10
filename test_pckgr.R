Sys.setenv(CODEBASE="~/Documents/")
src <- paste(Sys.getenv("CODEBASE"),"DESeq-Packager/DESeq-Packager.R", sep="")
source(src)

project_config = paste(Sys.getenv("CODEBASE"), "DESeq-Packager/project_config.yaml", sep="")
countDataSet <- DESeq_Packager(project_config, "data_source", "ensembl_gene_id", "FPKM")
save("countDataSet", file=paste(Sys.getenv("CODEBASE"),"DESeq-Packager/.RData", sep=""))


# after running DESeq_Packager, the user will conduct the DESeq analysis themselves
# you can try a sample DESeq analysis by uncommenting the code below
if(FALSE){
  
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
  
  hist(res$pval, breaks=100, col='skyblue', border ='slateblue', main = '')
  
  #filter for significant genes
  resSig = res[res$padj<0.1, ]
  head(resSig[order(resSig$pval), ])

}