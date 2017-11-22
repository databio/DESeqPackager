install.packages(c("readr", "data.table", "plyr"))
library("readr")
library("data.table")
library("plyr")
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq")
library("DESeq")

#What do I do with the working directory?
setwd("/Users/AG/R/expr_tsv")


#vector of files and their names
files <- c("RNA_EWS-FLI1_High_rep1.tsv", "RNA_EWS-FLI1_High_rep2.tsv", "RNA_EWS-FLI1_Low_rep1.tsv", "RNA_EWS-FLI1_Low_rep2.tsv")
type <- c("High1", "High2", "Low1", "Low2")


#function to read in the files
setup_datafiles <- function(file_vector, type){
  output <- vector(mode = "list", length = length(file_vector))
  i <- 1
  for(file in file_vector){
    dt <- fread(file)
    dt <- subset(dt, select = c("ensembl_gene_id", "FPKM"))
    dt <- dt[, FPKM:=as.integer(FPKM)]
    dt <- rename(dt, replace = c("FPKM" = paste(type[i], "FPKM", sep = "_")))
    output[[i]] <- dt
    names(output)[i] <- file
    i <- i+1
  }
  return(output)
}

#run the function on the given list
list_files <- setup_datafiles(files, type)


#function to merge the datatables in the list
merging <- function(list){
  countTable <- merge(list[[1]], list[[2]], by = "ensembl_gene_id")
  for(i in 3:length(list)){
    countTable <- merge(countTable, list[[i]])
  }
  return(countTable)
}

#run the merge function
countTable <- merging(list_files)
#remove the names of the genes
countTable <- countTable[,2:ncol(countTable)]


#set the conditions......this could be passed in as a parameter to an enclosing function
condition <- factor(c("knockout", "knockout", "control", "control"))

#run DESeq analysis
cds <- newCountDataSet(countTable, condition)
cds <- estimateSizeFactors(cds)

cds <- estimateDispersions(cds, fitType = "local")
plotDispEsts(cds)

res <- nbinomTest(cds, "knockout", "control")
head(res)
plotMA(res)

hist(res$pval, breaks =100, col='skyblue', border = 'slateblue', main = '')

#filter for significant genes
resSig = res[res$padj<0.1, ]
head(resSig[order(resSig$pval), ])