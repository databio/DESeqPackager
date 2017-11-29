install.packages(c("readr", "data.table", "plyr"))
library("readr")
library("data.table")
library("plyr")
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq")
library("DESeq")

#What do I do with the working directory?
setwd("/Users/AG/R/expr_tsv")


####parameters

#files <- c("RNA_EWS-FLI1_High_rep1.tsv", "RNA_EWS-FLI1_High_rep2.tsv", "RNA_EWS-FLI1_Low_rep1.tsv", "RNA_EWS-FLI1_Low_rep2.tsv")
files <- c("RNA_EWS-FLI1_High_DMSO_rep1.tsv", "RNA_EWS-FLI1_High_DMSO_rep2.tsv", "RNA_EWS-FLI1_High_rep1.tsv", "RNA_EWS-FLI1_High_rep2.tsv")
#type <- c("High1", "High2", "Low1", "Low2")
type <- c("HighDMSO1", "HighDMSO2", "High1", "High2")

gene_name_col <- "ensembl_gene_id"
relevant_data_col <- "FPKM"


####function to read in the files
setup_datafiles <- function(file_vector, type, gene_name_col, relevant_data_col){
  output <- vector(mode = "list", length = length(file_vector)) #creates empty list
  
  i <- 1
  for(file in file_vector){
    dt <- fread(file)
    
    dt <- subset(dt, select = c(gene_name_col, relevant_data_col)) #leave only the geneID and the relevant data (FPKM)
    
    dt[[2]] <- as.integer(dt[[2]]) #convert the FPKM data to integers, required for DESeq
    
    names(dt)[2] <- paste(type[i], relevant_data_col, sep = "_") #rename columns in datatable to the type of data (ex: HighDMSO1_FPKM)
    
    output[[i]] <- dt #add datatable to the list
    
    i <- i+1 #increment to access the next element in the type vector
  }
  return(output)
}

#run the function and form a list of the given files
list_files <- setup_datafiles(files, type, gene_name_col, relevant_data_col)

#function to merge the datatables in the list
merging <- function(list){
  countTable <- merge(list[[1]], list[[2]], by = gene_name_col)
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
condition2 <- factor(c("treatment", "treatment", "control", "control"))



#run DESeq analysis
cds <- newCountDataSet(countTable, condition)
cds2 <- newCountDataSet(countTable2, condition2)
cds2 <- estimateSizeFactors(cds2)

cds <- estimateDispersions(cds, fitType = "local")
plotDispEsts(cds)

res <- nbinomTest(cds, "knockout", "control")
head(res)
plotMA(res)

hist(res$pval, breaks =100, col='skyblue', border = 'slateblue', main = '')

#filter for significant genes
resSig = res[res$padj<0.1, ]
head(resSig[order(resSig$pval), ])