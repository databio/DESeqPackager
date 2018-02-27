#' packages multiple data tables into one data frame for DESeq analysis
#'
#' @param yaml file path to project's yaml
#' @param data_source title of the sample annotation sheet column with file names
#' @param gene_name_col name of the column with unique gene identifiers
#' @param relevant_data_col name of the column with the relevant data for DESeq
#' @return countTable, the data structure needed for DESeq
#' @export
DESeq_table_maker <- function(yaml, data_source, gene_name_col, relevant_data_col){
  
  #importing pepr project
  devtools::install_github("pepkit/pepr")
  library(pepr)
  p <- Project(file = yaml)
  sample_frame <- samples(p)
  files <- sample_frame[ , data_source]
  type <- sample_frame[ , sample_name]
  
  dt_list <- vector(mode="list",length=length(files))
  
  #checking for data table dependency and reading in files
  if (requireNamespace("data.table")) {
    print("in fread")
    for(i in 1:length(files)){
      sampleTable <- data.table::fread(files[i])
      sampleTable <- sampleTable[, .(get(gene_name_col), get(relevant_data_col))]
      sampleTable[[2]] <- as.integer(sampleTable[[2]])
      names(sampleTable)[1] <- gene_name_col
      names(sampleTable)[2] <- type[i]
      dt_list[[i]] <- sampleTable
    }
  } else {
    print("in read.table")
    for(i in 1:length(files)){
      sampleTable <- read.table(files[i], header = TRUE)
      sampleTable <- sampleTable[, c(gene_name_col, relevant_data_col)]
      sampleTable[[2]] <- as.integer(sampleTable[[2]]) #integers are required for DESeq
      names(sampleTable)[2] <- type[i]
      df_list[[i]] <- sampleTable
    }
  }
  print("done reading")
  #function to merge the datatables in the list
  countTable <- merge(dt_list[[1]], dt_list[[2]], by = gene_name_col)
  for(i in 3:length(dt_list)){
    countTable <- merge(countTable, dt_list[[i]])
  }
  
  #set the row names as the gene names, then remove the gene name column
  row.names(countTable) <- countTable[[1]]
  countTable[,1] <- NULL

  return(countTable)
}

countTable <- DESeq_table_maker("${PROCESSED}/DESeq-Packager/project_config.yaml", "data_source", "ensembl_gene_id", "FPKM")



#------------------------------DESeq Analysis-----------------------------------#

#the user will conduct the DESeq analysis themselves
#this code is for my testing purposes
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq")
library("DESeq")


condition <- factor(c("controlDMSO","controlDMSO","knockoutDMSO","knockoutDMSO"))
cds <- newCountDataSet(countTable, condition)
cds <- estimateSizeFactors(cds)

cds <- estimateDispersions(cds, fitType = "local")
plotDispEsts(cds)

res <- nbinomTest(cds, "controlDMSO", "knockoutDMSO")
plotMA(res)

hist(res$pval, breaks =100, col='skyblue', border = 'slateblue', main = '')

#filter for significant genes
resSig = res[res$padj<0.1, ]
head(resSig[order(resSig$pval), ])