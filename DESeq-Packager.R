source("https://bioconductor.org/biocLite.R")
biocLite("DESeq")
library("DESeq")

#' packages multiple data tables into one data frame for DESeq analysis
#'
#' @param yaml file path to project's yaml
#' @return countTable, the data structure needed to input into DESeq
#' @export
DESeq_table_maker <- function(yaml, data_source){
  
  #importing pepr project
  devtools::install_github("pepkit/pepr")
  library(pepr)
  p <- Project(file = yaml)
  sample_frame <- samples(p)
  files <- sample_frame[ , data_source]
  type <- sample_frame[ , sample_name]
  
  #prompting user to specify columns for DESeq Analysis
  gene_name_col <- readline(prompt= "Specify which column in each data table is the gene name column: ")
  relevant_data_col <- readline(prompt= "Specify which column contains the relevant data for DESeq analysis: ")
  
  print("Packaging...")
  
  #reading in the files
  df_list <- vector(mode="list", length=length(files))
  for(i in 1:length(files)){
    sampleTable <- read.table(files[i], header = TRUE)
    sampleTable <- sampleTable[, c(gene_name_col, relevant_data_col)]
    sampleTable[[2]] <- as.integer(sampleTable[[2]]) #integers are required for DESeq
    names(sampleTable)[2] <- type[i]
    df_list[[i]] <- sampleTable
  }
  
  #function to merge the datatables in the list
  countTable <- merge(df_list[[1]], df_list[[2]], by = gene_name_col)
  for(i in 3:length(df_list)){
    countTable <- merge(countTable, df_list[[i]])
  }

  #set the row names as the gene names, then remove the gene name column
  row.names(countTable) <- countTable[[1]]
  countTable[,1] <- NULL
  
  print("Finished! View the new data frame in the variable called 'countTable'")
  return(countTable)
}

countTable <- DESeq_table_maker("Documents/GitHub/DESeq-Packager/project_config.yaml", "data_source")



#-----------------DESeq Analysis-------------------------#
#set the conditions......this could be passed in as a parameter to an enclosing function
condition <- factor(c("controlDMSO","controlDMSO","knockoutDMSO","knockoutDMSO"))
#run DESeq analysis
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