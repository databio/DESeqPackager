#' Packages RNA sequencing results from multiple data files into one data frame for DESeq analysis using the PEP format
#' (read more about the PEP format at https://pepkit.github.io/)
#'
#' @param yaml file path to project_config yaml
#' @param data_source name of the column in the sample annotation sheet with the file names
#' @param gene_names name of the column in each data file with gene names
#' @param gene_counts name of the column in each data file with the count data for DESeq
#' @return countDataSet, the data frame needed for DESeq
#' @export
DESeq_Packager <- function(yaml, data_source, gene_names, gene_counts){
  if(!requireNamespace("pepr"))
    devtools::install_github("pepkit/pepr")
  p <- pepr::Project(file = yaml)
  sample_frame <- pepr::samples(p)
  sample_names <- sample_frame[["sample_name"]]
  files <- sample_frame[ , data_source]
  
  #checking for data table dependency and reading in files into a list
  if (!requireNamespace("data.table")) {
    install.packages("data.table")
  }
  library(data.table)
  
  print("reading in tables...")
  
  #create a table for a sample, then put it into a list
  dt_list <- vector(mode="list",length=length(files))
  for(i in 1:length(files)){
    sampleTable <- fread(files[i])
    sampleTable <- sampleTable[, c(gene_names, gene_counts), with=FALSE]
    sampleTable[[gene_counts]] <- lapply(sampleTable[[gene_counts]], as.integer)
    colnames(sampleTable)[1] <- gene_names
    colnames(sampleTable)[2] <- sample_names[i]
    dt_list[[i]] <- sampleTable
  }
  
  print("merging samples into one table...")
  
  #join each sample table from the list
  countDataSet <- merge(dt_list[[1]], dt_list[[2]])
  for(i in 3:length(dt_list)){
    countDataSet <- merge(countDataSet, dt_list[[i]])
  }
  
  #set the row names as the gene names, then remove the gene name column
  row.names(countDataSet) <- countDataSet[[1]]
  countDataSet[,1] <- NULL
  
  print("packaged!")
  save("countDataSet", file="~/Documents/Databio/DESeq-Packager/.RData")
  return(countDataSet)
}