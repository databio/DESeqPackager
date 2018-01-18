install.packages("data.table")
library("data.table")
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq")
library("DESeq")

#importing pepr project
devtools::install_github("pepkit/pepr")
library(pepr)
p <- Project(file = "project_config.yaml")
sample_frame <- samples(p)
files <- sample_frame[ , data_source]

#prompting user to rename files 
parameters <- c()
while(TRUE){
  s <- readline(prompt= "Rename each file with sample identifiers (Type DONE when finished): ")
  if(s=="DONE"){
    rm(s)
    break
  }
  parameters <- c(parameters, s)
}
type <- c()
for(column in parameters){
  type <- paste(type, sample_frame[, get(column)], sep = "")
}
rm(parameters, column)

#prompting user to specify columns for DESeq Analysis
gene_name_col <- readline(prompt= "Specify which column in each data table is the gene name column: ")
relevant_data_col <- readline(prompt= "Specify which column is the relevant data for DESeq analysis: ")

#reading in the files
df_list <- vector(mode="list", length=length(files))
for(i in 1:length(files)){
  sampleTable <- read.table(files[i], header = TRUE)
  sampleTable <- sampleTable[, c(gene_name_col, relevant_data_col)]
  sampleTable[[2]] <- as.integer(sampleTable[[2]]) #integers are required for DESeq
  names(sampleTable)[2] <- paste(type[i], relevant_data_col, sep = "_")
  df_list[[i]] <- sampleTable
}

#function to merge the datatables in the list
countTable <- merge(df_list[[1]], df_list[[2]], by=gene_name_col)
for(i in 3:length(df_list)){
  countTable <- merge(countTable, df_list[[i]])
}
#set the row names as the gene names, then remove the gene name column
row.names(countTable) <- countTable[[1]]
countTable[,1] <- NULL

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