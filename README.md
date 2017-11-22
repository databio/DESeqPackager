# DESeq-Packager
From Aaron Gu

The main function in the project is `setup_datafiles()`, which will eventually take in a pepr project. For now, it is coded to take in a vector of file names and a vector of the names of the samples.
After the datafiles are read in, they are merged into one countTable using the function `merging()`. The countTable is the essential data structure that can be run through the various DESeq analyses. I have included the DESeq functions for graphing the expression change (`plotDispEsts` and `plotMA`) as well as the table of genes with significant differential expression (`resSig`) in this code.
