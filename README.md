# DESeq-Packager
From Aaron Gu

`DESeq_Packager(yaml, data_source, gene_names, gene_counts)`

DESeq-Packager takes in an RNA sequencing project in PEP format (https://pepkit.github.io/docs/pepr/) and combines the data for each sample into a DESeq countDataSet structure.

## Quick Start

To try it on the example data provided, clone the Git repository, and run the first two lines of code in the test_pckgr.R file. In a few seconds, the countDataSet will be produced and ready for subsequent DESeq analysis.

Otherwise just copy the DESeq-Packager.R file into your working directory and make sure to include the `source("DESeq-Packager.R")` before you use the function.
