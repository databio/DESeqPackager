# DESeq-Packager
From Aaron Gu

`DESeq_Packager(pepr, data_source, gene_names, gene_counts)`

DESeq-Packager takes in an RNA-seq project in [PEP format](https://pepkit.github.io/docs/pepr/) and combines the data for each sample into a DESeq countDataSet structure. Best used with output from the [rnapipe](https://github.com/databio/rnapipe) pipeline.

## Quick Start

To try it on the example data provided, clone the Git repository. Open up the test_pckgr.R file and change the `CODEBASE` R environment variable to the location of your DESeq-Packager directory, and then run the next 4 lines of code. In a few seconds, DESeq-Packager will produce a the countDataSet ready for subsequent DESeq analysis.

To run the test from the command line, use the command `R < test_pckgr.R --no-save`. The countDataSet result should be saved to the .R_Data file.
