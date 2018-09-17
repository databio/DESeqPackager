# DESeq-Packager
From Aaron Gu

`DESeq_Packager(pepr, data_source, gene_names, gene_counts)`

DESeq-Packager takes in an RNA-seq project in [PEP format](https://pepkit.github.io/docs/pepr/) and combines the data for each sample into a DESeq countDataSet structure. Best used with output from the [rnapipe](https://github.com/databio/rnapipe) pipeline.

### Required Packages

- pepr
- data.table

## Quick Start

The test_pckgr.html file contains the instructions to use the DESeq-Packager function. In a few seconds, DESeq-Packager will produce a countDataSet ready for subsequent DESeq analysis.

