# DESeq-Packager
From Aaron Gu

`DESeq_Packager(pepr, data_source, gene_names, gene_counts)`

DESeq-Packager takes in an RNA-seq project in [PEP format](https://pepkit.github.io/docs/pepr/) and combines the data for each sample into a DESeq countDataSet structure. It's best used with output from the [rnapipe](https://github.com/databio/rnapipe) pipeline.

### Required Packages

- pepr
- data.table

## Quick Start

The quick start uses sample gene expression data available for download at `expr_tsv.tar.gz` from big.databio.org/example_data/

In a new R file, load in the function
```R
source("DESeq-Packager.R")
```

Construct a pepr project using the yaml
```R
p = pepr::Project(file="project_config.yaml")
```

Run DESeq-Packager
```R
countDataSet <- DESeq_Packager(p, "data_source", "ensembl_gene_id", "FPKM")
```