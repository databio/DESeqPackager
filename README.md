# DESeq-Packager
From Aaron Gu

```R
DESeq_Packager(pepr, data_source, gene_names, gene_counts)
```

DESeq-Packager takes in an RNA-seq project in [PEP format](https://pepkit.github.io/docs/pepr/) and combines the data for each sample into a DESeq countDataSet structure. It's best used with output from the [rnapipe](https://github.com/databio/rnapipe) pipeline.

### Required Packages

- pepr
- data.table

## Quick Start

The quick start uses sample gene expression data available for download at `expr_tsv.tar.gz` from http://big.databio.org/example_data/deseq_packager/

In a new R file, load in the function
```R
source("DESeq-Packager.R")
```

Construct a pepr project using the yaml located in this repository
```R
p = pepr::Project(file="project_config.yaml")
```

Run DESeq-Packager
```R
countDataSet <- DESeq_Packager(p, "data_source", "ensembl_gene_id", "FPKM")
```

### Running DESeq-Packager with BioConductor

[BiocProject](http://code.databio.org/BiocProject/index.html) is another way to read in biological data in the PEP format. It can load both project metadata and sample data in a single line of code.

```R
bpArgs = BiocProject(file="project_config.yaml", funcArgs=list(data_source="data_source", gene_names="ensembl_gene_id", gene_counts="FPKM"))

getData(bpArgs)
```
