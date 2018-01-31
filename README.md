# DESeq-Packager
From Aaron Gu

DESeq-Packager takes in a pepr project (https://pepkit.github.io/docs/pepr/) and puts all of the data from the files into a countTable data structure, which is necessary for DESeq.

To use DESeq-Packager, create a YAML file for your project annotation sheet, then specify the yaml file and the data source column name to the DESeq-Packager function. As an example: `countTable <- DESeq_table_maker("Documents/GitHub/DESeq-Packager/project_config.yaml", "data_source")` 
