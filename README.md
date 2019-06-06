# DESEq2 Analysis wrapper

## Usage

```
   deseq2analysis::analyseFromConfig("test.parameters.yaml", getwd())
```

It expects as input data in the yaml config (see example/test.parameters.yaml)
 * a tabular file with counts in columns per sample, genes in rows : countstable
 * a tabular sample description file : groupingtable
 * a tabular file specifying the comparisons done with the groupingtable groups : comparisonstable
 * a tabular file with the ENSEMBL info (geneid, gene_biotype, some_biotype, gene_name, gc, length) : ensembltable 
   the ensembltable can be created with a python script, and will also be made optional   

## Example
   there is a very small example/test dataset in example, and it has to be
   run from the example folder
   
```
   cd example
   ./rundemo.R
```


## Important
   deseq2analysis::analyseFromConfig the knitdir has to specified correctly
   in a singularity container its most often not getwd() !
   its also better to specify all paths in the config file as absolute paths
  
