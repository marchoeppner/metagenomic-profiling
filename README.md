![](images/ikmb_bfx_logo.png)

# IKMB Metagenomic profiling pipeline

## Overview

This pipelines analyses short reads and identifies the most likely species in the respective sample. 

## Executing the pipeline

This pipeline requires Nextflow >= 0.30.2. All dependencies are provided through Bioconda and pre-built containers. 
Please make sure that conda/miniconda2 or Singularity are available before starting the pipeline. 

Metaphlan requires a reference databases that is *not* included with the Bioconda packages. On RZCluster, these are
available automatically through the included config file. 

To run the pipeline, do:

`nextflow run marchoeppner/metagenomic-profile --reads '/path/to/reads/' ` 




