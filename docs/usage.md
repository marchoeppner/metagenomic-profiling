![](../images/ikmb_bfx_logo.png)

# Usage

This pipeline requires Nextflow >= 19.01.2. All dependencies are provided through Bioconda and pre-built containers.
Please make sure that conda/miniconda2 or Singularity are available before starting the pipeline.

Metaphlan3 requires a reference databases that is *not* included with the Bioconda packages. On the Kiel MedCluster, these are
available automatically through the included config file. 

To run the pipeline, do:

`nextflow run marchoeppner/metagenomic-profile --reads '/path/to/reads/*_R{1,2}_001.fastq.gz'`

## Other options

### `--figures`
Produce graphical representations of the data (heatmap using hclust2 and taxonomic representation with GraphLan2).

### `--genome`
This options needs to be pre-configured! On the IKMB Medcluster, valid options are:

* human
* mouse
* chimp

Will align reads against the respective genome to allow users to gauge the degree of "host" contribution to the overall data. 


