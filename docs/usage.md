![](../images/ikmb_bfx_logo.png)

# Usage

This pipeline requires Nextflow >= 19.01.2. All dependencies are provided through Bioconda and pre-built containers.

On Medcluster, please do the following to load the required dependencies

```bash
module load singularity nextflow
```

Metaphlan3 requires a reference databases that is *not* included with the Bioconda packages. On the Kiel MedCluster, these are
available automatically through the included config file. 

To run the pipeline, do:

```bash
nextflow run marchoeppner/metagenomic-profiling --reads '/path/to/reads/*_R{1,2}_001.fastq.gz'
```

Note, that the wildcard patterns are used to group paired-end reads together. If the naming of your files does not comply with this example, you may have to use a different wildcard pattern. For Illumina data generated at the CCGA, the example shown should work. 

## Other options

### `--figures`
Produce graphical representations of the data (heatmap using hclust2 and taxonomic representation with GraphLan2).

### `--genome`
This options needs to be pre-configured! On the IKMB Medcluster, valid options are:

* human
* mouse
* chimp

Will align reads against the respective genome to allow users to gauge the degree of "host" contribution to the overall data. 

### `--virus`
Also check for viruses not included in Metaphlan. This uses Kraken2 and the kraken2_db variable to point to a compatible database. For the Medcluster, this database was built from RefSeq viruses. 

### `--rapid`
This option performs a pre-filtering of reads against a pre-defined bloom filter to remove potential host reads. This option is useful if you are trying to analyze very deep metagenomes with potentially a lot of host contamination (dozens of GB in raw data). For typical metagenomes (2-6GB), this option is not recommended. 
 

