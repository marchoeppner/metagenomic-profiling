// Pipeline variables

OUTDIR = params.outdir 

def helpMessage() {
  log.info"""
  =================================================================
   IKMB | Metaphlan2 Pipeline | v${workflow.manifest.version}
  =================================================================
  Usage:
  The typical command for running the pipeline is as follows:
  nextflow run marchoeppner/metagenomic-profile --reads '/path/to/*_R{1,2}_001.fastq.gz' 
  Mandatory arguments:
  --reads 	The path to the folder containing PE metagenomic reads (1 per sample)

  Optonal arguments:
  --ref         The path to the reference genome for mapping stats (must be accompanied by a BWA index!)
  --genome	Instead of --ref, use a pre-configured genome sequence by its common name (only RZCluster)
  --email 	An eMail adress to which reports are sent
  -profile      The nextflow execution profile to use

  """.stripIndent()
}

// Show help message
if (params.help){
	helpMessage()
	exit 0
}

// sanity checks 

if (params.metaphlan_pkl) {
	METAPHLAN_PKL=file(params.metaphlan_pkl)
	if (!METAPHLAN_PKL.exists()) exit 1, "Could not find the Metaphlan PKL file - please check the path"
} else {
	exit 1, "No Metaphlan PKL file was specified, aborting..."
}

if (params.metaphlan_db) {
	METAPHLAN_DB=params.metaphlan_db
	db_path = file(METAPHLAN_DB)
	if (!db_path.exists()) exit 1, "Could not find your Metaphlan DB - please check the path"
} else {
	exit 1, "No Metaphlan database was specified, aborting..."
}

if (params.genome) {
	if (!params.genomes) {
		exit 1, "Specified a genome name for host mapping, but no genomes are configured for your profile...exiting."
	} else if (!params.genomes.containsKey(params.genome)) {
		exit 1, "Specified unknown name for the host genome...valid options are: ${params.genomes.keySet()}"
	}
}

if (params.ref) {
	REF = file(params.ref)
	index_file = file(params.ref + ".bwt")
	if (!REF.exists()) exit 1, "Could not find the specified reference genome - please check the path"
	if  (!index_file.exists()) exit 1, "Found genome reference, but seems to be missing the BWA index files"
} else if (params.genome) {
	REF = file(params.genomes[ params.genome])
        if (!REF.exists()) exit 1, "Could not find the specified reference genome - please check the path"
} else {
	exit 1, "Must provide the path to a refrence genome and BWA index"
}

// Logging and reporting

params.version = workflow.manifest.version

// Header log info 

log.info "=========================================" 
log.info "METAPHLAN2 P I P E L I N E"
log.info "IKMB pipeline version v${params.version}" 
log.info "Nextflow Version: $workflow.nextflow.version" 
log.info "Command Line: $workflow.commandLine" 
if (params.genome) {
	log.info "Host genome: ${params.genome}"
} else if (params.ref) {
	log.info "Host genome: ${params.ref}"
}
if (workflow.containerEngine) {
	log.info "Container Engine: ${workflow.containerEngine}"
}
log.info "=========================================" 

// Starting the workflow

Channel.fromFilePairs(params.reads)
	.ifEmpty {exit 1, "Could not find the specified input reads $params.reads"}
	.set { inputFastp }

process runFastp {

	publishDir "${OUTDIR}/${sampleID}/fastp"

	input:
	set val(sampleID), file(reads) from inputFastp

	output:
	set val(sampleID), file(left),file(right) into inputBwa,inputMetaphlan
	set file(html),file(json) into fastp_results

	script:
	left = reads[0].getBaseName() + "_trimmed.fastq.gz"
	right = reads[1].getBaseName() + "_trimmed.fastq.gz"
	json = reads[0].getBaseName() + ".fastp.json"
	html = reads[0].getBaseName() + ".fastp.html"

	"""
		fastp --in1 ${reads[0]} --in2 ${reads[1]} --out1 $left --out2 $right -w ${task.cpus} -j $json -h $html --length_required 35
	"""
}

if (params.ref) {
	process runBwa {

	   publishDir "${OUTDIR}/${sampleID}/Host", mode: 'copy'

	   input:
	   set sampleID,file(left),file(right) from inputBwa

	   output:
	   set sampleID,file(bam) into alignedBam
	   file(stats) into BamStats

	   file(samtools_version) into version_samtools

	   script:

	   bam = sampleID + ".host_mapped.bam"
	   stats = sampleID + "_bwa_stats.txt"

	   samtools_version = "v_samtools.txt"

	   """
        	samtools --version &> $samtools_version
		bwa mem -M -t ${task.cpus} ${REF} $left $right | samtools sort -m 4G -O BAM - > $bam
		samtools stats $bam > $stats
	
	   """
	}
}

process runMetaphlan {

   publishDir "${OUTDIR}/${sampleID}/Metaphlan2", mode: 'copy'

   input:
   set sampleID,file(left_reads),file(right_reads) from inputMetaphlan

   output:
   file(metaphlan_out) into outputMetaphlan
   file "v_metaphlan.txt" into version_metaphlan

   script:

   metaphlan_out = sampleID + ".out"

   """
     metaphlan2.py --version &> v_metaphlan.txt
     metaphlan2.py --bowtie2db $METAPHLAN_DB --nproc ${task.cpus} --input_type fastq <(zcat $left_reads $right_reads ) > $metaphlan_out

   """

}

process runMergeAbundance {

	publishDir "${OUTDIR}/Metaphlan2", mode: 'copy'

	input:
	file(results) from outputMetaphlan.collect()

	output:
	file(abundances) into abundanceMetaphlan

	script:
	abundances = "metaphlan_abundances.txt"

	"""
		merge_metaphlan_tables.py ${results.join(" ")} > $abundances
	"""
}

process runBuildHeatmap {

	publishDir "${OUTDIR}/Metaphlan2", mode: 'copy'

	input:
	file(abundance) from abundanceMetaphlan

	output:
	file(heatmap) into outputHeatmap

	script:
	heatmap = "metaphlan.abundances.png"

	"""
		grep -E "(s__)|(^ID)" $abundance | grep -v "t__" | sed 's/^.*s__//g' > merged_abundance_table_species.txt
		hclust2.py -i merged_abundance_table_species.txt -o $heatmap --ftop 25 \
			--f_dist_f braycurtis \
			--s_dist_f braycurtis \
			--cell_aspect_ratio 0.5 \
			-l --flabel_size 6 \
			--slabel_size 6 \
			--max_flabel_len 100 \
			--max_slabel_len 100 \
			--minv 0.1 \
			--dpi 300
	"""
}

process runMultiQCFastq {

    tag "Generating fastq level summary and QC plots"
    publishDir "${OUTDIR}/Summary/Fastqc", mode: 'copy'

    input:
    file('*') from fastp_results.flatten().toList()

    output:
    file("fastq_multiqc*") into runMultiQCFastqOutput

    script:

    """
    multiqc -n fastq_multiqc *.json *.html
    """
}

workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="
}
