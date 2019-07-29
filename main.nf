// Pipeline variables

OUTDIR = params.outdir 

METAPHLAN_PKL=file(params.metaphlan_pkl)
METAPHLAN_DB=params.metaphlan_db
REF = file(params.ref)

inputFile=file(params.samples)

// Logging and reporting

params.version = "1.0" 
// Header log info 

log.info "=========================================" 
log.info "IKMB pipeline version v${params.version}" 
log.info "Nextflow Version: $workflow.nextflow.version" 
log.info "Command Line: $workflow.commandLine" 
log.info "=========================================" 


// Starting the workflow

Channel.fromFilePairs(params.reads)
	.ifEmpty {exit 1; "Could not find the specified input reads $params.reads"}
	.set { inputFastp }

process runFastp {

	tag "${indivID}|${sampleID}"

	input:
	set val(sampleID), fastqR1, fastqR2 from inputFastp

	output:
	set val(sampleID), file(left),file(right) into fastpOutput
	set file(html),file(json) into fastp_results

	script:
	left = file(fastqR1).getBaseName() + "_trimmed.fastq.gz"
	right = file(fastqR2).getBaseName() + "_trimmed.fastq.gz"
	json = file(fastqR1).getBaseName() + ".fastp.json"
	html = file(fastqR1).getBaseName() + ".fastp.html"

	"""
		fastp --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right -w ${task.cpus} -j $json -h $html --length_required 35
	"""
}

process runBwa {

   tag "${patientID}|${sampleID}"
   publishDir "${OUTDIR}/${patientID}/${sampleID}/Host", mode: 'copy'

   input:
   set sampleID,file(left),file(right) from inputBwa

   output:
   set sampleID,file(bam) into alignedBam
   file(stats) into BamStats

   file(samtools_version) into version_samtools

   script:

   bam = sampleID + ".bam"
   stats = sampleID + "_bwa_stats.txt"

   samtools_version = "v_samtools.txt"

   """
        samtools --version &> $samtools_version
	bwa mem -M -t ${task.cpus} ${REF} $left $right | /opt/samtools/1.9/bin/samtools sort -O BAM - > $bam
	samtools stats $bam > $stats
	
   """

}

inputMerge = alignedBam.groupTuple(by: [0,1])

process runMergeBam {

   input:
   set sampleID,file(bams) from inputMerge

   output:
   set sampleID,file(merged_bam) into mergedBam

   script:
   merged_bam = patientID + "_" + sampleID + ".merged.bam"
   
   """
	samtools merge $merged_bam $bams
   """

}

// We extract the reads not mapping to the host genome
process extractUnmapped {

   tag "${patientID}|${sampleID}"
   publishDir "${OUTDIR}/${patientID}/${sampleID}/Host", mode: 'copy'

   input:
   set sampleID,file(bam) from mergedBam

   output:
   set sampleID,file(left),file(right) into inputMetaphlan
  
   script:
   left = sampleID + "_R1.fastq.gz"
   right = sampleID + "_R2.fastq.gz"

   """
	samtools fastq -f 4 -1 $left -2 $right $bam
   """

}

process runMetaphlan {

   tag "${patientID}|${sampleID}"
   publishDir "${OUTDIR}/${patientID}/${sampleID}/Metaphlan2", mode: 'copy'

   input:
   set sampleID,file(left_reads),file(right_reads) from inputMetaphlan

   output:
   file(metaphlan_out) into outputMetaphlan
   file "v_metaphlan.txt" into version_metaphlan

   script:

   metaphlan_out = sampleID + "_metaphlan_report.txt"

   """
     metaphlan2.py --version &> v_metaphlan.txt
     metaphlan2.py --bowtie2db $METAPHLAN_DB --nproc ${task.cpus} --input_type fastq <(zcat $left_reads $right_reads ) > $metaphlan_out

   """

}

process runMergeAbundance {

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
    file('*') from trimgalore_fastqc_reports.flatten().toList()

    output:
    file("fastq_multiqc*") into runMultiQCFastqOutput

    script:

    """
    multiqc -n fastq_multiqc *.zip *.html
    """
}

workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="
}
