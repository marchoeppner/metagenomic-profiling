// Pipeline variables

OUTDIR = params.outdir 

def helpMessage() {
  log.info"""
  =================================================================
   IKMB | Metaphlan3 Pipeline | v${workflow.manifest.version}
  =================================================================
  Usage:
  The typical command for running the pipeline is as follows:
  nextflow run marchoeppner/metagenomic-profile --reads '/path/to/*_R{1,2}_001.fastq.gz' 
  Mandatory arguments:
  --reads 		The path to the folder containing PE metagenomic reads (1 per sample)

  Optonal arguments:
  --virus		Run Kraken on virus database
  --ref         	The path to the reference genome for mapping stats (must be accompanied by a BWA index!)
  --genome		Instead of --ref, use a pre-configured genome sequence by its common name (only RZCluster)
  --email 		An eMail adress to which reports are sent
  --figures 		Create overview graphics from the result (default: false). Only recommended for smaller sample sizes. 
  -profile      	The nextflow execution profile to use

  """.stripIndent()
}

// Show help message
if (params.help){
	helpMessage()
	exit 0
}

// sanity checks 

if (params.metaphlan_db) {
	METAPHLAN_DB=params.metaphlan_db
	db_path = file(METAPHLAN_DB)
	if (!db_path.exists()) exit 1, "Could not find your Metaphlan DB - please check the path"
} else {
	exit 1, "No Metaphlan database was specified, aborting..."
}
if (params.virus && params.kraken2_db) {
	KRAKEN2_DB=params.kraken2_db
	db_path = file(KRAKEN2_DB)
	if (!db_path.exists()) exit 1, "Could not find your KrakenDB - please check the path"
	params.kraken = true
} else if (params.virus) {
        exit 1, "No Kraken database was specified, aborting..."
}

if (params.genome) {
	if (!params.genomes) {
		exit 1, "Specified a genome name for host mapping, but no genomes are configured for your profile...exiting."
	} else if (!params.genomes.containsKey(params.genome)) {
		exit 1, "Specified unknown name for the host genome...valid options are: ${params.genomes.keySet()}"
	}
}

if (!params.kneaddata_db) {
	exit 1, "No Kneadata DB directory defined (--kneaddata_db)"
}

if (params.ref) {
	REF = file(params.ref)
	index_file = file(params.ref + ".bwt")
	if (!REF.exists()) exit 1, "Could not find the specified reference genome - please check the path"
	if  (!index_file.exists()) exit 1, "Found genome reference, but seems to be missing the BWA index files"
} else if (params.genome) {
	REF = file(params.genomes[ params.genome])
        if (!REF.exists()) exit 1, "Could not find the specified reference genome - please check the path"
}

// Logging and reporting

params.version = workflow.manifest.version

// Header log info 

log.info "=========================================" 
log.info "METAPHLAN3 P I P E L I N E"
log.info "IKMB pipeline version v${params.version}" 
log.info "Nextflow Version: 	$workflow.nextflow.version" 
log.info "=== Inputs =============================="
log.info "Kneaddata DB:		${params.kneaddata_db}"
log.info "Metaphlan DB:		${params.metaphlan_db}"
log.info "Reads:			${params.reads}"
if (params.genome) {
	log.info "Host genome: 		${params.genome}"
} else if (params.ref) {
	log.info "Host genome: 		${params.ref}"
}
log.info "=========================================="
log.info "Command Line:         $workflow.commandLine"
if (workflow.containerEngine) {
	log.info "Container Engine: 	${workflow.containerEngine}"
}
log.info "=========================================" 

// Starting the workflow

Channel.fromFilePairs(params.reads )
	.ifEmpty {exit 1, "Could not find the specified input reads $params.reads"}
	.into { inputKneaddata ; inputFastQC ; inputBwa  }

process runFastQC {

	label 'fastqc'

	publishDir "${OUTDIR}/${sampleID}/FastQC", mode: 'copy'

	input:
	set val(sampleID),file(reads) from inputFastQC

	output:
	file "*_fastqc.{zip,html}" into fastqc_results

	script:
	"""
	fastqc --quiet --threads $task.cpus $reads
	"""

}

process runKneaddata {

	label 'kneaddata'

        publishDir "${OUTDIR}/${sampleID}/Kneaddata", mode: 'copy'

        //scratch true

        input:
        set val(sampleID),file(reads) from inputKneaddata

        output:
        set val(sampleID),file("${outdir}/${left}"),file("${outdir}/${right}") into (inputMetaphlan,inputKraken)

        script:
        left = sampleID + "_R1_001_kneaddata_paired_1.fastq.gz"
        right = sampleID + "_R1_001_kneaddata_paired_2.fastq.gz"

        outdir = "output"

        """
                kneaddata --input ${reads[0]} --input ${reads[1]} \
                        -t ${task.cpus} \
                        --reference-db ${params.kneaddata_db} \
			--run-trf \
                        --output $outdir \
                        --trimmomatic /usr/local/share/trimmomatic-0.39-1

                cd $outdir
                for i in \$(echo *paired_*.fastq); do gzip \$i ; done;
        """

}

if ( params.ref || params.genome ) {
	process runBwa {

	   publishDir "${OUTDIR}/${sampleID}/Host", mode: 'copy'

	   scratch true

	   input:
	   set sampleID,file(reads) from inputBwa

	   output:
	   file(stats) into BamStats

	   file(samtools_version) into version_samtools

	   script:

	   bam = sampleID + ".host_mapped.bam"
	   stats = sampleID + ".txt"

	   samtools_version = "v_samtools.txt"

	   """
        	samtools --version &> $samtools_version
		bwa mem -M -t ${task.cpus} ${REF} $reads | samtools sort -m 4G -O BAM - > $bam
		samtools stats $bam > $stats
		rm $bam
	
	   """
	}

} else {
	BamStats = Channel.empty()
}

if (params.kraken) {
	process runKraken2 {

        	label 'kraken'

	        publishDir "${OUTDIR}/${sampleID}/Kraken/", mode: 'copy'

	        input:
        	set val(sampleID),file(left),file(right) from inputKraken

	        output:
        	set val(sampleID),file(report) into KrakenReport
	        file(kraken_log)

        	script:
	        report = sampleID + ".kraken2_report.txt"
        	kraken_log = sampleID + ".kraken2.log"

	        """
        	        kraken2 --db $KRAKEN2_DB --threads ${task.cpus} --output $kraken_log --report $report $left $right
	        """
	}

	process Kraken2Yaml {

	        input:
        	file(reports) from KrakenReport.collect()

	        output:
        	file(report_yaml) into KrakenYaml

	        script:

        	report_yaml = "kraken_report_mqc.yaml"
	        """	
        	        kraken2yaml.pl --outfile $report_yaml
	        """

	}

} else {
	KrakenYaml = Channel.empty()
}

process runMetaphlan {

   label 'metaphlan'

   publishDir "${OUTDIR}/${sampleID}/Metaphlan3", mode: 'copy'

   input:
   set val(sampleID),file(left_reads),file(right_reads) from inputMetaphlan

   output:
   file(metaphlan_out) into outputMetaphlan
   set val(sampleID),file(sam_out) into outputMetaphlanBowtie
   file "v_metaphlan.txt" into version_metaphlan

   script:

   metaphlan_out = sampleID + ".out"
   bowtie_out = sampleID + "bowtie2.txt"
   sam_out = sampleID + ".sam.bz2"
   """
     metaphlan --version &> v_metaphlan.txt
     zcat $left_reads > left.fq
     zcat $right_reads > right.fq
     metaphlan left.fq,right.fq --bowtie2db $METAPHLAN_DB --samout $sam_out --bowtie2out $bowtie_out --nproc ${task.cpus} -o $metaphlan_out --input_type fastq
     rm *.fq

   """

}

process runSample2Markers {

	label 'metaphlan'

	publishDir "${OUTDIR}/Strainphlan2/Markers", mode: 'copy'

	when:
	params.strainphlan

	input:
	set val(sampleID),file(sam_out) from outputMetaphlanBowtie

	output:
	set val(sampleID),file("*.markers") into SampleMarkers

	script:

	"""
		sample2markers.py -n ${task.cpus} -i $sam_out -o .
	"""
}

process runMergeAbundance {

	publishDir "${OUTDIR}/Metaphlan3", mode: 'copy'

	input:
	file(results) from outputMetaphlan.collect()

	output:
	file(abundances) into (inputHeatmap, inputGraphlan)

	script:
	abundances = "metaphlan_abundances.txt"

	"""
		merge_metaphlan_tables.py ${results.join(" ")} > $abundances
	"""
}

process runGraphlan {

	publishDir "${OUTDIR}/Metaphlan3/Phylogeny", mode: 'copy'

	when:
	params.figures
	
	input:
	file(abundances) from inputGraphlan

	output:
	set file(pyhlo_png),file(phylo_svg),file(phylo_xml),file(phylo_annot),file(phylo_legend),file(phylo_svg_legend),file(phylo_svg_annot) into outputGraphlan

	script:
	pyhlo_png = "metaphlan.phylogeny.png"
        phylo_svg = "metaphlan.phylogeny.svg"
	phylo_svg_annot = "metaphlan.phylogeny_annot.svg"
	phylo_svg_legend = "metaphlan.phylogeny_legend.svg"
	phylo_xml = "metaphlan.phylogeny.xml"
	phylo_annot = "metaphlan.phylogeny_annot.png"
	phylo_legend = "metaphlan.phylogeny_legend.png"

	"""
		export2graphlan.py --skip_rows 1,2 \
			-i $abundances \
			--tree merged_abundance.tree.txt \
			--annotation merged_abundance.annot.txt \
			--most_abundant 100 \
			--abundance_threshold 1 \
			--least_biomarkers 10 \
			--annotations 5,6 \
			--external_annotations 7 \
			--min_clade_size 1

		graphlan_annotate.py --annot merged_abundance.annot.txt merged_abundance.tree.txt $phylo_xml

		graphlan.py --dpi ${params.dpi}  $phylo_xml $pyhlo_png --external_legends
		graphlan.py --dpi ${params.dpi}  $phylo_xml $phylo_svg --external_legends
	"""

}

process runBuildHeatmap {

	publishDir "${OUTDIR}/Metaphlan3/Heatmap", mode: 'copy'

	when:
	params.figures 

	input:
	file(abundance) from inputHeatmap

	output:
	file(heatmap) into outputHeatmap

	script:
	heatmap = "metaphlan.abundances.png"
	"""
		grep -E "(s__)|(^ID)" $abundance | grep -v "t__" | sed 's/^.*s__//g' > merged_abundance_table_species.txt
		hclust2.py -i merged_abundance_table_species.txt -o $heatmap \
		--ftop 25 \
		--f_dist_f braycurtis \
		--s_dist_f braycurtis \
		--cell_aspect_ratio 0.5 \
		-l --flabel_size 6 \
		--slabel_size 6 --max_flabel_len 100 \
		--max_slabel_len 100 \
		--minv 0.1 --dpi 300
	"""
}

process runMultiQC {

	publishDir "${OUTDIR}/MultiQC", mode: 'copy'

	label 'multiqc'

	input:
	file ('*') from fastqc_results.collect()
	file ('*') from BamStats.collect()
	file('*') from KrakenYaml.collect()
	output:
	file("multiqc_report.html")

	script:

	"""
		cp $params.logo .
                cp $baseDir/assets/multiqc_config.yaml multiqc_config.yaml

		multiqc .
	"""
}

workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="
}
