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

if (!params.genome) {
	exit 1, "No Host genome was specified...valid options are: ${params.genomes.keySet()}"
} else if (!params.genomes) {
	exit 1, "Specified a genome name for host mapping, but no genomes are configured for your profile...exiting."
} else if (!params.genomes.containsKey(params.genome)) {
	exit 1, "Specified unknown name for the host genome...valid options are: ${params.genomes.keySet()}"
} else {

	log.info "Using ${params.genome} as host species...setting relevant options"

	bowtie_base = params.genomes[params.genome].bowtie_index
	bloom_index = params.genomes[params.genome].bloom_index

	host_genome = Channel.fromPath("${bowtie_base}*")
}

if (params.rapid) {
	log.info "Running in rapid mode - this may produce less accurate results!"
}

BLOOMFILTER_HOST = bloom_index

// Logging and reporting

params.version = workflow.manifest.version

// Header log info 

log.info "=========================================" 
log.info "METAPHLAN3 P I P E L I N E"
log.info "IKMB pipeline version v${params.version}" 
log.info "Nextflow Version: 	$workflow.nextflow.version" 
log.info "=== Inputs =============================="
log.info "Metaphlan DB:		${params.metaphlan_db}"
log.info "Reads:			${params.reads}"
log.info "Host genome: 		${params.genome}"
if (params.rapid) {
	log.info "Host filter:		${params.bloomfilter_host}"
}
log.info "=========================================="
log.info "Command Line:         	$workflow.commandLine"
if (workflow.containerEngine) {
	log.info "Container Engine: 	${workflow.containerEngine}"
}
log.info "=========================================" 

// Starting the workflow

Channel.fromFilePairs(params.reads , flat: true )
	.ifEmpty {exit 1, "Could not find the specified input reads $params.reads"}
	.into { Reads ; inputFastQC }


process runFastQC {

	label 'fastqc'

	publishDir "${OUTDIR}/${sampleID}/FastQC", mode: 'copy'

	input:
	set val(sampleID),file(left),file(right) from inputFastQC

	output:
	file "*_fastqc.{zip,html}" into fastqc_results

	script:
	"""
	fastqc --quiet --threads $task.cpus $left $right
	"""

}

if (params.rapid) {

	process BloomFilter {

		label 'biobloom'

		scratch true

		input:
		set val(sampleID),file(left),file(right) from Reads

		output:
		set val(sampleID),file(clean_reads) into inputReformat

		script:
		analysis = sampleID + ".Host"
		clean_reads = analysis + ".filtered.fastq.gz"

                """
                        biobloomcategorizer -p $analysis --gz_output -d -n -e -s 0.01 -t ${task.cpus} -f "$BLOOMFILTER_HOST" $left $right | gzip > $clean_reads
                """

        }

        // *************************
        // Take the interlaved non-host reads and produce sane PE data
        // *************************
        process runDeinterlave {

                label 'bbmap'

		scratch true

                publishDir "${OUTDIR}/${id}/Bloomfilter/Host", mode: 'copy'

                input:
                set val(id),file(reads) from inputReformat

                output:
                set val(id),file(left),file(right) into inputBBduk

                script:
                left = id + "_R1_001.bloom_non_host.fastq.gz"
                right = id + "_R2_001.bloom_non_host.fastq.gz"

                """
                        reformat.sh in=$reads out1=$left out2=$right addslash int
                """

        }

} else {

	inputBBduk = Reads

}

process trimReads {

	label 'bbmap'

	scratch true

	input:
	set val(sampleID),file(left),file(right) from inputBBduk
	
	output:
	set val(sampleID),file(left_trimmed),file(right_trimmed) into filterPEReads
	set val(sampleID),file(unpaired) into filterSEReads
	path bbduk_adapter_stats
		
	script:
	bbduk_adapter_stats = sampleID + ".bbduk.adapter.stats"

	left_trimmed = left.getBaseName() + "_trimmed.fastq.gz"
	right_trimmed = right.getBaseName() + "_trimmed.fastq.gz"

	unpaired = sampleID + "_unpaired.fastq.gz"

	"""
		bbduk.sh stats=$bbduk_adapter_stats threads=${task.cpus} in=${left} in2=${right} out1=${left_trimmed} out2=${right_trimmed} outs=$unpaired ref=${params.adapters} ktrim=r k=23 mink=11 hdist=1 minlength=${params.min_read_length} tpe tbo
	"""
}

process cleanPEReads {

	label 'bbmap'

	scratch true

	input:
	set val(sampleID),file(left),file(right) from filterPEReads

	output:
	set val(sampleID),file(left_clean),file(right_clean) into pe_reads_clean

	script:
	left_clean = left.getBaseName() + "_clean.fastq.gz"
	right_clean = right.getBaseName() + "_clean.fastq.gz"
	artifact_stats = sampleID + ".bbduk.artifacts.stats"
	
	"""
		 bbduk.sh stats=$artifact_stats threads=${task.cpus} in=${left} in2=${right} k=31 ref=artifacts,phix ordered cardinality out1=${left_clean} out2=${right_clean} minlength=${params.min_read_length}
	"""

}

process cleanSEReads {

	label 'bbmap'

	scratch true

	input:
	set val(sampleID),file(unpaired) from filterSEReads

	output:
	set val(sampleID),file(unpaired_clean) into se_reads_clean

	script:

	unpaired_clean = unpaired.getBaseName() + "_clean.fastq.gz"

	"""
		bbduk.sh threads=${task.cpus} in=${unpaired}  k=31 ref=artifacts,phix ordered cardinality out1=${unpaired_clean} minlength=${params.min_read_length}

	"""

}

reads_for_mapping = pe_reads_clean.join(se_reads_clean)

process filterReads {

	label 'bowtie2'

	publishDir "${params.outdir}/${sampleID}/reads_clean", mode: 'copy'

	input:
	set val(sampleID),file(left),file(right),file(unpaired) from reads_for_mapping
	path(bwt_files) from host_genome.collect()

	output:
	set val(sampleID),path(left_clean),path(right_clean),path(unpaired_clean) into inputKraken,inputMetaphlan
	path(bowtie_log) into bowtie_log

	script:
	left_clean = sampleID + ".clean.R1.fastq.gz"
	right_clean = sampleID + ".clean.R2.fastq.gz"
	unpaired_clean = sampleID + ".clean.unpaired.fastq.gz"
	bowtie_log = sampleID + ".txt"
	index = "genome"
	"""
		bowtie2 -x $index -1 $left -2 $right -U $unpaired -S /dev/null --no-unal -p ${task.cpus} --un-gz $unpaired_clean  --un-conc-gz ${sampleID}.clean.R%.fastq.gz 2> $bowtie_log
	"""

}

if (params.kraken) {
	process runKraken2 {

        	label 'kraken'

	        publishDir "${OUTDIR}/${sampleID}/Kraken/", mode: 'copy'

	        input:
        	set val(sampleID),file(left),file(right),file(unpaired) from inputKraken

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
   set val(sampleID),file(left_reads),file(right_reads),file(unpaired) from inputMetaphlan

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
