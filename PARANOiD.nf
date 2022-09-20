#!/usr/bin/env nextflow

import java.nio.file.*
import groovy.io.FileType

params.split_fastq_by = 1000000
params.speed = false

//essential inputs
input_reads = Channel.fromPath( params.reads )			//FASTQ file(s) containing reads
reference = Channel.fromPath( params.reference )		//FASTA file containing reference sequence(s)
barcode_file = Channel.fromPath( params.barcodes )		//TSV file containing experiment names and the corresponding experiemental barcode sequence

reference.into { reference_to_mapping; reference_to_extract_transcripts; reference_to_extract_sequences; reference_to_pureCLIP } 

params.barcode_pattern = "NNNNNXXXXXXNNNN" 				//STRING containing barcode pattern -> N = random barcode; X = experimental barcode
val_barcode_pattern = Channel.from( params.barcode_pattern )

params.domain = "pro" 									//STRING decides if bowtie2 or STAR is used -> pro = bowtie2; eu = STAR
params.annotation = 'NO_FILE'							//GFF/GTF file containing the annotation belonging to the reference

params.output = "./output/"								//PATH leading to directory in which output is supposed to be stored
params.min_length = 30									//INT minimum length of reads required after preprocessing
params.min_qual = 20									//INT minimum quality of nucleotide required to retain it
params.min_percent_qual_filter = 90						//INT minimum percentage of nucleotides in a read required to have their quality above params.min_qual to rtain the read
params.barcode_mismatches = 1 							//INT maximum number of mismatches in the experimental barcode sequence to still allocate the read

params.mapq = 2 										//INT minimum mapq-score needed to retain an alignment -> soon to be deprecated
params.imgDir = "$PWD/work/images"						//PATH directory to store singularity images

params.map_to_transcripts = false 						//BOOLEAN decides if top X sequences from the reference are presented in the output
params.number_top_transcripts = 10 						//INT number of top transcripts presented when params.map_to_transcripts = true

params.merge_replicates = false

params.peak_calling = false 							//BOOLEAN decides if peak calling via pureclip takes place after normal processing

// RNA species or regions
params.gene_id = "ID"									//STRING name of gene_id used within the annotation file
params.color_barplot = "#69b3a2"						//STRING color used for barplots
//params.rna_species = 'lnc_RNA,miRNA,mRNA,ncRNA,rRNA,snoRNA,snRNA,tRNA'
params.rna_species = '3_prime_UTR,transcript,5_prime_UTR'	//STRING RNA-{species,regions} used for the distribution of CL-sites
Channel.from(params.rna_species.split(',')).into{ rna_species_to_feature_counts; rna_species_to_distribution }
if(params.annotation != 'NO_FILE'){
	annotation_to_RNA_species_distribution = Channel.fromPath( params.annotation )
}


//parameters for peak distance
params.peak_distance = false
params.percentile = 90									//INT percentile that decides which cl-sites are considered when calculating distances and extracting sequences
params.distance = 50 									//INT maximum distance to check for distances between cl-sites

//params for sequence extraction
params.sequence_extraction = false
params.seq_len = 10											//INT length to both sides of cl-sites from which nucleotides are recovered 
params.sequence_format_txt = false 					//BOOLEAN if false sequence are extracted in txt format; if true sequences are extracted in fasta format

// Check if the speed mode was being used. If so the fastq input file will be split every ${params.split_fastq_by} reads to greatly enhance the preprocessing speed
input_reads.into { input_reads_QC; input_reads_processing }
if (params.speed) {
	input_reads_processing.splitFastq( by: params.split_fastq_by, file:true ).set{input_reads_processing}
} 


sjdbGTFfile = file(params.annotation) 
 
//FastQC v0.11.9
process quality_control {
	
	input:
	file query from input_reads_QC

	output:
	file "${query.baseName}*" into out1, collect_statistics_qc1

	"""
	fastqc ${query} -o .
	"""
}

process adapter_removal {

	input:
	file query from input_reads_processing

	output:
	file "${query}_trimmed.fq" into reads_qualityFilter
	file "${query}_trimming_report.txt" into collect_statistics_adapter

	"""
	trim_galore --cores ${task.cpus} --basename ${query} -o . --length ${params.min_length} ${query} --quality 0
	"""
}

//FASTX Toolkit 0.0.14
process quality_filter {

	input:
	file query from reads_qualityFilter

	output:
	file "${query.baseName}.qual-filter.fastq" into fastq_quality_filter_to_barcode_extraction, fastq_quality_filter_to_quality_control_2
	file 'summary-quality-filter.txt' into collect_statistics_quality_filter


	"""
	fastq_quality_filter -v -q ${params.min_qual} -p ${params.min_percent_qual_filter} -i ${query} -o ${query.baseName}.qual-filter.fastq > summary-quality-filter.txt
	"""
}

//FastQC v0.11.9
process quality_control_2 {
	
	input:
	file query from fastq_quality_filter_to_quality_control_2.flatten().toList()

	output:
	file "quality-control-2*" into qc_2_out, collect_statistics_qc2

	"""
	cat ${query} > quality-control-2.fastq
	fastqc quality-control-2.fastq -o .
	"""
}

process extract_rnd_barcode {

	input:
	file query from fastq_quality_filter_to_barcode_extraction

	output:
	file "${query.baseName}.rndBarcode.fastq" into fastq_barcode_extraction_to_split_exp_barcode
	file("${query.simpleName}.log") into collect_statistics_rnd_barcode_extraction

	"""
	umi_tools extract --stdin ${query} --bc-pattern ${params.barcode_pattern} --log ${query.simpleName}.log --stdout ${query.baseName}.rndBarcode.fastq
	"""
}

process check_barcode_file {

	input:
	file "barcodes" from barcode_file

	output:
	file "checkedBarcodes" into checked_barcodes

	"""
	awk '{gsub(/\\W/,"_",\$1); print}' barcodes > checkedBarcodes
	"""
}

//FASTX Toolkit 0.0.14
process split_exp_barcode {

	input:
	set file(query), file(barcodes) from fastq_barcode_extraction_to_split_exp_barcode.combine(checked_barcodes)

	output:
	file "barcode_split.log" into (collect_statistics_split,log_split_experimental_barcode)
	file "output/*.fastq" into fastq_split_barcode_to_remove_exp_barcode

	script:
	if(params.barcode_mismatches == 0)
		"""
		mkdir output
		cat ${query} | fastx_barcode_splitter.pl --bcfile ${barcodes} --bol --exact --prefix ./output/ --suffix .fastq > log_split.txt
		"""
	else
		"""
		mkdir output
		cat ${query} | fastx_barcode_splitter.pl --bcfile ${barcodes} --bol --mismatches ${params.barcode_mismatches} --prefix ./output/ --suffix .fastq > barcode_split.log
		"""
}

process get_length_exp_barcode {

	input:
	val(pattern) from val_barcode_pattern

	output:
	env(exp_barcode_length) into int_length_exp_barcode

	"""
	exp_barcode_length=\$((\$(echo -n ${pattern} | sed 's/[^Xx]//g' | wc -m) + 1)) 
	"""
}

//FASTX Toolkit 0.0.14
process remove_exp_barcode {
	tag {query.simpleName}

	input:
	set file(query), val(length_exp_barcode) from fastq_split_barcode_to_remove_exp_barcode.flatten().filter{ it.size() > 0 }.combine(int_length_exp_barcode)

	output:
	file "*.preprocessed.fastq" into fastq_remove_barcode_to_collect

	"""
	fastx_trimmer -f ${length_exp_barcode} -i ${query} -o ${query.baseName}.preprocessed.fastq
	"""
}

fastq_remove_barcode_to_collect
	.map{file -> tuple(file.name - ~/\.[\w.]+.fastq$/,file)}
	.groupTuple()
	.set{ fastq_collect_preprocessed_to_merge }

process merge_preprocessed_reads {
	tag {name}

	input:
	set val(name), file("query") from fastq_collect_preprocessed_to_merge

	output:
	file("${name}.fastq") into fastq_merge_preprocessed_to_alignment

	"""
	cat ${query} > ${name}.fastq
	"""
}

if ( params.domain == 'pro' || params.map_to_transcripts == true){
	//bowtie2 version 2.3.5.1
	process build_index_bowtie {

		input:
		file ref from reference_to_mapping

		output:
		set file("${ref}"), file("${ref}.*") into bowtie_index_build_to_mapping

		"""
		bowtie2-build ${ref} ${ref}
		"""
	}

	process mapping_bowtie{
		tag {query.simpleName}
		publishDir "${params.output}/alignments", mode: 'copy'

		input:
		set file(ref), file(index) from bowtie_index_build_to_mapping.first()
		file query from fastq_merge_preprocessed_to_alignment

		output:
		file "${query.baseName}.bam" into (bam_mapping_to_filter_empty, bam_mapping_to_output)
		file "${query.simpleName}.statistics.txt" into collect_statistics_mapping

		"""
		bowtie2 --no-unal -q -p ${task.cpus} -U ${query} -x ${ref} 2> ${query.simpleName}.statistics.txt | samtools view -bS - > ${query.baseName}.bam
		"""
	}
} else if ( params.domain == 'eu' ) {
	//Version 2.7.3a
	process build_index_STAR {

		input:
		file referenceGenome from reference_to_mapping
		file gtf from sjdbGTFfile

		output:
		file index into star_index_build_to_mapping

		script:
		if(params.annotation == 'NO_FILE')
			"""
			mkdir index
			STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ./index --genomeFastaFiles ${referenceGenome} 
			"""
		else
			"""
			mkdir index
			STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ./index --genomeFastaFiles ${referenceGenome} --sjdbGTFfile ${gtf}
			"""
	}

	process mapping_STAR{
		tag {query.simpleName}
		publishDir "${params.output}/alignments", mode: 'copy'

		input:
		set file(query), file(indexDir) from fastq_merge_preprocessed_to_alignment.combine(star_index_build_to_mapping)

		output:
		file("${query.baseName}.Aligned.sortedByCoord.out.bam") into (bam_mapping_to_filter_empty, bam_mapping_to_output)
		file("${query.baseName}.Log.*") into collect_statistics_mapping

		"""
		STAR --runThreadN ${task.cpus} --genomeDir ${indexDir} --readFilesIn ${query} --outFileNamePrefix ${query.baseName}. --alignEndsType Extend5pOfRead1 --outSAMtype BAM SortedByCoordinate
		"""
	}
}

process filter_empty_bams{
	tag {query.simpleName}
	echo true

	input:
	file query from bam_mapping_to_filter_empty

	output:
	file "${query.baseName}.filtered.bam" optional true into bam_filter_empty_to_split

	"""
	if [[ \$(samtools view ${query} | wc -l) == 0 ]]; then
		echo "File ${query.simpleName} contains no alignable reads"
	else
		mv ${query} ${query.baseName}.filtered.bam
	fi
	"""
}
if( params.map_to_transcripts == false && params.speed ) {
	process split_bam_by_chromosome{
		tag {query.simpleName}

		input:
		file(query) from bam_filter_empty_to_split

		output:
		file("*.bam") into bam_split_to_sort

		"""
		bamtools split -in ${query} -reference
		"""
	}
} else {
	bam_filter_empty_to_split.set{ bam_split_to_sort }
}

process sort_bam{
	tag {query.simpleName}

	input:
	file query from bam_split_to_sort.flatten()

	output:
	set file("${query.baseName}.sorted.bam"), file("${query.baseName}.sorted.bam.bai")  into bam_sort_to_deduplicate

	"""
	samtools sort ${query} -o ${query.baseName}.sorted.bam
	samtools index ${query.baseName}.sorted.bam
	"""
}

process deduplicate{
	tag {query.simpleName}

	input:
	set file(query), file(index) from bam_sort_to_deduplicate

	output:
	file "${query.baseName}.deduplicated.bam" into (bam_depuplicate_to_sort,bam_deduplicate_to_index)
	file "${query.baseName}.deduplicated.log*" into log_deduplicate_to_collect_statistics

	"""
	umi_tools dedup -I ${query} --output-stats ${query.baseName}.deduplicated.log -S ${query.baseName}.deduplicated.bam
	"""
}

bam_depuplicate_to_sort
	.map{file -> tuple(file.name - ~/\.[\w.]+.bam$/,file)}
	.groupTuple()
	.set{ bam_dedup_sort_to_merge }

process merge_deduplicated_bam {
	tag {name}

	publishDir "${params.output}/alignments/deduplicated", mode: 'copy'

	input:
	set val(name), file(query) from bam_dedup_sort_to_merge

	output:
	file("${name}.bam") into (bam_merge_to_calculate_crosslinks, bam_merge_to_extract_transcripts, bam_merge_to_pureCLIP)

	"""
	samtools merge -c -p ${name}.bam ${query}
	"""
}
if (params.map_to_transcripts == true){
	process count_hits {
		tag {bam.simpleName}

		publishDir "${params.output}/transcripts/overview-hits", mode: 'copy'

		input:
		file bam from bam_merge_to_extract_transcripts

		output:
		file("${bam.baseName}.hits.tsv") into tsv_count_to_get_hits

		"""
		samtools view ${bam} | cut -f3 | sort | uniq -c | sort -nr > ${bam.baseName}.hits.tsv
		"""
	}

	process get_top_hits {

		publishDir "${params.output}/transcripts/overview-hits", mode: 'copy'

		input:
		file tsv from tsv_count_to_get_hits.flatten().toList()

		output:
		file("transcript-targets-top${params.number_top_transcripts}.txt") into (txt_get_hits_to_extract_alignments,txt_get_hits_to_extract_sequences)

		"""
		head -${params.number_top_transcripts} -q *.tsv  | rev | cut -f1 -d' ' | rev | sort | uniq > transcript-targets-top${params.number_top_transcripts}.txt
		"""
	}

	process index_alignments {
		tag {bam.simpleName}

		input:
		file bam from bam_deduplicate_to_index

		output:
		set file("${bam.baseName}.sorted.bam"), file("${bam.baseName}.sorted.bam.bai") into bam_index_to_extract_alignments

		"""
		samtools sort ${bam} -o ${bam.baseName}.sorted.bam
		samtools index ${bam.baseName}.sorted.bam
		"""
	}

	process extract_top_alignments {
		tag {bam.simpleName}

		input:
		set file(bam), file(bai), file(txt_alignments) from bam_index_to_extract_alignments.combine(txt_get_hits_to_extract_alignments)

		output:
		file("${bam.simpleName}_filtered_top${params.number_top_transcripts}.bam") into bam_extract_alignments_to_calc_crosslink

		"""
		samtools view -hb ${bam} `cat ${txt_alignments}` > ${bam.simpleName}_filtered_top${params.number_top_transcripts}.bam
		"""
	}

	process remove_newlines {

		input:
		file ref from reference_to_extract_transcripts

		output:
		file("${ref.baseName}.removed_newlines.fna") into fasta_rm_newline_to_extract_sequences

		"""
		awk '!/^>/ {printf "%s", \$0; n = "\\n" } /^>/ { print n \$0; n = "" } END { printf "%s", n }' ${ref} > ${ref.baseName}.removed_newlines.fna
		"""
	}

	process extract_top_transcript_sequences {
		tag {txt_sequences.simpleName}

		publishDir "${params.output}/transcripts", mode: 'copy'

		input:
		set file(txt_sequences), file(ref) from txt_get_hits_to_extract_sequences.combine(fasta_rm_newline_to_extract_sequences)

		output:
		file("${ref.simpleName}.top${params.number_top_transcripts}_transcripts.fna") into output_top_transcripts


		"""
		egrep -A1 -f ${txt_sequences} ${ref} > ${ref.simpleName}.top${params.number_top_transcripts}_transcripts.fna
		"""
	}

	bam_merge_to_calculate_crosslinks.mix(bam_extract_alignments_to_calc_crosslink)
	.set{collected_bam_files}
} else {
	bam_merge_to_calculate_crosslinks
	.set{collected_bam_files}
}

process calculate_crosslink_sites{
	tag {query.simpleName}
	publishDir "${params.output}/raw-wig-files", mode: 'copy', pattern: "${query.simpleName}_{forward,reverse}.wig"

	input:
	file query from collected_bam_files

	output:
	file "${query.simpleName}.wig2" into wig_calculate_crosslink_to_group_samples
	file "${query.simpleName}_{forward,reverse}.wig" into wig_to_output

	"""
	create-wig-from-bam.py --input ${query} --mapq ${params.mapq} --output ${query.simpleName}.wig2
	wig2-to-wig.py --input ${query.simpleName}.wig2 --output ${query.simpleName}
	"""
}

//From here on only further analyses

if( params.merge_replicates == true ){

	//groups files according to their experiment
	wig_calculate_crosslink_to_group_samples
	.map{file -> tuple(file.name - ~/_rep_\d*(_filtered_top)?\d*.wig2$/,file)} 
	.groupTuple()
	.set{grouped_samples}

	process merge_wigs{
		tag {name}

		publishDir "${params.output}/merged-wig-files", mode: 'copy', pattern: "${name}_{forward,reverse}.wig"

		input:
		set name, file(query) from grouped_samples

		output:
		file "${name}.wig2" into collected_wig_files
		file "${name}_{forward,reverse}.wig" into output

		"""
		merge-wig.py --wig ${query} --output ${name}.wig2
		wig2-to-wig.py --input ${name}.wig2 --output ${name}
		"""
	}
} else {
	wig_calculate_crosslink_to_group_samples
		.set{collected_wig_files}
}

// Generate one channel per postprocessing analysis
collected_wig_files.into{ collected_wig_2_to_RNA_species_distribution; collected_wig_2_to_sequence_extraction; collected_wig_2_to_peak_distance }


if (params.peak_calling == true){
	process index_for_peak_calling {
		tag{query.simpleName}

		input:
		file(query) from bam_merge_to_pureCLIP

		output:
		set file("${query.simpleName}.sorted.bam"), file("${query.simpleName}.sorted.bam.bai") into bambai_index_to_peak_calling

		"""
		samtools sort ${query} -o ${query.simpleName}.sorted.bam
		samtools index ${query.simpleName}.sorted.bam
		"""
	}

	//TODO: think about what to do with mtc. It's only a temporary solution
	process pureCLIP {
		tag{bam.simpleName}
		errorStrategy 'ignore' //TODO: is supposed to be only temporal. Need to find a solution for: ERROR: Emission probability became 0.0! This might be due to artifacts or outliers.

		publishDir "${params.output}/peak_calling", mode: 'copy', pattern: "${bam.simpleName}.pureCLIP_crosslink_sites.bed"

		input:
		set file(bam), file(bai), file(ref) from bambai_index_to_peak_calling.combine(reference_to_pureCLIP)

		output:
		file("${bam.simpleName}.pureCLIP_crosslink_sites.bed")
		file("${bam.simpleName}.pureCLIP_crosslink_sites.params") into params_peak_calling_to_collect_statistics

		"""
		pureclip -i ${bam} -bai ${bai} -g ${ref} -nt ${task.cpus} -o ${bam.simpleName}.pureCLIP_crosslink_sites.bed -mtc 5000 -mtc2 5000 -ld
		"""
	}
}

if (/*params.rna_species == true &&*/ params.annotation != 'NO_FILE'){
	process wig_to_bam {
		tag{query.simpleName}

		input:
		file(query) from collected_wig_2_to_RNA_species_distribution

		output:
		file("${query.baseName}.bam") into bam_convert_to_feature_counts

		"""
		wig-to-bam.py --input ${query} --output ${query.baseName}.bam
		"""

	}

	process feature_counts {

		input:
		set file(query), val(rna_species), file(annotation) from bam_convert_to_feature_counts.combine(rna_species_to_feature_counts).combine(annotation_to_RNA_species_distribution)

		output:
		file("${query.simpleName}.${rna_species}.tsv") into tsv_feature_counts_to_sort


		"""
		featureCounts -T ${task.cpus} -t ${rna_species} -g ${params.gene_id} -a ${annotation} -R CORE -M -o ${query.simpleName} ${query}
		mv ${query}.featureCounts ${query.simpleName}.${rna_species}.tsv
		"""
	}

	tsv_feature_counts_to_sort
		.map{file -> tuple(file.name - ~/\.[\w.]+.tsv$/,file)}
		.groupTuple()
		.set{tsv_sort_to_calculate_distribution}

	process get_RNA_species_distribution {
		echo true

		publishDir "${params.output}/RNA_subtypes", mode: 'copy'

		input:
		val(species) from rna_species_to_distribution.collect()
		set val(name), file(query) from tsv_sort_to_calculate_distribution

		output:
		file("${name}.species_distribution.tsv") into (output_rna_species_tsv,tsv_rna_species_distribution_to_barplot)

		script:
		species_as_string = species.join(' ')
		"""
		calc-RNA-species-distribution.py --input ${query} --rna_species ${species_as_string} --output ${name}.species_distribution.tsv
		"""
	}

	process generate_RNA_species_barplot {
		publishDir "${params.output}/RNA_subtypes", mode: 'copy'

		input:
		file(query) from tsv_rna_species_distribution_to_barplot

		output:
		file("${query.baseName}.png") into output_rna_species_png

		"""
		RNA_species_barcharts.R --input ${query} --output ${query.baseName}.png --color "${params.color_barplot}"
		"""
	}
}
if (params.sequence_extraction == true) {
	process sequence_extraction {

		publishDir "${params.output}/extracted_sequences", mode: 'copy', pattern: "*.extracted-sequences.*"

		input:
		set file(query),file(reference) from collected_wig_2_to_sequence_extraction.combine(reference_to_extract_sequences)

		output:
		file "*.extracted-sequences.*" into extracted_sequences_to_output
		
		script:
		if(params.sequence_format_txt == true)
			"""
			wig2-to-wig.py --input ${query} --output ${query.baseName}
			extract-sequences-around-cross-link-sites.py --input ${query.baseName}*.wig --reference ${reference} --output ${query.baseName}.extracted-sequences.txt --length ${params.seq_len} --percentile ${params.percentile} 
			"""
		else
			"""
			wig2-to-wig.py --input ${query} --output ${query.baseName}
			extract-sequences-around-cross-link-sites.py --input ${query.baseName}*.wig --reference ${reference} --output ${query.baseName}.extracted-sequences.fasta --length ${params.seq_len} --percentile ${params.percentile} --outfmt_fasta
			"""
	}
}

if (params.peak_distance == true) {
	process calculate_peak_distance {

		publishDir "${params.output}/peak_distance", mode: 'copy', pattern: "${query.baseName}.peak-distance.{tsv,png}"

		input:
		file query from collected_wig_2_to_peak_distance

		output:
		file "${query.baseName}.peak-distance.{tsv,png}" into tsv_peak_distance_into_output

		"""
		wig2-to-wig.py --input ${query} --output ${query.baseName}
		peak-distance.py --input ${query.baseName}_{forward,reverse}.wig --output ${query.baseName}.peak-distances.tsv --percentile ${params.percentile} --distance ${params.distance}
		plot-distances.R --input ${query.baseName}.peak-distances.tsv --output ${query.baseName}.peak-distance.png
		"""
	}
}

process generate_barcode_barplot {
	publishDir "${params.output}/statistics", mode: 'copy'

	input:
	file(query) from log_split_experimental_barcode.first()

	output:
	file("${query.baseName}.png") into output_experimental_barcode_distribution_png

	"""
	plot_experimental_barcode_distribution.R --logs ${query} --output ${query.baseName} --type png --color "${params.color_barplot}"
	"""
}


/*
process multiqc{
	publishDir "${params.output}/statistics", mode: 'move'

	input:
	file adapter from collect_statistics_adapter.first().flatten().toList()
	file qual from collect_statistics_quality_filter.first().flatten().toList()
	file qc1 from collect_statistics_qc1.first().flatten().toList()
	file qc2 from collect_statistics_qc2.first().flatten().toList()
	file split from collect_statistics_split.first().flatten().toList()
	file mapping from collect_statistics_mapping.first().flatten().toList()
	file deduplication from log_deduplicate_to_collect_statistics.first().flatten().toList()
	//file(params) from params_peak_calling_to_collect_statistics.flatten().toList().first


	output:
	file "multiqc_*" into multiqc_to_output


	"""
	multiqc .
	"""
}
*/