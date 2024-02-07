/*
 * Extracts random barcodes from reads and puts them into the header
 */
process extract_rnd_barcode {
	tag {query.simpleName}

	input:
	path(query)

	output:
	path("${query.baseName}.rndBarcode.fastq"), emit: fastq_rnd_barcode_extracted
	path("${query.simpleName}.log"), 			emit: report_rnd_barcode_extraction
	path("${task.process}.version.txt"), 		emit: version

	"""
	umi_tools extract \
		--stdin ${query} \
		--bc-pattern ${params.barcode_pattern} \
		--log ${query.simpleName}.log \
		--stdout ${query.baseName}.rndBarcode.fastq

	echo -e "${task.process}\tumi-tools\t\$(umi_tools --version | cut -f3 -d' ')" > ${task.process}.version.txt
	"""
}

/*
 * Checks barcode file on correctness: Throws error when -> Duplicate experiment names, duplicate barcodes, non-nucleotides in barcode, invalid whitespaces, wrong number of columns, incorrect number of input files
 */
process check_barcode_file {

	input:
	path("barcodes")

	output:
	path("checkedBarcodes"), 				emit: barcodes
	path("${task.process}.version.txt"),	emit: version

	"""
	check_barcode_file.py barcodes > checkedBarcodes

	echo -e "${task.process}\tcustom-script\t1.0.0" > ${task.process}.version.txt
	"""
}

/*
 * Splits reads according to the experimental barcode. Generates one read file for each experiment plus one for unmatched reads
 */
process split_exp_barcode {
	tag {query.simpleName}
	publishDir "${params.output}/statistics/barcode-distribution", mode: 'copy', pattern: 'barcode_split.txt'

	input:
	tuple path(query), path(barcodes)

	output:
    path "output/*.fastq", 					emit: fastq_split_experimental_barcode
	path "barcode_split.log", 				emit: report_split_experimental_barcode 
	path("${task.process}.version.txt"),	emit: version

	script:
	def local_mismatches = params.barcode_mismatches == 0 ? '--exact' : '--mismatches ' + params.barcode_mismatches

	"""
	mkdir output

	cat ${query} \
		| fastx_barcode_splitter.pl \
		--bcfile ${barcodes} \
		--bol \
		${local_mismatches} \
		--prefix ./output/ \
		--suffix .fastq \
		> barcode_split.log

	echo -e "${task.process}\tcustom-script\tno_version_available" > ${task.process}.version.txt
	"""
}

/*
 * Extracts length of experimental barcode
 */
process get_length_exp_barcode {

	input:
	val(pattern)

	output:
	env(exp_barcode_length), emit: var_length

	"""
	exp_barcode_length=\$((\$(echo -n ${pattern} | sed 's/[^Xx]//g' | wc -m) + 1)) 
	"""
}

/*
 * Removes experimental barcode from 5' end of reads
 */
process remove_exp_barcode {
	tag {query.simpleName}

	input:
	tuple path(query), val(length_exp_barcode)

	output:
	path("*.preprocessed.fastq"), 			emit: fastq_trimmed
	path("${task.process}.version.txt"),	emit: version

	"""
	fastx_trimmer \
		-f ${length_exp_barcode} \
		-i ${query} \
		-o ${query.baseName}.preprocessed.fastq

	echo -e "${task.process}\tfastx_trimmer\t\$(fastx_trimmer -h | head -2 | tail -1 | rev | cut -d ' ' -f 5- | rev)" > ${task.process}.version.txt
	"""
}

/*
 * Merges several FASTQ files belonging to the same experiment into a single file
 */
process merge_preprocessed_reads {
	tag {name}

	input:
	tuple val(name), path("query")

	output:
	path("${name}.fastq"), emit: fastq_merged

	"""
	cat ${query} > ${name}.fastq
	"""
}

/*
 * Generates a plot showing the ampunt of reads being assigned to their experiments
 */
process generate_barcode_barplot {
	publishDir "${params.output}/statistics/barcode-distribution", mode: 'copy'

	input:
	path(query)

	output:
	path("${query.baseName}.png"), 			emit: png_to_output_dir
	path("${task.process}.version.txt"),	emit: version

	"""
	plot_experimental_barcode_distribution.R \
		--logs ${query} \
		--output ${query.baseName} \
		--type png \
		--color "${params.color_barplot}"

	echo -e "${task.process}\tplot_experimental_barcode_distribution.R\tcustom_script" > ${task.process}.version.txt
	"""
}
