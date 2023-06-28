/*
 * Extracts random barcodes from reads and puts them into the header
 * Params: params.barcode_pattern -> Pattern used to extract random barcode
 * Input: [FASTQ] 1 or more read file(s)
 * Output: fastq_rnd_barcode_extracted -> [FASTQ] read file(s) with random barcode in the appropriate header
 *		report_rnd_barcode_extraction -> [LOG] report of random barcode extraction
 */
process extract_rnd_barcode {
	tag {query.simpleName}

	input:
	path(query)

	output:
	path("${query.baseName}.rndBarcode.fastq"), emit: fastq_rnd_barcode_extracted
	path("${query.simpleName}.log"), emit: report_rnd_barcode_extraction

	"""
	umi_tools extract --stdin ${query} --bc-pattern ${params.barcode_pattern} --log ${query.simpleName}.log --stdout ${query.baseName}.rndBarcode.fastq
	"""
}

/*
 * Checks barcode file on correctness: Throws error when -> Duplicate experiment names, duplicate barcodes, non-nucleotides in barcode, invalid whitespaces, wrong number of columns, incorrect number of input files
 * Input: [TSV] 1 barcode file
 * Output: [TSV] Corrected barcode file 
 */
process check_barcode_file {

	input:
	path("barcodes")

	output:
	path("checkedBarcodes")

	"""
	check_barcode_file.py barcodes > checkedBarcodes
	"""
}

/*
 * Splits reads according to the experimental barcode. Generates one read file for each experiment plus one for unmatched reads
 * Params: params.barcode_mismatches -> allowed number of mismatches allowed to still assign a read to an experiment
 * Input: tuple of [FASTQ] read file and [TSV] barcode file
 * Output: fastq_split_experimental_barcode -> [FASTQ] Several read files
 *	report_split_experimental_barcode -> [LOG] Report showing amount of reads assigned to each experiment
 */
process split_exp_barcode {
	tag {query.simpleName}
	publishDir "${params.output}/statistics/barcode-distribution", mode: 'copy', pattern: 'barcode_split.txt'

	input:
	tuple path(query), path(barcodes)

	output:
    path "output/*.fastq", emit: fastq_split_experimental_barcode
	path "barcode_split.log", emit: report_split_experimental_barcode 

	script:
	if(params.barcode_mismatches == 0)
		"""
		mkdir output
		cat ${query} | fastx_barcode_splitter.pl --bcfile ${barcodes} --bol --exact --prefix ./output/ --suffix .fastq > barcode_split.txt
		"""
	else
		"""
		mkdir output
		cat ${query} | fastx_barcode_splitter.pl --bcfile ${barcodes} --bol --mismatches ${params.barcode_mismatches} --prefix ./output/ --suffix .fastq > barcode_split.log
		"""
}

/*
 * Extracts length of experimental barcode
 * Input: [STRING] Barcode pattern
 * Output: [INTEGER] Length of experimental barcode + 1
 */
process get_length_exp_barcode {

	input:
	val(pattern)

	output:
	env(exp_barcode_length)

	"""
	exp_barcode_length=\$((\$(echo -n ${pattern} | sed 's/[^Xx]//g' | wc -m) + 1)) 
	"""
}

/*
 * Removes experimental barcode from 5' end of reads
 * Input: Tuple of [FASTQ] 1 read file and [INTEGER] length of experimental barcode + 1
 * Output: [FASTQ] Read file with experimental barcode removed
 */
process remove_exp_barcode {
	tag {query.simpleName}

	input:
	tuple path(query), val(length_exp_barcode)

	output:
	path "*.preprocessed.fastq"

	"""
	fastx_trimmer -f ${length_exp_barcode} -i ${query} -o ${query.baseName}.preprocessed.fastq
	"""
}

/*
 * Merges several FASTQ files belonging to the same experiment into a single file
 * Input: Tuple of [STRING] experiment name and [FASTQ] several read files
 * Output: [FASTQ] 1 read file
 */
process merge_preprocessed_reads {
	tag {name}

	input:
	tuple val(name), path("query")

	output:
	path("${name}.fastq")

	"""
	cat ${query} > ${name}.fastq
	"""
}

/*
 * Generates a plot showing the ampunt of reads being assigned to their experiments
 * Input: [LOG] reports generated by the process split_exp_barcode
 * Output: [PNG] Plot showing experimental barcode distribution
 * Pushes output to ${params.output}/statistics
 */
process generate_barcode_barplot {
	publishDir "${params.output}/statistics/barcode-distribution", mode: 'copy'

	input:
	path(query)

	output:
	path("${query.baseName}.png")

	"""
	plot_experimental_barcode_distribution.R --logs ${query} --output ${query.baseName} --type png --color "${params.color_barplot}"
	"""
}
