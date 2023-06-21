/*
 * Generates BAM files from WIG files. For peak noted an entry is made in the BAM file
 * Input: [WIG] Peak file 
 * Output: [BAM] Alignment file depicting peak locations 
 */
process wig_to_bam {
	tag {query.simpleName}

	input:
	path(query)

	output:
	path("${query.baseName}.bam")

	"""
	wig-to-bam.py --input ${query} --output ${query.baseName}.bam
	"""
}

/*
 * Counts amount of times a specific feature is present. An entry is generated for every entry present in the BAM file (-R CORE)
 * Input: Tuple of [BAM] Peaks depicted by alignment file, [STR] one RNA subtype to be checked and [GFF3]|[GTF] an annotation file
 * Params: params.gene_id -> Name used for the gene_id int he annotation file
 * Output: [TSV] List of if entries belong to the current RNA subtype or not  
 */
process feature_counts {
	tag {query.simpleName}

	input:
	tuple path(query), val(rna_subtypes), path(annotation)

	output:
	path("${query.simpleName}.${rna_subtypes}.tsv")


	"""
	featureCounts -T ${task.cpus} -t ${rna_subtypes} -g ${params.gene_id} -a ${annotation} -R CORE -M -o ${query.simpleName} ${query}
	mv ${query}.featureCounts ${query.simpleName}.${rna_subtypes}.tsv
	"""
}

/*
 * Determines the RNA subtype distribution. Currently, if a peak got assigned to 2 subtypes it's declared ambiguous
 * Input: [STR] All RNA subtypes used for the analysis
 *			Tuple of [STR] Sample name and [TSV] several lists describing assignment to RNA subtypes 
 * Output: 	tsv_subtype_distribution 	-> [TSV] Distribution of peaks to all RNA subtypes  
 *			report_errors 				-> [LOG] Report for RNA subtype distribution
 */
process get_RNA_subtypes_distribution {
	tag {name}

	publishDir "${params.output}/RNA_subtypes", mode: 'copy'

	input:
	val(subtypes)
	tuple val(name), path(query)

	output:
	path("${name}.subtypes_distribution.tsv"), emit: tsv_subtype_distribution
	path("${name}.subtype.log"), emit: report_errors, optinoal: true

	script:
	subtypes_as_string = subtypes.join(' ')
	"""
	calc-RNA-subtypes-distribution.py --input ${query} --rna_subtypes ${subtypes_as_string} --output ${name}.subtypes_distribution.tsv > ${name}.subtype.log
	"""
}

/*
 * Generates a bar chart depicting the RNA subtype distribution calculated in get_RNA_subtypes_distribution
 * Input: [TSV] File containing the RNA subtype distribution
 * Params: params.color_barplot -> Hexcode for color used by the plot
 * Output: [PNG] Bar chart depicting RNA subtype distribtion  
 */
process generate_RNA_subtypes_barplot {
	publishDir "${params.output}/RNA_subtypes", mode: 'copy'

	input:
	path(query)

	output:
	path("${query.baseName}.png")

	"""
	RNA_subtypes_barcharts.R --input ${query} --output ${query.baseName} --type png --color "${params.color_barplot}"
	"""
}

/*
 * Collects warning from RNA subtype distribution and summarizes them in a single file
 * Input: [LOG] Warning from RNA subtype distribtion 
 * Output: [TXT] Collected warnings from RNA subtype distribution 
 */
process collect_subtype_analysis_errors {
	publishDir "${params.output}/statistics", mode: 'copy', pattern: 'subtype-analysis-warnings.txt'

	input:
	path(query)

	output:
	path("subtype-analysis-warnings.txt"), optional: true

	"""
	for i in ${query}
	do
		if [[ ! \$(cat \$i | wc -l) == 0 ]]; then
			echo "File: \$i" >> subtype-analysis-warnings.txt
			cat \$i >> subtype-analysis-warnings.txt
		fi
	done
	"""
}