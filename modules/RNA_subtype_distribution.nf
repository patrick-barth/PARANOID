/*
 * Generates BAM files from WIG files. For peak noted an entry is made in the BAM file
 */
process wig_to_bam {
	tag {query.simpleName}

	input:
	path(query)

	output:
	path("${query.baseName}.bam"), 			emit: bam
	path("${task.process}.version.txt"), 	emit: version

	"""
	wig-to-bam.py \
		--input ${query} \
		--output ${query.baseName}.bam

	echo -e "${task.process}\twig-to-bam.py\tcustom_script" > ${task.process}.version.txt
    echo -e "${task.process}\tpython\t\$(python --version | rev | cut -d' ' -f1 | rev)" >> ${task.process}.version.txt
	"""
}

/*
 * Counts amount of times a specific feature is present. An entry is generated for every entry present in the BAM file (-R CORE)
 */
process feature_counts {
	tag {query.simpleName}

	input:
	tuple path(query), val(rna_subtypes), path(annotation)

	output:
	path("${query.simpleName}.${rna_subtypes}.tsv"), 	emit: features
	path("${task.process}.version.txt"), 				emit: version

	"""
	featureCounts \
		-T ${task.cpus} \
		-t ${rna_subtypes} \
		-g ${params.gene_id} \
		-a ${annotation} \
		-R CORE \
		-M \
		-o ${query.simpleName} \
		${query}

	mv ${query}.featureCounts ${query.simpleName}.${rna_subtypes}.tsv

	echo -e "${task.process}\tfeatureCounts\t\$(featureCounts -v 2>&1 | head -2 | tail -1 |cut -d' ' -f2)" > ${task.process}.version.txt
	"""
}

/*
 * Determines the RNA subtype distribution. Currently, if a peak got assigned to 2 subtypes it's declared ambiguous
 */
process get_RNA_subtypes_distribution {
	tag {name}

	publishDir "${params.output}/RNA_subtypes", mode: 'copy'

	input:
	val(subtypes)
	tuple val(name), path(query)

	output:
	path("${name}.subtype_distribution.tsv"), 	emit: tsv_subtype_distribution
	path("${name}.subtype.log"), 				emit: report_errors, optional: true
	path("${name}.ambiguous.tsv"), 				emit: tsv_ambiguous_peaks, optional: true
	path("${task.process}.version.txt"), 		emit: version

	script:
	subtypes_as_string = subtypes.join(' ')
	"""
	calc-RNA-subtypes-distribution.py \
		--input ${query} \
		--rna_subtypes ${subtypes_as_string} \
		--output ${name} \
		--report_ambiguous \
		> ${name}.subtype.log

	echo -e "${task.process}\tcalc-RNA-subtypes-distribution.py\tcustom_script" > ${task.process}.version.txt
    echo -e "${task.process}\tpython\t\$(python --version | rev | cut -d' ' -f1 | rev)" >> ${task.process}.version.txt
	"""
}

/*
 * Generates a bar chart depicting the RNA subtype distribution calculated in get_RNA_subtypes_distribution  
 */
process generate_RNA_subtypes_barplot {
	publishDir "${params.output}/RNA_subtypes", mode: 'copy'

	input:
	path(query)

	output:
	path("${query.baseName}.png"), 			emit: png_to_output_dir
	path("${task.process}.version.txt"), 	emit: version

	"""
	RNA_subtypes_barcharts.R \
		--input ${query} \
		--output ${query.baseName} \
		--type png \
		--color "${params.color_barplot}"

	echo -e "${task.process}\tRNA_subtypes_barcharts.R\tcustom_script" > ${task.process}.version.txt
    echo -e "${task.process}\tR\t\$(R --version | head -1 | cut -d' ' -f3)" >> ${task.process}.version.txt
	"""
}

/*
 * Collects warning from RNA subtype distribution and summarizes them in a single file
 */
process collect_subtype_analysis_errors {
	publishDir "${params.output}/statistics", mode: 'copy', pattern: 'subtype-analysis-warnings.txt'

	input:
	path(query)

	output:
	path("subtype-analysis-warnings.txt"), optional: true, emit: report

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