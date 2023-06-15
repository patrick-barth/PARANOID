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

process get_RNA_subtypes_distribution {
	tag {name}

	publishDir "${params.output}/RNA_subtypes", mode: 'copy'

	input:
	val(subtypes)
	tuple val(name), path(query)

	output:
	path("${name}.subtypes_distribution.tsv"), emit: tsv_subtype_distribution
	path("${name}.subtype.log"), emit: report_errors, optinoal:true

	script:
	subtypes_as_string = subtypes.join(' ')
	"""
	calc-RNA-subtypes-distribution.py --input ${query} --rna_subtypes ${subtypes_as_string} --output ${name}.subtypes_distribution.tsv > ${name}.subtype.log
	"""
}

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