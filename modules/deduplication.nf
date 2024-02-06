/*
 * Sorts and indexes BAM files 
 */
process sort_bam{
	tag {query.simpleName}

	input:
	path(query)

	output:
	tuple path("${query.baseName}.sorted.bam"), path("${query.baseName}.sorted.bam.bai"), emit: index_and_alignment

	"""
	samtools sort ${query} -o ${query.baseName}.sorted.bam
	samtools index ${query.baseName}.sorted.bam
	"""
}

/*
 * Removes PCR-duplicates by comparing the random barcode (which was saved in the header of each read) and the 
 *  alignment position. If both are identical for 2 or more reads all but one are removed.
 */
process deduplicate{
	publishDir "${params.output}/statistics/PCR-deduplication", mode: 'copy', pattern: "${query.baseName}.deduplicated.log*"
	tag {query.simpleName}
	memory { 40.GB + 5.B * query.size() }

	input:
	tuple path(query), path(index)

	output:
	path("${query.baseName}.deduplicated.bam"), emit: bam_deduplicated
	path("${query.baseName}.deduplicated.log*"), emit: report_deduplicated

	"""
	umi_tools dedup \
		--random-seed=42 \
		-I ${query} \
		--output-stats ${query.baseName}.deduplicated.log \
		-S ${query.baseName}.deduplicated.bam
	"""
}

/*
 * Merges several BAM files into a single one
 */
process merge_deduplicated_bam {
	tag {name}

	input:
	tuple val(name), path(query)

	output:
	path("${name}.bam"), emit: bam_merged

	"""
	samtools merge \
		-c \
		-p ${name}.bam \
		${query}
	"""
}
