/*
 * Sorts and indexes BAM files
 * Input: [BAM] Aligned sequences 
 * Output: Tuple of [BAM] Sorted aligned sequences and [BAI] Index file  
 */
process sort_bam{
	tag {query.simpleName}

	input:
	path(query)

	output:
	tuple path("${query.baseName}.sorted.bam"), path("${query.baseName}.sorted.bam.bai")

	"""
	samtools sort ${query} -o ${query.baseName}.sorted.bam
	samtools index ${query.baseName}.sorted.bam
	"""
}

/*
 * Removes PCR-duplicates by comparing the random barcode (which was saved in the header of each read) and the 
 *  alignment position. If both are identical for 2 or more reads all but one are removed.
 * Input: Tuple of [BAM] Sorted aligned sequences and [BAI] Index file
 * Output: bam_deduplicated -> [BAM] Aligned sequences without PCR duplicates 
 *		report_deduplicated -> [LOG] Report of deduplication
 */
process deduplicate{
	publishDir "${params.output}/statistics/PCR-deduplication", mode: 'copy', pattern: "${query.baseName}.deduplicated.log*"
	tag {query.simpleName}
	memory { 40.GB + 5.B * query.size() }
	queue 'chaos'

	input:
	tuple path(query), path(index)

	output:
	path("${query.baseName}.deduplicated.bam"), emit: bam_deduplicated
	path("${query.baseName}.deduplicated.log*"), emit: report_deduplicated

	"""
	umi_tools dedup --random-seed=42 -I ${query} --output-stats ${query.baseName}.deduplicated.log -S ${query.baseName}.deduplicated.bam
	"""
}

/*
 * Merges several BAM files into a single one
 * Input: Tuple of [STR] Name of outptu file and [BAM] several BAM files  
 * Output: [BAM] Aligned sequences
 */
process merge_deduplicated_bam {
	tag {name}

	input:
	tuple val(name), path(query)

	output:
	path("${name}.bam")

	"""
	samtools merge -c -p ${name}.bam ${query}
	"""
}
