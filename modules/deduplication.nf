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

process deduplicate{
	publishDir "${params.output}/statistics/PCR-deduplication", mode: 'copy', pattern: "${query.baseName}.deduplicated.log*"
	tag {query.simpleName}
	memory { 20.GB + 1.B * query.size() }

	input:
	tuple path(query), path(index)

	output:
	path("${query.baseName}.deduplicated.bam"), emit: bam_deduplicated
	path("${query.baseName}.deduplicated.log*"), emit: report_deduplicated

	"""
	umi_tools dedup --random-seed=42 -I ${query} --output-stats ${query.baseName}.deduplicated.log -S ${query.baseName}.deduplicated.bam
	"""
}

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