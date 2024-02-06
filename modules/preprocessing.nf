process quality_control {
	tag {query.simpleName}
	
	input:
	path(query)

	output:
	path("${query.baseName}*"), emit: summary

	"""
	fastqc ${query} -o .
	"""
}

process quality_control_2 {
	tag {query.simpleName}
	
	input:
	path(query)

	output:
	path("quality-control-2*"), emit: summary

	"""
	cat ${query} > quality-control-2.fastq
	fastqc quality-control-2.fastq -o .
	"""
}

process adapter_removal {
	tag {query.simpleName}

	input:
	path(query)

	output:
	path("${query}_trimmed.fq"), emit: fastq_trimmed
	path("${query}_trimming_report.txt"), emit: report_trimming

	"""
	trim_galore \
		--cores ${task.cpus} \
		--basename ${query} \
		-o . \
		--length ${params.min_length} \
		${query} \
		--quality 0
	"""
}

process quality_filter {
	tag {query.simpleName}
	publishDir "${params.output}/statistics", mode: 'copy', pattern: "summary-quality-filter.txt"

	input:
	path(query)

	output:
	path("${query.baseName}.qual-filter.fastq"), emit: fastq_quality_filtered 
	path('summary-quality-filter.txt'), emit: report_quality_filter 


	"""
	fastq_quality_filter \
		-v \
		-q ${params.min_qual} \
		-p ${params.min_percent_qual_filter} \
		-i ${query} \
		-o ${query.baseName}.qual-filter.fastq \
		> summary-quality-filter.txt
	"""
}