process quality_control {
	tag {query.simpleName}
	
	input:
	path(query)

	output:
	path("${query.baseName}*"), 			emit: summary
	path("${task.process}.version.txt"), 	emit: version

	"""
	fastqc ${query} -o .

	echo -e "${task.process}\tFastQC\t\$(fastqc --version | rev | cut -f 1 -d' ' | rev)" > ${task.process}.version.txt
	"""
}

process quality_control_2 {
	tag {query.simpleName}
	
	input:
	path(query)

	output:
	path("quality-control-2*"), 			emit: summary
	path("${task.process}.version.txt"), 	emit: version

	"""
	cat ${query} > quality-control-2.fastq
	fastqc quality-control-2.fastq -o .

	echo -e "${task.process}\tFastQC\t\$(fastqc --version | rev | cut -f 1 -d' ' | rev)" > ${task.process}.version.txt
	"""
}

process adapter_removal {
	tag {query.simpleName}

	input:
	path(query)

	output:
	path("${query}_trimmed.fq"), 			emit: fastq_trimmed
	path("${query}_trimming_report.txt"), 	emit: report_trimming
	path("${task.process}.version.txt"), 	emit: version

	"""
	trim_galore \
		--cores ${task.cpus} \
		--basename ${query} \
		-o . \
		--length ${params.min_length} \
		${query} \
		--quality 0

	echo -e "${task.process}\ttrim_galore\t\$(trim_galore -v | head -4 | tail -1 | sed -e 's/^[ \t]*//' | rev | cut -f 1 -d' ' | rev)" > ${task.process}.version.txt
	echo -e "${task.process}\tcutadapt\t\$(cutadapt --version)" >> ${task.process}.version.txt
	"""
}

process quality_filter {
	tag {query.simpleName}
	publishDir "${params.output}/statistics", mode: 'copy', pattern: "summary-quality-filter.txt"

	input:
	path(query)

	output:
	path("${query.baseName}.qual-filter.fastq"), 	emit: fastq_quality_filtered 
	path('summary-quality-filter.txt'), 			emit: report_quality_filter
	path("${task.process}.version.txt"), 			emit: version


	"""
	fastq_quality_filter \
		-v \
		-q ${params.min_qual} \
		-p ${params.min_percent_qual_filter} \
		-i ${query} \
		-o ${query.baseName}.qual-filter.fastq \
		> summary-quality-filter.txt

		echo -e "${task.process}\tfastq_quality_filter\t\$(fastq_quality_filter -h | head -2 | tail -1 | rev | cut -d ' ' -f 5- | rev)" > ${task.process}.version.txt
	"""
}